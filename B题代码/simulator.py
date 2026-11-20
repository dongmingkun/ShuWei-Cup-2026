# -*- coding: utf-8 -*-
"""B题事件驱动时间轴仿真器: 行驶/转位/5s联合观测/复位 状态机.

规则口径(GC/TS 全套):
- 恒速 3 m/s, 行程时间 = 绕障最短路长 / 3(转向与避让不耗时但路径必须绕行);
- 行进期云台锁定正东水平; 停稳后才转位; 转位/复位时长 = max(|Δpan|,|Δtilt|)/30°/s;
- 同一警报的两车都完成指向才启动 5s 联合观测 => 警报开始时刻 = 两车 ready 时刻的 max;
- 观测毕立即复位到正东水平(已在则复位时长 0), 复位毕立即驶向下一任务.

声明的建模约定(论文假设章须写明):
- 云台 pan 轴按最短角差旋转(卷绕到 (-180,180]); tilt 无卷绕;
- 车辆初始/复位后云台姿态为 (pan=0 正东, tilt=0 水平), 出发时刻为 t0.
"""
import math
from dataclasses import dataclass

from geometry import (CAMERA_HEIGHT, SPEED, TARGETS, VEHICLE_INIT,
                      shortest_path, is_station_feasible)

GIMBAL_RATE = 30.0    # °/s, 两轴并行
OBS_DURATION = 5.0    # s, 联合观测时长


def wrap180(a: float) -> float:
    """卷绕到 (-180, 180]."""
    a = (a + 180.0) % 360.0 - 180.0
    return 180.0 if a == -180.0 else a


def pointing(station, target):
    """站到目标的指向角 (pan°自正东逆时针∈[0,360), tilt°俯仰)."""
    dx = target.x - station[0]
    dy = target.y - station[1]
    pan = math.degrees(math.atan2(dy, dx)) % 360.0
    tilt = math.degrees(math.atan2(target.H - CAMERA_HEIGHT, math.hypot(dx, dy)))
    return pan, tilt


def rotate_duration(cur_pan: float, cur_tilt: float,
                    dst_pan: float, dst_tilt: float) -> float:
    """云台转位/复位时长 = max(|Δpan|(最短角差), |Δtilt|) / 30°/s."""
    d_pan = abs(wrap180(dst_pan - cur_pan))
    d_tilt = abs(dst_tilt - cur_tilt)
    return max(d_pan, d_tilt) / GIMBAL_RATE


@dataclass
class Task:
    vehicle: int
    alert: str
    station: tuple        # (x, y) 站点(=相机盘心水平坐标)


@dataclass
class SimulationResult:
    events: dict                      # vehicle -> [event dict] (按时间序)
    alert_start: dict                 # alert -> 联合观测开始时刻
    alert_end: dict
    max_start: float                  # Q2(2) 目标值 min max(start)
    makespan: float                   # 最后一辆车复位完成时刻

    def timeline_str(self, vehicle: int) -> str:
        lines = [f"== V{vehicle} 时间轴 =="]
        for e in self.events[vehicle]:
            lines.append(f"[{e['t_start']:8.3f} -> {e['t_end']:8.3f}] "
                         f"{e['kind']:<7} {e['alert'] or '':<4} {e['detail']}")
        return "\n".join(lines)


def _topo_order(alerts, edges):
    """Kahn 拓扑排序; edges = {(prev_alert, next_alert)}; 有环则抛错."""
    indeg = {a: 0 for a in alerts}
    succ = {a: [] for a in alerts}
    for u, v in edges:
        succ[u].append(v)
        indeg[v] += 1
    queue = sorted(a for a in alerts if indeg[a] == 0)
    order = []
    while queue:
        a = queue.pop(0)
        order.append(a)
        for w in succ[a]:
            indeg[w] -= 1
            if indeg[w] == 0:
                queue.append(w)
        queue.sort()
    if len(order) != len(alerts):
        raise ValueError(f"任务顺序存在循环依赖, 无法调度: {sorted(alerts)} / 边 {sorted(edges)}")
    return order


def simulate(tasks, t0: float = 0.0, init_positions=None, enforce: bool = True):
    """按各车任务队列执行事件仿真.

    tasks: list[Task], 每车任务按列表顺序执行; 每个警报必须恰好两辆不同车各 1 个 Task;
    t0: 出发时刻(问题2(2)取 0; 问题4 重调度取决策时点, 各车云台须已在正东水平);
    enforce=True 时站点违反三约束直接抛错.
    """
    if init_positions is None:
        init_positions = dict(VEHICLE_INIT)
    veh_tasks, alert_tasks = {}, {}
    for t in tasks:
        veh_tasks.setdefault(t.vehicle, []).append(t)
        alert_tasks.setdefault(t.alert, []).append(t)
    for a, ts in alert_tasks.items():
        vs = [t.vehicle for t in ts]
        if len(ts) != 2 or vs[0] == vs[1]:
            raise ValueError(f"警报 {a} 必须恰好由两辆不同车观测, 实际 {vs}")

    bad = [f"V{t.vehicle}@{t.alert} 站点 {t.station} 违反三约束"
           for t in tasks if not is_station_feasible(t.station, TARGETS[t.alert])]
    if bad:
        if enforce:
            raise ValueError("; ".join(bad))

    edges = set()
    for v, ts in veh_tasks.items():
        for i in range(1, len(ts)):
            edges.add((ts[i - 1].alert, ts[i].alert))
    order = _topo_order(alert_tasks.keys(), edges)

    state = {v: {"pos": tuple(map(float, init_positions[v])), "t": float(t0),
                 "events": []} for v in veh_tasks}
    alert_start, alert_end = {}, {}

    for a in order:
        tgt = TARGETS[a]
        plans = []
        for t in alert_tasks[a]:
            st = state[t.vehicle]
            depart = st["t"]
            res = shortest_path(st["pos"], t.station)
            if res is None:
                raise ValueError(f"V{t.vehicle} 前往 {t.alert} 站点 {t.station} 不可达")
            length, path = res
            arrive = depart + length / SPEED
            pan, tilt = pointing(t.station, tgt)
            slew = rotate_duration(0.0, 0.0, pan, tilt)
            plans.append({"v": t.vehicle, "depart": depart, "arrive": arrive,
                          "slew": slew, "ready": arrive + slew, "pan": pan,
                          "tilt": tilt, "length": length, "path": path,
                          "station": tuple(map(float, t.station))})
        start = max(p["ready"] for p in plans)
        end = start + OBS_DURATION
        alert_start[a], alert_end[a] = start, end
        for p in plans:
            st = state[p["v"]]
            ev = st["events"]
            ev.append({"kind": "depart", "alert": a, "t_start": p["depart"],
                       "t_end": p["depart"],
                       "detail": f"自{tuple(round(c, 3) for c in st['pos'])}出发, "
                                 f"路径长 {p['length']:.3f} m"})
            ev.append({"kind": "arrive", "alert": a, "t_start": p["depart"],
                       "t_end": p["arrive"],
                       "detail": f"到站 {tuple(round(c, 3) for c in p['station'])}"})
            ev.append({"kind": "slew", "alert": a, "t_start": p["arrive"],
                       "t_end": p["ready"],
                       "detail": f"pan={p['pan']:.2f}° tilt={p['tilt']:.2f}°, "
                                 f"转位 {p['slew']:.3f} s"})
            if p["ready"] < start - 1e-12:
                ev.append({"kind": "wait", "alert": a, "t_start": p["ready"],
                           "t_end": start, "detail": "等待联合观测同步"})
            ev.append({"kind": "obs", "alert": a, "t_start": start,
                       "t_end": end, "detail": "联合观测 5 s"})
            reset = rotate_duration(p["pan"], p["tilt"], 0.0, 0.0)
            ev.append({"kind": "reset", "alert": a, "t_start": end,
                       "t_end": end + reset, "detail": f"复位 {reset:.3f} s"})
            st["t"] = end + reset
            st["pos"] = p["station"]

    max_start = max(alert_start.values())
    makespan = max(st["t"] for st in state.values())
    return SimulationResult(events={v: st["events"] for v, st in state.items()},
                            alert_start=alert_start, alert_end=alert_end,
                            max_start=max_start, makespan=makespan)
