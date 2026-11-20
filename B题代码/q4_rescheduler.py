# -*- coding: utf-8 -*-
"""Q4: 扰动重调度与策略比较 (P0 保持原方案 / P1 事件触发重优化 / P2 规则式贪心).

运行: python -X utf8 -B q4_rescheduler.py
输出: q4_output/ 下 q4_summary.json 与图件.

口径(题面 Q4-01..Q4-18 + 台账, 框架 §4.5):
- 基线 = Q2(2) 最优方案; 扰动 t=5 s; 评估窗 [0,60] s;
- 扰动前观测已开始的警报视为已完成、不可重排(Q4-08); 已派出但未开始观测的车可改派,
  重排自当前位置/姿态继续(Q4-08 台账#10);
- eta = k/N: N = 窗内应观测警报总数(含新警报), k = 成功启动联合观测的警报数
  (两车到站且完成指向); t_last = 最晚观测开始时刻, 有 pending 警报无法在窗内启动则
  t_last = 60; 字典序: 先 eta 最大, 再 t_last 最小;
- 通信中断 [5, t_rec]: 中断期按原指令继续, 恢复后方接新指令 -> 决策时点 = t_rec(Q4-05);
- 无先知: P0/P1 在 t=0 均不知扰动; P1 仅在扰动实际发生后触发;
- 新警报 T4 = O + R(cosθ, sinθ), H=50 (BG20 极坐标口径);
- 失效车退场不设障(Q4-04); 与 Q2(2) 一致的受限搜索: 每车至多承担 2 个任务槽,
  行驶 = shortest_path 绕障最短路/3 (有限图上界), 同警报两站方位角差 >= 1°.

策略:
- P0 保持原方案(解析求值: 锁定警报计入, 失效车未开始的警报其配对车若含失效车则失败);
- P1 决策时点确定性重优化: 结构枚举(警报-车辆对指派 + 复用车序) -> 贪心粗筛 ->
  best-response 精化, 目标 (k 最大, t_last 最小);
- P2 规则式: 每轮选"最早可完成"的警报, 在当前空闲车对中按联合 ready 最小指派(可解释基线).
"""
import json
import math
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations, permutations
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geometry
from geometry import SPEED, TARGETS, VEHICLE_INIT, Target, shortest_path
from simulator import (GIMBAL_RATE, OBS_DURATION, Task, _topo_order, pointing,
                       rotate_duration, simulate)
from q2_2_scheduler import AZ_MIN, alert_gamma, build_alert_data, travel_scalar, wrap180

OUT = Path(__file__).with_name("q4_output")
OUT.mkdir(exist_ok=True)
BASE_JSON = Path(__file__).with_name("q2_output") / "q2_2_summary.json"
WINDOW = 60.0
BAND_PROXY = 120      # 每槽先按欧氏距离取前 120 候选做精确行程
BAND_TOP = 80         # 静态 ready 带宽
SCREEN_TOP = 60       # 结构粗筛保留数
REFINE_TOP = 6        # 精化结构数
REFINE_CAND = 20      # 精化每槽候选数
REFINE_ROUNDS = 1
_ORIG_PAIRS = {"T1": (3, 5), "T2": (2, 4), "T3": (3, 1)}   # Q2(2) 基线指派

# ---------------------------------------------------------------- 通用工具
_TRAVEL_CACHE = {}


def _cached_travel(pos, station):
    key = (round(pos[0], 3), round(pos[1], 3),
           round(station[0], 3), round(station[1], 3))
    if key not in _TRAVEL_CACHE:
        _TRAVEL_CACHE[key] = travel_scalar(pos, station)
    return _TRAVEL_CACHE[key]


def point_along(path, dist):
    """折线上自起点行进 dist 米后的位置."""
    remain = dist
    for (x0, y0), (x1, y1) in zip(path[:-1], path[1:]):
        seg = math.hypot(x1 - x0, y1 - y0)
        if remain <= seg + 1e-12:
            f = 0.0 if seg == 0 else remain / seg
            return (x0 + f * (x1 - x0), y0 + f * (y1 - y0))
        remain -= seg
    return (path[-1][0], path[-1][1])


# ---------------------------------------------------------------- 基线重建与 τ 状态
def rebuild_base():
    """重演 Q2(2) 最优方案, 返回 (res, plan)."""
    d = json.load(open(BASE_JSON, encoding="utf-8"))
    best = d["best"]
    stn = {tuple(k.split("@")): tuple(v) for k, v in best["stations"].items()}
    stn = {(int(v[1:]), a): s for (v, a), s in stn.items()}
    struct = best["structure"]
    k = struct["double_veh"]
    a1, a2 = struct["ordered_alerts"]
    a3 = [a for a in ("T1", "T2", "T3") if a not in (a1, a2)][0]
    seqs = {k: [(k, a1), (k, a2)],
            struct["a1_partner"]: [(struct["a1_partner"], a1)],
            struct["a2_partner"]: [(struct["a2_partner"], a2)],
            struct["a3_vehicles"][0]: [(struct["a3_vehicles"][0], a3)],
            struct["a3_vehicles"][1]: [(struct["a3_vehicles"][1], a3)]}
    tasks = [Task(v, a, stn[(v, a)]) for v, pairs in seqs.items() for (v, a) in pairs]
    res = simulate(tasks)
    assert abs(res.max_start - best["max_start_s"]) < 1e-3
    return res, stn


def capture_states(res_base, plan, tau, failed=None):
    """各车在 tau 时刻的状态: 位置/云台姿态/最早可接受新指令时刻.

    在途: 沿最短路折线插值, 云台锁定正东; 转位中: 姿态按比例插值;
    等待/观测中: 姿态为指向角, 观测中的车 free_t = 观测毕 + 复位时长.
    """
    states = {}
    for v in range(1, 6):
        if v == failed:
            continue
        pos_prev = VEHICLE_INIT[v]
        st = None
        for e in res_base.events[v]:
            if e["kind"] == "depart":
                continue
            t0, t1 = e["t_start"], e["t_end"]
            station = plan[(v, e["alert"])]
            if tau >= t1 - 1e-9:
                pos_prev = station
                continue
            pan, tilt = pointing(station, TARGETS[e["alert"]])
            if e["kind"] == "arrive":          # 在途
                length, path = shortest_path(pos_prev, station)
                st = {"pos": point_along(path, (tau - t0) * SPEED),
                      "pan": 0.0, "tilt": 0.0, "free_t": tau}
            elif e["kind"] == "slew":          # 转位中(可中断重指)
                f = (tau - t0) / (t1 - t0)
                st = {"pos": station, "pan": wrap180(f * wrap180(pan)),
                      "tilt": f * tilt, "free_t": tau}
            elif e["kind"] == "wait":          # 已指向等待联合同步
                st = {"pos": station, "pan": pan, "tilt": tilt, "free_t": tau}
            elif e["kind"] == "obs":           # 锁定观测中, 完成观测+复位后释放
                st = {"pos": station, "pan": pan, "tilt": tilt,
                      "free_t": t1 + rotate_duration(pan, tilt, 0.0, 0.0)}
            elif e["kind"] == "reset":         # 复位中(可中断重指)
                f = (tau - t0) / (t1 - t0)
                st = {"pos": station, "pan": wrap180((1 - f) * wrap180(pan)),
                      "tilt": (1 - f) * tilt, "free_t": tau}
            break
        if st is None:                          # tau 晚于全部事件
            st = {"pos": pos_prev, "pan": 0.0, "tilt": 0.0, "free_t": tau}
        states[v] = st
    locked = {a: {"start": s, "end": res_base.alert_end[a]}
              for a, s in res_base.alert_start.items() if s <= tau + 1e-9}
    return states, locked


# ---------------------------------------------------------------- 续仿真
def simulate_resume(tasks, states, tau):
    """异构初态续仿真: 与 simulator.simulate 同语义, 但各车初值取自 states.

    首任务转位自当前姿态; 之后(复位毕)自正东水平. 返回 dict 形式结果.
    """
    veh_tasks, alert_tasks = {}, {}
    for t in tasks:
        veh_tasks.setdefault(t.vehicle, []).append(t)
        alert_tasks.setdefault(t.alert, []).append(t)
    for a, ts in alert_tasks.items():
        vs = [t.vehicle for t in ts]
        if len(ts) != 2 or vs[0] == vs[1]:
            raise ValueError(f"警报 {a} 指派非法: {vs}")
    edges = set()
    for v, ts in veh_tasks.items():
        for i in range(1, len(ts)):
            edges.add((ts[i - 1].alert, ts[i].alert))
    order = _topo_order(alert_tasks.keys(), edges)

    st = {v: {"pos": tuple(map(float, states[v]["pos"])),
              "t": max(float(states[v]["free_t"]), float(tau)),
              "pan": float(states[v]["pan"]), "tilt": float(states[v]["tilt"]),
              "first": True, "events": []} for v in veh_tasks}
    alert_start, alert_end = {}, {}
    for a in order:
        tgt = TARGETS[a]
        plans = []
        for t in alert_tasks[a]:
            s = st[t.vehicle]
            tv = _cached_travel(s["pos"], t.station)
            if tv is None:
                return None                      # 不可达 -> 该结构不可行
            arrive = s["t"] + tv
            pan, tilt = pointing(t.station, tgt)
            base = (s["pan"], s["tilt"]) if s["first"] else (0.0, 0.0)
            slew = rotate_duration(base[0], base[1], pan, tilt)
            plans.append({"v": t.vehicle, "depart": s["t"], "arrive": arrive,
                          "slew": slew, "ready": arrive + slew, "pan": pan,
                          "tilt": tilt, "station": t.station})
        start = max(p["ready"] for p in plans)
        end = start + OBS_DURATION
        alert_start[a], alert_end[a] = start, end
        for p in plans:
            s = st[p["v"]]
            ev = s["events"]
            ev.append({"kind": "depart", "alert": a, "t_start": p["depart"],
                       "t_end": p["depart"], "detail": ""})
            ev.append({"kind": "arrive", "alert": a, "t_start": p["depart"],
                       "t_end": p["arrive"], "detail": ""})
            ev.append({"kind": "slew", "alert": a, "t_start": p["arrive"],
                       "t_end": p["ready"], "detail": ""})
            if p["ready"] < start - 1e-12:
                ev.append({"kind": "wait", "alert": a, "t_start": p["ready"],
                           "t_end": start, "detail": ""})
            ev.append({"kind": "obs", "alert": a, "t_start": start,
                       "t_end": end, "detail": ""})
            reset = rotate_duration(p["pan"], p["tilt"], 0.0, 0.0)
            ev.append({"kind": "reset", "alert": a, "t_start": end,
                       "t_end": end + reset, "detail": ""})
            s["t"] = end + reset
            s["pos"] = p["station"]
            s["pan"] = s["tilt"] = 0.0
            s["first"] = False
    return {"events": {v: s["events"] for v, s in st.items()},
            "alert_start": alert_start, "alert_end": alert_end,
            "makespan": max(s["t"] for s in st.values())}


# ---------------------------------------------------------------- 结构枚举
def _finish_spec(assign):
    """按指派生成车辆任务序的全部分合法组合(每车至多 2 任务, 序不产生循环)."""
    per_v = {}
    for a, vs in assign.items():
        for v in vs:
            per_v.setdefault(v, []).append(a)
    doublers = {v: sorted(al) for v, al in per_v.items() if len(al) == 2}
    order_lists = [{}]
    for v, als in doublers.items():
        new = []
        for base in order_lists:
            for seq in permutations(als):
                d = dict(base)
                d[v] = list(seq)
                new.append(d)
        order_lists = new
    out = []
    for orders in order_lists:
        edges = set()
        for v, seq in orders.items():
            edges.add((seq[0], seq[1]))
        try:
            proc = _topo_order(list(assign.keys()), edges)
        except ValueError:
            continue                                 # 序矛盾
        out.append({"assign": dict(assign), "orders": orders, "proc": proc})
    return out


def gen_specs(A, V):
    """警报-车辆指派结构枚举(每车 ≤2 槽, 每警报恰 2 辆不同车)."""
    specs = []
    if len(A) == 1:
        for p in combinations(V, 2):
            specs += _finish_spec({A[0]: list(p)})
    elif len(A) == 2:
        a1, a2 = A
        for p1 in combinations(V, 2):
            for p2 in combinations(V, 2):
                cnt = Counter(p1) + Counter(p2)
                if any(c > 2 for c in cnt.values()):
                    continue
                specs += _finish_spec({a1: list(p1), a2: list(p2)})
    else:                                            # 3 警报
        a1, a2, a3 = A
        if len(V) == 5:                              # 恰一复用车 (Q2(2) 360 族)
            for k in V:
                others = [v for v in V if v != k]
                for x, y in permutations(A, 2):
                    rem = [a for a in A if a not in (x, y)][0]
                    for p1 in others:
                        for p2 in others:
                            if p2 == p1:
                                continue
                            rest = [v for v in others if v not in (p1, p2)]
                            specs += _finish_spec({x: [k, p1], y: [k, p2],
                                                   rem: list(rest)})
        else:                                        # 4 车 3 警报: 双复用
            for p_a in combinations(V, 2):
                for p_b in combinations(V, 2):
                    for p_c in combinations(V, 2):
                        cnt = Counter(p_a) + Counter(p_b) + Counter(p_c)
                        if any(c > 2 for c in cnt.values()):
                            continue
                        specs += _finish_spec({a1: list(p_a), a2: list(p_b),
                                               a3: list(p_c)})
    # 去重: 同一 (assign, orders) 视为同结构
    seen, uniq = set(), []
    for s in specs:
        key = (tuple(sorted((a, tuple(vs)) for a, vs in s["assign"].items())),
               tuple(sorted((v, tuple(seq)) for v, seq in s["orders"].items())))
        if key not in seen:
            seen.add(key)
            uniq.append(s)
    return uniq


# ---------------------------------------------------------------- 候选带与贪心
def build_bands(states, V, A, AD, tau):
    """bands[(v, a)] = [(ready_est, 候选索引)] 升序, 精确绕障行程."""
    bands = {}
    for a in A:
        ad = AD[a]
        arr = ad["arr"]
        for v in V:
            s = states[v]
            t0 = max(s["free_t"], tau)
            d = np.hypot(arr[:, 0] - s["pos"][0], arr[:, 1] - s["pos"][1])
            idxs = np.argsort(d)[:BAND_PROXY]
            lst = []
            for i in idxs:
                tv = _cached_travel(s["pos"], arr[i])
                if tv is None:
                    continue
                if s["free_t"] > tau + 1e-9:         # 完成锁定任务后自正东出发
                    sle = ad["slew"][i]
                else:                                # 决策时点就地改派
                    sle = rotate_duration(s["pan"], s["tilt"],
                                          ad["pan"][i], ad["tilt"][i])
                lst.append((t0 + tv + sle, int(i)))
            lst.sort()
            bands[(v, a)] = lst[:BAND_TOP]
    return bands


def _joint_pick(l1, l2, ad, k=40):
    """两排序候选表中最小 max(ready) 且方位角差达标的组合."""
    arr = ad["arr"]
    best = None
    for r1, i1 in l1[:k]:
        s1 = arr[i1]
        for r2, i2 in l2[:k]:
            if best is not None and r2 > best[0]:
                break
            if alert_gamma(ad, s1, arr[i2]) < AZ_MIN:
                continue
            m = max(r1, r2)
            if best is None or m < best[0]:
                best = (m, i1, i2)
    return best


def greedy_spec_plan(spec, states, AD, bands, tau):
    plan = {}
    for a in spec["proc"]:
        v1, v2 = spec["assign"][a]
        ad = AD[a]
        lists = []
        for v in (v1, v2):
            prev = [aa for aa, vs in spec["assign"].items()
                    if aa != a and v in vs and (v, aa) in plan]
            if prev:                                  # 复用车第二段: 自第一站出发
                s1 = plan[(v, prev[0])]
                d = np.hypot(ad["arr"][:, 0] - s1[0], ad["arr"][:, 1] - s1[1])
                near = set(int(i) for i in np.argsort(d)[:BAND_PROXY])
                lst = []
                for r0, i in bands[(v, a)]:
                    if i in near:
                        tv = _cached_travel(s1, ad["arr"][i])
                        if tv is not None:
                            lst.append((tv + ad["slew"][i], i))
                if not lst:                           # 回退: 带内候选逐一精确试
                    for r0, i in bands[(v, a)][:40]:
                        tv = _cached_travel(s1, ad["arr"][i])
                        if tv is not None:
                            lst.append((tv + ad["slew"][i], i))
                lst.sort()
                if not lst:
                    return None
                lists.append(lst)
            else:
                lists.append(bands[(v, a)])
        pick = _joint_pick(lists[0], lists[1], ad)
        if pick is None:
            return None
        _, i1, i2 = pick
        plan[(v1, a)] = tuple(map(float, ad["arr"][i1]))
        plan[(v2, a)] = tuple(map(float, ad["arr"][i2]))
    return plan


def repair_plan(spec, plan, states, AD, bands, tau):
    """方位角差不足时, 用带宽内候选替换某一槽(按 ready 代理择优)."""
    for a in spec["proc"]:
        ad = AD[a]
        v1, v2 = spec["assign"][a]
        s1, s2 = plan[(v1, a)], plan[(v2, a)]
        if alert_gamma(ad, s1, s2) >= AZ_MIN:
            continue
        best = None
        for v, other, mine in ((v1, s2, s1), (v2, s1, s2)):
            for r0, i in bands[(v, a)][:40]:
                s = tuple(map(float, ad["arr"][i]))
                if alert_gamma(ad, s, other) < AZ_MIN:
                    continue
                key = max(r0, 0.0)
                if best is None or key < best[0]:
                    best = (key, v, s)
        if best is not None:
            plan[(best[1], a)] = best[2]
    return plan


def spec_tasks(spec, plan):
    tasks = []
    seq_of = {v: [] for v in range(1, 6)}
    for v, seq in spec["orders"].items():
        seq_of[v] = list(seq)
    singles = []
    for a, vs in spec["assign"].items():
        for v in vs:
            if a not in seq_of[v]:
                singles.append((v, a))
    for v, seq in seq_of.items():
        for a in seq:
            tasks.append(Task(v, a, plan[(v, a)]))
    for v, a in singles:
        tasks.append(Task(v, a, plan[(v, a)]))
    return tasks


def evaluate_spec(spec, plan, states, tau, locked, N):
    tasks = spec_tasks(spec, plan)
    res = simulate_resume(tasks, states, tau)
    if res is None:
        return None
    k = sum(1 for s in res["alert_start"].values() if s <= WINDOW)
    full = k + len(locked) == N                    # N 含锁定与不可覆盖警报
    starts = {a: s for a, s in res["alert_start"].items()}
    all_starts = list(starts.values()) + [lk["start"] for lk in locked.values()]
    t_last = max(all_starts) if full else WINDOW
    return {"k": k + len(locked), "eta": (k + len(locked)) / N,
            "t_last": t_last, "starts": starts, "res": res,
            "full": full, "plan": plan}


def better(cur, cand):
    if cand is None:
        return cur
    if cur is None:
        return cand
    key_c, key_n = (cur["k"], -cur["t_last"]), (cand["k"], -cand["t_last"])
    return cand if key_n > key_c else cur


def refine_spec(spec, plan, states, AD, tau, locked, N, cur, bands):
    plan = {kk: tuple(map(float, vv)) for kk, vv in plan.items()}
    best_plan, best_eval = plan, cur
    for _ in range(REFINE_ROUNDS):
        improved = False
        for (v, a) in list(plan.keys()):
            ad = AD[a]
            partner_v = [u for u in spec["assign"][a] if u != v][0]
            arr = ad["arr"]
            d = np.hypot(arr[:, 0] - plan[(v, a)][0],
                         arr[:, 1] - plan[(v, a)][1])
            near = [int(i) for i in np.argsort(d)[:REFINE_CAND]]
            band_idx = [i for _, i in bands.get((v, a), [])[:REFINE_CAND]]
            for i in dict.fromkeys(near + band_idx):
                s_new = tuple(map(float, arr[i]))
                if s_new == plan[(v, a)]:
                    continue
                if alert_gamma(ad, s_new, plan[(partner_v, a)]) < AZ_MIN:
                    continue
                trial = dict(plan)
                trial[(v, a)] = s_new
                ev = evaluate_spec(spec, trial, states, tau, locked, N)
                if ev is not None and better(best_eval, ev) is ev:
                    best_plan, best_eval = trial, ev
                    improved = True
        if not improved:
            break
        plan = best_plan
    return best_plan, best_eval


# ---------------------------------------------------------------- 三策略
def p0_eval(locked, failed, has_new, base_starts):
    N = 3 + (1 if has_new else 0)
    k = len(locked)
    started = [lk["start"] for lk in locked.values()]
    for a, vs in _ORIG_PAIRS.items():
        if a in locked:
            continue
        if failed is None or failed not in vs:
            k += 1
            started.append(base_starts[a])
    t_last = max(started) if k == N else WINDOW
    return {"k": k, "eta": k / N, "t_last": t_last,
            "note": "保持原方案: 新警报不响应; 失效车未开始的警报若其配对含失效车则失败"}


def _proc_from(orders, assign):
    edges = set()
    for v, seq in orders.items():
        edges.add((seq[0], seq[1]))
    try:
        return _topo_order(list(assign.keys()), edges)
    except ValueError:
        return list(assign.keys())


def _dyn_list(v, a, s, ad, bands, tau):
    arr = ad["arr"]
    t0 = max(s["free_t"], tau)
    d = np.hypot(arr[:, 0] - s["pos"][0], arr[:, 1] - s["pos"][1])
    out = []
    for r0, i in bands.get((v, a), [])[:30]:
        tv = _cached_travel(s["pos"], arr[i])
        if tv is None:
            continue
        sle = (ad["slew"][i] if s["free_t"] > tau + 1e-9
               else rotate_duration(s["pan"], s["tilt"], ad["pan"][i], ad["tilt"][i]))
        out.append((t0 + tv + sle, int(i)))
    out.sort()
    return out


# ---------------------------------------------------------------- 场景求解
def solve_scenario(sc):
    """单场景: 返回 P0/P1/P2 的 (eta, t_last) 与细节. 进程池 worker."""
    try:
        return _solve_scenario_impl(sc)
    except Exception as exc:                         # 单场景失败不拖垮整批
        import traceback
        return {"id": sc["id"], "type": sc["type"], "params": {
                    k2: v for k2, v in sc.items() if k2 not in ("id", "type")},
                "error": f"{exc.__class__.__name__}: {exc}",
                "traceback": traceback.format_exc()}


def _solve_scenario_impl(sc):
    has_new = "a" in sc["type"]
    N = 3 + (1 if has_new else 0)
    failed = sc.get("failed")
    tau = sc.get("t_rec", 5.0) if "c" in sc["type"] else 5.0
    if has_new:
        geometry.TARGETS["T4"] = Target("T4", float(sc["theta"]), float(sc["R"]), 50.0)

    res_base, plan_base = rebuild_base()
    states, locked = capture_states(res_base, plan_base, tau, failed)
    A = [a for a in ("T1", "T2", "T3") if a not in locked]
    if has_new:
        A.append("T4")
    V = list(states.keys())

    out = {"id": sc["id"], "type": sc["type"], "params": {
        k2: v for k2, v in sc.items() if k2 not in ("id", "type")}}
    out["tau_decision"] = tau
    out["locked"] = sorted(locked)
    out["open_alerts"] = A
    out["vehicles"] = V
    out["N"] = N

    p0 = p0_eval(locked, failed, has_new, res_base.alert_start)
    out["P0"] = p0

    if not A:                                        # 无可重排对象
        for p in ("P1", "P2"):
            out[p] = dict(p0)
            out[p]["note"] = "决策时点前全部警报已启动, 无可重排对象"
        out["states"] = _states_ser(states)
        return out

    AD = {a: build_alert_data(a) for a in A}
    for a in A:
        if AD[a]["n"] == 0:
            out.setdefault("structural", []).append(f"{a} 可观测域为空")
    AD = {a: ad for a, ad in AD.items() if ad["n"] > 0}
    A = [a for a in A if a in AD]
    if not A:
        for p in ("P1", "P2"):
            out[p] = dict(p0)
            out[p]["note"] = "待观测警报均不可观测"
        out["states"] = _states_ser(states)
        return out

    bands = build_bands(states, V, A, AD, tau)

    # P1: 结构枚举 -> 粗筛 -> 精化
    specs = gen_specs(A, V)

    def proxy(sp):
        s = 0.0
        for a, vs in sp["assign"].items():
            r = [bands[(v, a)][0][0] for v in vs if bands.get((v, a))]
            s += max(r) if r else 1e9
        return s

    specs.sort(key=proxy)
    best = None
    evaluated = []
    for sp in specs[:SCREEN_TOP]:
        pl = greedy_spec_plan(sp, states, AD, bands, tau)
        if pl is None:
            continue
        pl = repair_plan(sp, pl, states, AD, bands, tau)
        ev = evaluate_spec(sp, pl, states, tau, locked, N)
        if ev is None:
            continue
        evaluated.append((sp, ev))
        best = better(best, ev)
    evaluated.sort(key=lambda x: (x[1]["k"], -x[1]["t_last"]))
    for sp, ev in evaluated[:REFINE_TOP]:
        pl2, ev2 = refine_spec(sp, dict(ev["plan"]), states, AD, tau, locked, N,
                               ev, bands)
        if ev2 is not None:
            evaluated.append((sp, ev2))
            best = better(best, ev2)
    if best is None:
        out["P1"] = dict(p0)
        out["P1"]["note"] = "未找到任何可行重排(结构均不可行)"
    else:
        sp_best = next(sp for sp, ev in evaluated if ev is best)
        out["P1"] = {"k": best["k"], "eta": best["eta"], "t_last": best["t_last"],
                     "starts": {a: round(s, 3) for a, s in best["starts"].items()},
                     "assign": {a: list(vs) for a, vs in sp_best["assign"].items()},
                     "orders": {v: list(q) for v, q in sp_best["orders"].items()},
                     "stations": {f"V{v}@{a}": [round(s[0], 3), round(s[1], 3)]
                                  for (v, a), s in best["plan"].items()}}

    # P2: 规则式(记录轮次以定复用车序)
    vst = {v: dict(states[v]) for v in V}
    remaining = list(A)
    assign, stations, round_of = {}, {}, {}
    ridx = 0
    while remaining:
        bestp = None
        for a in remaining:
            ad = AD[a]
            for v1, v2 in combinations(V, 2):
                l1 = _dyn_list(v1, a, vst[v1], ad, bands, tau)
                l2 = _dyn_list(v2, a, vst[v2], ad, bands, tau)
                pick = _joint_pick(l1, l2, ad, k=15)
                if pick is None:
                    continue
                m, i1, i2 = pick
                if bestp is None or m < bestp[0]:
                    bestp = (m, a, v1, v2, i1, i2)
        if bestp is None:
            break
        m, a, v1, v2, i1, i2 = bestp
        ridx += 1
        round_of[a] = ridx
        assign[a] = (v1, v2)
        stations[(v1, a)] = tuple(map(float, AD[a]["arr"][i1]))
        stations[(v2, a)] = tuple(map(float, AD[a]["arr"][i2]))
        for v, i in ((v1, i1), (v2, i2)):
            s = tuple(map(float, AD[a]["arr"][i]))
            pan, tilt = pointing(s, TARGETS[a])
            end = m + OBS_DURATION
            vst[v] = {"pos": s, "pan": 0.0, "tilt": 0.0,
                      "free_t": end + rotate_duration(pan, tilt, 0.0, 0.0)}
        remaining.remove(a)
    if not assign:
        out["P2"] = dict(p0)
        out["P2"]["note"] = "规则式无可行指派"
    else:
        orders = {}
        for v in V:
            als = [a for a in assign if v in assign[a]]
            if len(als) == 2:
                orders[v] = sorted(als, key=lambda a: round_of[a])
        spec2 = {"assign": assign, "orders": orders,
                 "proc": _proc_from(orders, assign)}
        tasks2 = []
        for v in V:
            seq = orders.get(v, [])
            tasks2 += [Task(v, a, stations[(v, a)]) for a in seq]
            tasks2 += [Task(v, a, stations[(v, a)]) for a in assign
                       if v in assign[a] and a not in seq]
        res2 = simulate_resume(tasks2, states, tau)
        if res2 is None:
            out["P2"] = dict(p0)
            out["P2"]["note"] = "规则式方案不可行"
        else:
            k2 = sum(1 for s in res2["alert_start"].values() if s <= WINDOW)
            allst = list(res2["alert_start"].values()) + \
                [lk["start"] for lk in locked.values()]
            out["P2"] = {"k": k2 + len(locked), "eta": (k2 + len(locked)) / N,
                         "t_last": max(allst) if k2 + len(locked) == N else WINDOW,
                         "starts": {a: round(s, 3) for a, s in res2["alert_start"].items()},
                         "assign": {a: list(vs) for a, vs in assign.items()}}
    out["states"] = _states_ser(states)
    return out


def _states_ser(states):
    return {f"V{v}": {"pos": [round(s["pos"][0], 3), round(s["pos"][1], 3)],
                      "pan_deg": round(s["pan"], 2), "tilt_deg": round(s["tilt"], 2),
                      "free_t": round(s["free_t"], 3)} for v, s in states.items()}


# ---------------------------------------------------------------- 场景集
SECTORS = [15 + 30 * i for i in range(12)]


def build_scenarios():
    scs = []
    nid = 0
    for th in SECTORS:                                # (a) 36
        for R in (30.0, 40.0, 50.0):
            nid += 1
            scs.append({"id": f"a-{nid:02d}", "type": "a", "theta": th, "R": R})
    for f in range(1, 6):                             # (b) 5
        nid += 1
        scs.append({"id": f"b-{nid:02d}", "type": "b", "failed": f})
    for tr in (10.0, 20.0, 30.0):                     # (c) 3
        nid += 1
        scs.append({"id": f"c-{nid:02d}", "type": "c", "t_rec": tr})
    for th in SECTORS:                                # (a+c) 24
        for tr in (10.0, 30.0):
            nid += 1
            scs.append({"id": f"a+c-{nid:02d}", "type": "a+c",
                        "theta": th, "R": 40.0, "t_rec": tr})
    for f in range(1, 6):                             # (b+c) 5
        nid += 1
        scs.append({"id": f"b+c-{nid:02d}", "type": "b+c", "failed": f, "t_rec": 10.0})
    hard = {f: hard_sector(f) for f in range(1, 6)}   # (a+b) 5
    for f in range(1, 6):
        nid += 1
        scs.append({"id": f"a+b-{nid:02d}", "type": "a+b", "failed": f,
                    "theta": hard[f], "R": 40.0})
    for f in (2, 3):                                  # (a+b+c) 2
        nid += 1
        scs.append({"id": f"a+b+c-{nid:02d}", "type": "a+b+c", "failed": f,
                    "theta": hard[f], "R": 40.0, "t_rec": 10.0})
    for tr in (7.5, 12.5, 15.0, 25.0):                # (a+c) t_rec 扫线加密 4
        nid += 1
        scs.append({"id": f"a+c-sweep-{nid:02d}", "type": "a+c",
                    "theta": 15.0, "R": 40.0, "t_rec": tr})
    for tr in (15.0, 25.0):                           # (c) 扫线加密 2
        nid += 1
        scs.append({"id": f"c-sweep-{nid:02d}", "type": "c", "t_rec": tr})
    return scs


def hard_sector(failed):
    """对失效场景逐扇区探测新警报最早可完成时刻, 取最难扇区."""
    res_base, plan_base = rebuild_base()
    states, locked = capture_states(res_base, plan_base, 5.0, failed)
    V = list(states.keys())
    worst, worst_val = None, -1.0
    for th in SECTORS:
        geometry.TARGETS["T4"] = Target("T4", float(th), 40.0, 50.0)
        ad = build_alert_data("T4")
        if ad["n"] == 0:
            val = 1e9
        else:
            arr = ad["arr"]
            top = {}
            for v in V:
                s = states[v]
                d = np.hypot(arr[:, 0] - s["pos"][0], arr[:, 1] - s["pos"][1])
                cand = []
                for i in np.argsort(d)[:20]:
                    tv = travel_scalar(s["pos"], arr[i])
                    if tv is None:
                        continue
                    sle = (ad["slew"][i] if s["free_t"] > 5.0 + 1e-9 else
                           rotate_duration(s["pan"], s["tilt"],
                                           ad["pan"][i], ad["tilt"][i]))
                    cand.append((max(s["free_t"], 5.0) + tv + sle, int(i)))
                cand.sort()
                top[v] = cand[:3]
            vals = []
            for v1, v2 in combinations(V, 2):
                pick = _joint_pick(top[v1], top[v2], ad, k=3)
                if pick is not None:
                    vals.append(pick[0])
            val = min(vals) if vals else 1e9
        if val > worst_val:
            worst, worst_val = th, val
    return worst


# ---------------------------------------------------------------- 主程序
def main():
    scs = build_scenarios()
    print(f"== Q4 场景总数 {len(scs)} ==")
    results = []
    workers = min(6, os.cpu_count() or 1)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, r in enumerate(ex.map(solve_scenario, scs)):
            results.append(r)
            if (i + 1) % 10 == 0:
                print(f"  {i + 1}/{len(scs)} 完成")
    by_id = {r["id"]: r for r in results}

    # 汇总
    summary = {"window_s": WINDOW, "metrics_rule": "eta=k/N 先最大, t_last 后最小; "
               "有 pending 警报未启动则 t_last=60", "scenarios": results}
    with open(OUT / "q4_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    _fig_a_grid(by_id)
    _fig_b(by_id)
    _fig_trec(by_id)
    _fig_value(results)
    _fig_gantt_examples(by_id)
    _print_digest(results)
    print(f"输出已写入 {OUT}")


def _print_digest(results):
    print("\n== 策略对比摘要 ==")
    for r in results:
        if "error" in r:
            print(f"[{r['id']:>12}] ERROR: {r['error']}")
            continue
        if r["type"] == "a":
            continue
        p0, p1, p2 = r["P0"], r["P1"], r["P2"]
        print(f"[{r['id']:>12}] P0 η={p0['eta']:.3f} t={p0['t_last']:.2f} | "
              f"P1 η={p1['eta']:.3f} t={p1['t_last']:.2f} | "
              f"P2 η={p2['eta']:.3f} t={p2['t_last']:.2f}")


def _fig_a_grid(by_id):
    a_ids = [k for k in by_id if k.startswith("a-") and "+c" not in k]
    grid = np.full((3, 12), np.nan)
    grid_eta = np.zeros((3, 12))
    for k in a_ids:
        r = by_id[k]
        th = float(r["params"]["theta"]); R = float(r["params"]["R"])
        i = {30.0: 0, 40.0: 1, 50.0: 2}[R]
        j = SECTORS.index(th)
        grid[i, j] = r["P1"]["t_last"]
        grid_eta[i, j] = r["P1"]["eta"]
    fig, ax = plt.subplots(figsize=(11.5, 4.2))
    im = ax.imshow(grid, cmap="RdYlGn_r", aspect="auto", vmin=10, vmax=60)
    ax.set_xticks(range(12)); ax.set_xticklabels([f"{int(s)}°" for s in SECTORS])
    ax.set_yticks(range(3)); ax.set_yticklabels(["R=30", "R=40", "R=50"])
    for i in range(3):
        for j in range(12):
            if np.isnan(grid[i, j]):
                ax.text(j, i, "n/a", ha="center", fontsize=7)
            else:
                txt = f"{grid[i, j]:.1f}" if grid_eta[i, j] == 1 else \
                    f"η={grid_eta[i, j]:.2f}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=7)
    ax.set_xlabel("new alert azimuth θ (sector center)")
    ax.set_title("Q4 (a): P1 t_last (s) over new-alert grid; η<1 cells shown as η; "
                 "P0 constant: η=0.75, t_last=60", fontsize=10)
    fig.colorbar(im, ax=ax, label="t_last (s)")
    fig.tight_layout()
    fig.savefig(OUT / "q4_fig_a_grid.png", dpi=300)
    plt.close(fig)


def _fig_b(by_id):
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    labels, p0v, p1v, p2v = [], [], [], []
    for k in sorted(k for k in by_id if k.startswith("b-") and "+" not in k):
        r = by_id[k]
        labels.append(f"V{r['params']['failed']} fails")
        p0v.append(r["P0"]["t_last"]); p1v.append(r["P1"]["t_last"])
        p2v.append(r["P2"]["t_last"])
    x = np.arange(len(labels))
    ax.bar(x - 0.25, p0v, 0.25, label="P0 retain", color="#bdc3c7")
    ax.bar(x, p1v, 0.25, label="P1 re-opt", color="#e67e22")
    ax.bar(x + 0.25, p2v, 0.25, label="P2 rule", color="#5dade2")
    for i, k in enumerate(sorted(k for k in by_id if k.startswith("b-") and "+" not in k)):
        r = by_id[k]
        for dx, p in ((-0.25, r["P0"]), (0, r["P1"]), (0.25, r["P2"])):
            ax.text(x[i] + dx, p["t_last"] + 1,
                    f"η={p['eta']:.2f}", ha="center", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("t_last (s)  (60 = some alert never starts)")
    ax.set_ylim(0, 70)
    ax.legend(fontsize=8)
    ax.set_title("Q4 (b) single-vehicle failure: t_last by strategy "
                 "(T1 locked → always counted completed)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "q4_fig_b_bars.png", dpi=300)
    plt.close(fig)


def _fig_trec(by_id):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    # 左: 纯 (c) 扫线
    cks = sorted((k for k, r in by_id.items() if r["type"] == "c"),
                 key=lambda k: float(by_id[k]["params"]["t_rec"]))
    axes[0].plot([float(by_id[k]["params"]["t_rec"]) for k in cks],
                 [by_id[k]["P1"]["t_last"] for k in cks],
                 "o-", color="#e67e22", label="P1")
    axes[0].plot([float(by_id[k]["params"]["t_rec"]) for k in cks],
                 [by_id[k]["P2"]["t_last"] for k in cks],
                 "^-", color="#5dade2", label="P2")
    axes[0].axhline(11.046, color="#7f8c8d", ls="--", lw=1, label="P0 = 11.046 s")
    axes[0].set_title("Q4 (c) pure outage: t_rec sweep")
    axes[0].set_xlabel("t_rec (s)"); axes[0].set_ylabel("t_last (s)")
    axes[0].set_ylim(8, 65)
    axes[0].legend(fontsize=8)
    # 右: (a+c) θ=15°, R=40 扫线
    acks = sorted((k for k, r in by_id.items() if r["type"] == "a+c"
                   and float(r["params"]["theta"]) == 15.0
                   and float(r["params"]["R"]) == 40.0),
                  key=lambda k: float(by_id[k]["params"]["t_rec"]))
    axes[1].plot([float(by_id[k]["params"]["t_rec"]) for k in acks],
                 [by_id[k]["P1"]["t_last"] for k in acks],
                 "o-", color="#1f618d", label="P1")
    axes[1].plot([float(by_id[k]["params"]["t_rec"]) for k in acks],
                 [by_id[k]["P0"]["t_last"] for k in acks],
                 "s--", color="#bdc3c7", label="P0")
    axes[1].set_title("Q4 (a+c): recovery delay vs performance (θ=15°, R=40)")
    axes[1].set_xlabel("t_rec (s)"); axes[1].set_ylabel("t_last (s)")
    axes[1].set_ylim(8, 65)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "q4_fig_trec.png", dpi=300)
    plt.close(fig)


def _fig_value(results):
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    colors = {"a": "#e67e22", "b": "#1f618d", "c": "#7f8c8d",
              "a+c": "#8e44ad", "b+c": "#16a085", "a+b": "#c0392b",
              "a+b+c": "#2c3e50"}
    for r in results:
        if "error" in r:
            continue
        de = r["P1"]["eta"] - r["P0"]["eta"]
        dt = r["P1"]["t_last"] - r["P0"]["t_last"]
        ax.scatter(dt, de, color=colors.get(r["type"], "k"), s=28, alpha=0.8)
    ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("Δt_last = t_last(P1) − t_last(P0)  (s)")
    ax.set_ylabel("Δη = η(P1) − η(P0)")
    ax.set_title("Q4 rescheduling value boundary (P1 vs P0, one point per scenario)",
                 fontsize=10)
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=t)
               for t, c in colors.items()]
    ax.legend(handles=handles, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "q4_fig_value.png", dpi=300)
    plt.close(fig)


def _fig_gantt_examples(by_id):
    KIND_COLOR = {"arrive": "#5dade2", "slew": "#f5b041", "wait": "#bdc3c7",
                  "obs": "#58d68d", "reset": "#af7ac5"}
    picks = []
    for k in by_id:
        if k.startswith("a-") and by_id[k]["params"]["theta"] == 195.0 \
                and float(by_id[k]["params"]["R"]) == 50.0:
            picks.append((k, "P1"))
    for k in by_id:
        if k.startswith("b-") and by_id[k]["params"]["failed"] == 3 and "+" not in k:
            picks.append((k, "P1"))
    fig, axes = plt.subplots(len(picks), 1,
                             figsize=(11.5, 4.0 * max(len(picks), 1)))
    if len(picks) == 1:
        axes = [axes]
    for ax, (k, pol) in zip(axes, picks):
        r = by_id[k]
        if "error" in r or "assign" not in r.get("P1", {}):
            continue
        # 重演该策略时间轴
        tau = r["tau_decision"]
        if "T4" in r.get("open_alerts", []):
            geometry.TARGETS["T4"] = Target("T4", float(r["params"]["theta"]),
                                            float(r["params"]["R"]), 50.0)
        res_base, plan_base = rebuild_base()
        states, locked = capture_states(res_base, plan_base, tau,
                                        r["params"].get("failed"))
        assign = {a: tuple(vs) for a, vs in r["P1"]["assign"].items()}
        orders = {int(v): list(q) for v, q in r["P1"]["orders"].items()}
        stations = {}
        for key, s in r["P1"]["stations"].items():
            v, a = key.split("@")
            stations[(int(v[1:]), a)] = tuple(s)
        singles = [(v, a) for a, vs in assign.items() for v in vs
                   if a not in orders.get(v, [])]
        tasks = [Task(v, a, stations[(v, a)]) for v, q in orders.items()
                 for a in q] + [Task(v, a, stations[(v, a)]) for v, a in singles]
        res = simulate_resume(tasks, states, tau)
        if res is None:
            continue
        title = f"{k} — P1  η={r['P1']['eta']:.2f}, t_last={r['P1']['t_last']:.2f} s"
        vids = sorted(res["events"])
        locked_names = set(r["locked"])
        for row, v in enumerate(vids):
            # 扰动前实际执行段: 锁定警报(含 T1 观测尾段/复位)全程绘制,
            # 其余(被改派中断的原任务)裁剪到决策时点 τ
            for e in res_base.events.get(v, []):
                if e["kind"] == "depart":
                    continue
                full = e["alert"] in locked_names
                t1 = e["t_end"] if full else min(e["t_end"], tau)
                t0 = e["t_start"]
                if t1 <= t0 + 1e-9:
                    continue
                ax.barh(row, t1 - t0, left=t0, height=0.55,
                        color=KIND_COLOR.get(e["kind"], "#ccc"),
                        edgecolor="k", linewidth=0.4, alpha=0.55)
                if e["kind"] == "obs" and t1 - t0 > 1.2:
                    ax.text((t0 + t1) / 2, row, e["alert"], ha="center",
                            va="center", fontsize=7, fontweight="bold", alpha=0.8)
            for e in res["events"][v]:
                if e["kind"] == "depart":
                    continue
                ax.barh(row, e["t_end"] - e["t_start"], left=e["t_start"],
                        height=0.55, color=KIND_COLOR.get(e["kind"], "#ccc"),
                        edgecolor="k", linewidth=0.4)
                if e["kind"] == "obs":
                    ax.text((e["t_start"] + e["t_end"]) / 2, row, e["alert"],
                            ha="center", va="center", fontsize=8, fontweight="bold")
        ax.axvline(tau, color="k", ls=":", lw=1.2)
        ax.text(tau + 0.1, len(vids) - 0.6, f"τ_d={tau:.0f}s", fontsize=8)
        ax.set_yticks(range(len(vids)))
        ax.set_yticklabels([f"V{v}" for v in vids])
        ax.set_title(title, fontsize=10)
        ax.set_ylim(-0.8, len(vids) - 0.3)
    axes[-1].set_xlabel("time (s)")
    fig.tight_layout()
    fig.savefig(OUT / "q4_fig_gantt_examples.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
