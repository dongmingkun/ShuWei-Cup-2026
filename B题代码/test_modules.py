# -*- coding: utf-8 -*-
"""几何引擎与仿真器单元测试: python test_modules.py (全部通过打印 ALL TESTS PASSED)."""
import math
import sys
from types import SimpleNamespace

from geometry import (BUILDINGS, TARGETS, VEHICLE_INIT, clearance_ok,
                      dist_point_to_rect, elev_min_radius, elevation_ok,
                      is_station_feasible, los_ok, segment_hits_box,
                      segment_rect_distance, shortest_path, travel_time)
from simulator import (OBS_DURATION, GIMBAL_RATE, Task, _topo_order,
                       rotate_duration, simulate)

TOL = 1e-6
PASS = 0


def check(name, cond):
    global PASS
    if not cond:
        print(f"FAIL: {name}")
        sys.exit(1)
    PASS += 1
    print(f"ok  {name}")


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol


# ---------- 1. 目标坐标换算(逆时针公式) ----------
t1, t2, t3 = TARGETS["T1"], TARGETS["T2"], TARGETS["T3"]
check("T1=(84.6410,70.0)", close(t1.x, 84.6410161, 1e-4) and close(t1.y, 70.0, 1e-9))
check("T2=(11.0289,72.5)", close(t2.x, 11.0288569, 1e-4) and close(t2.y, 72.5, 1e-9))
check("T3=(59.8351,13.2948)", close(t3.x, 59.8351237, 1e-4) and close(t3.y, 13.2948186, 1e-4))

# ---------- 2. 车辆几何约束(贴边允许) ----------
check("clearance (50,50) ok", clearance_ok((50, 50)))
check("clearance B1内部 (23,21) 拒", not clearance_ok((23, 21)))
check("clearance B1贴边 (20.5,21) 许", clearance_ok((20.5, 21)))
check("clearance 边界外 (1.4,50) 拒", not clearance_ok((1.4, 50)))
check("clearance 边界贴边 (1.5,50) 许", clearance_ok((1.5, 50)))
check("dist_point_to_rect", close(dist_point_to_rect(20.5, 21, 22, 34, 16, 26), 1.5))

# ---------- 3. 俯仰约束 ----------
r80 = (t1.H - 2) / math.tan(math.radians(80))
p_east = (t1.x + r80, t1.y)
check("elevation r=r80 恰80° 可行", elevation_ok(p_east, t1))
check("elevation r=r80-0.1 拒", not elevation_ok((t1.x + r80 - 0.1, t1.y), t1))
check("elev_min_radius(T3)≈9.3453", close(elev_min_radius(t3), 53 / math.tan(math.radians(80)), 1e-9))

# ---------- 4. 视线遮挡(切触即遮挡) ----------
B1 = BUILDINGS[0]
box1 = (B1.x0, B1.x1, B1.y0, B1.y1, 0.0, B1.h)
check("遮挡: 穿B1 (15,20)->(40,20,H5)", not los_ok((15, 20), SimpleNamespace(x=40, y=20, H=5)))
check("无遮挡: (15,20)->(40,50,H5)", los_ok((15, 20), SimpleNamespace(x=40, y=50, H=5)))
check("切触面 y=26 算遮挡", segment_hits_box((20, 26, 2), (40, 26, 2), box1))
check("面外 1e-6 不遮挡", not segment_hits_box((20, 26 + 1e-6, 2), (40, 26 + 1e-6, 2), box1))
check("楼顶上方 z=25 不遮挡", not segment_hits_box((20, 21, 25), (36, 21, 25), box1))
check("楼体内部高度 z=10 遮挡", segment_hits_box((20, 21, 10), (36, 21, 10), box1))
check("T1 站(92.22,70) 无遮挡", los_ok((t1.x + r80, t1.y), t1))

# ---------- 5. 绕障最短路 ----------
d_straight = math.dist((88, 68), (95, 95))
L, path = shortest_path((88, 68), (95, 95))
check("直连无障: 长度=直线距", close(L, d_straight, 1e-6) and len(path) == 2)
check("travel_time=(L/3)", close(travel_time((88, 68), (95, 95)), d_straight / 3, 1e-9))
L3, _ = shortest_path((15, 14), (50, 93))
check("(15,14)->(50,93) 擦B1角但真实净距1.539>=1.5 直连可行",
      close(L3, math.dist((15, 14), (50, 93)), 1e-6))
check("(15,14)->(50,93) 对B1净距恰>=1.5(贴边口径)",
      segment_rect_distance((15, 14), (50, 93), B1.x0, B1.x1, B1.y0, B1.y1)
      >= 1.5 - 1e-9)

straight2 = math.dist((15, 14), (38, 20))
L2, path2 = shortest_path((15, 14), (38, 20))
check("绕B1: 长度>直线", L2 > straight2 + 0.5)
check("绕B1: 长度合理上界", L2 < straight2 + 30)
check("绕B1: 各顶点满足净距", all(clearance_ok(p) for p in path2))
check("绕B1: 各段真实净距>=1.5", all(
    segment_rect_distance(path2[i], path2[i + 1],
                          b.x0, b.x1, b.y0, b.y1) >= 1.5 - 1e-9
    for i in range(len(path2) - 1) for b in BUILDINGS))

Lc, _ = shortest_path((20.8, 14.8), (20, 14))
check("B1西南角: 合法直连可达且=直线距",
      Lc is not None and close(Lc, math.dist((20.8, 14.8), (20, 14)), 1e-9))
Lf, pf = shortest_path((20.8, 14.8), (60, 60))
check("B1西南角膨胀带内合法点可接入路网远行",
      Lf is not None and all(clearance_ok(p) for p in pf)
      and all(segment_rect_distance(pf[i], pf[i + 1],
                                    b.x0, b.x1, b.y0, b.y1) >= 1.5 - 1e-9
              for i in range(len(pf) - 1) for b in BUILDINGS))

try:
    shortest_path((23, 21), (50, 50))
    check("起点在楼内应抛错", False)
except ValueError:
    check("起点在楼内应抛错", True)

# ---------- 6. 云台转位(最短角差卷绕) ----------
check("pan 0->180 = 6s", close(rotate_duration(0, 0, 180, 0), 180 / GIMBAL_RATE))
check("pan 0->270 最短=90 = 3s", close(rotate_duration(0, 0, 270, 0), 90 / GIMBAL_RATE))
check("pan 270->0 最短=90 = 3s", close(rotate_duration(270, 0, 0, 0), 90 / GIMBAL_RATE))
check("pan180+tilt80 取长轴=6s", close(rotate_duration(0, 0, 180, 80), 6.0))
check("复位归零 = max(pan,tilt)/30", close(rotate_duration(90, 80, 0, 0), 3.0))
check("已正东水平复位=0", close(rotate_duration(0, 0, 0, 0), 0.0))

# ---------- 7. 仿真器: T3 双车联合观测(手推数值) ----------
r3 = elev_min_radius(t3)
S_a = (t3.x, t3.y - r3)   # V1 从南侧站, 直连无障碍
S_b = (t3.x, t3.y + r3)   # V2 从北侧站, 直连无障碍
check("场景站点可行", is_station_feasible(S_a, t3) and is_station_feasible(S_b, t3))
check("V1->S_a 直连无障", close(shortest_path(VEHICLE_INIT[1], S_a)[0],
                                math.dist(VEHICLE_INIT[1], S_a), 1e-6))
check("V2->S_b 直连无障", close(shortest_path(VEHICLE_INIT[2], S_b)[0],
                                math.dist(VEHICLE_INIT[2], S_b), 1e-6))

res = simulate([Task(1, "T3", S_a), Task(2, "T3", S_b)])
tr1 = math.dist(VEHICLE_INIT[1], S_a) / 3
tr2 = math.dist(VEHICLE_INIT[2], S_b) / 3
slew = max(90.0, 80.0) / GIMBAL_RATE          # 两站指向 pan=±90°, tilt=80°
exp_ready1, exp_ready2 = tr1 + slew, tr2 + slew
exp_start = max(exp_ready1, exp_ready2)
check("V1 ready=行程+转位", close(res.events[1][2]["t_end"], exp_ready1, 1e-9))
check("V2 ready=行程+转位", close(res.events[2][2]["t_end"], exp_ready2, 1e-9))
check("联合观测开始=max(ready)", close(res.alert_start["T3"], exp_start, 1e-9))
check("观测结束=+5s", close(res.alert_end["T3"], exp_start + OBS_DURATION, 1e-9))
check("早到者有 wait 事件", any(e["kind"] == "wait" for e in res.events[2]))
check("复位=3s(pan90/270最短)", close(
    [e for e in res.events[1] if e["kind"] == "reset"][0]["t_end"]
    - [e for e in res.events[1] if e["kind"] == "reset"][0]["t_start"], 3.0, 1e-9))
check("makespan=观测结束+复位", close(res.makespan, exp_start + OBS_DURATION + 3.0, 1e-9))
check("max_start=警报开始", close(res.max_start, res.alert_start["T3"], 1e-9))

# ---------- 8. 拓扑排序与异常 ----------
check("拓扑: 链 T1->T2->T3", _topo_order({"T1", "T2", "T3"},
      {("T1", "T2"), ("T2", "T3")}) == ["T1", "T2", "T3"])
try:
    _topo_order({"T1", "T2"}, {("T1", "T2"), ("T2", "T1")})
    check("循环依赖应抛错", False)
except ValueError:
    check("循环依赖应抛错", True)

try:
    simulate([Task(1, "T3", (23, 21)), Task(2, "T3", S_b)])
    check("站点入楼应抛错", False)
except ValueError:
    check("站点入楼应抛错", True)

try:
    simulate([Task(1, "T3", S_a)])
    check("单车观测应抛错", False)
except ValueError:
    check("单车观测应抛错", True)

print(f"\nALL TESTS PASSED ({PASS} checks)")
