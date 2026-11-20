# -*- coding: utf-8 -*-
"""B题公共几何引擎: 三约束判定(车辆净距/视线遮挡/云台俯仰) + 绕障最短路.

口径(与 B题-题意逐句翻译.md 一致):
- 方位角按题面公式 (50+R*cosθ, 50+R*sinθ) 的标准数学逆时针约定;
- 视线段与建筑棱柱任一棱/面/顶点切触即遮挡(闭集判定);
- 车辆圆盘半径1.5m: 盘心距边界与建筑表面 >= 1.5m(允许贴边, 不入内部);
- 相机位于盘心、高2m; 云台俯仰可行域 [-10°, 80°](闭区间).
"""
import heapq
import math
from dataclasses import dataclass

FIELD_SIZE = 100.0
VEHICLE_RADIUS = 1.5
SPEED = 3.0                 # m/s
CAMERA_HEIGHT = 2.0         # m
ELEV_MIN = -10.0            # deg
ELEV_MAX = 80.0             # deg
DISPATCH = (50.0, 50.0)
EPS = 1e-9


@dataclass(frozen=True)
class Building:
    x0: float
    x1: float
    y0: float
    y1: float
    h: float


BUILDINGS = [
    Building(22, 34, 16, 26, 20),   # B1
    Building(63, 79, 59, 71, 28),   # B2
    Building(10, 24, 65, 76, 12),   # B3
    Building(72, 86, 22, 34, 24),   # B4
    Building(41, 51, 38, 46, 32),   # B5
]

VEHICLE_INIT = {1: (15.0, 14.0), 2: (88.0, 10.0), 3: (14.0, 86.0),
                4: (88.0, 68.0), 5: (50.0, 93.0)}


@dataclass(frozen=True)
class Target:
    name: str
    theta_deg: float
    R: float
    H: float

    @property
    def x(self) -> float:
        return 50.0 + self.R * math.cos(math.radians(self.theta_deg))

    @property
    def y(self) -> float:
        return 50.0 + self.R * math.sin(math.radians(self.theta_deg))


TARGETS = {
    "T1": Target("T1", 30.0, 40.0, 45.0),
    "T2": Target("T2", 150.0, 45.0, 35.0),
    "T3": Target("T3", 285.0, 38.0, 55.0),
}


def elev_min_radius(target) -> float:
    """俯仰约束给出的最小水平距离 r >= (H-2)/tan(80°)."""
    return (target.H - CAMERA_HEIGHT) / math.tan(math.radians(ELEV_MAX))


# ---------------------------------------------------------------- 约束1: 车辆几何

def dist_point_to_rect(px: float, py: float, x0: float, x1: float,
                       y0: float, y1: float) -> float:
    dx = max(x0 - px, 0.0, px - x1)
    dy = max(y0 - py, 0.0, py - y1)
    return math.hypot(dx, dy)


def clearance_ok(p, margin: float = VEHICLE_RADIUS) -> bool:
    """盘心 p 处圆盘不越场界、不入建筑内部(距离恰为 margin 的贴边允许)."""
    x, y = p
    if not (margin - EPS <= x <= FIELD_SIZE - margin + EPS
            and margin - EPS <= y <= FIELD_SIZE - margin + EPS):
        return False
    for b in BUILDINGS:
        if dist_point_to_rect(x, y, b.x0, b.x1, b.y0, b.y1) < margin - EPS:
            return False
    return True


# ---------------------------------------------------------------- 约束2: 视线遮挡

def segment_hits_box(p0, p1, box, eps: float = EPS) -> bool:
    """3D 线段 p0->p1 与闭长方体 box=(xmin,xmax,ymin,ymax,zmin,zmax) 相交(含切触)则 True.

    slab 法; eps 用于把数值上贴面的情形判为接触. 遮挡口径要求"切触即遮挡",
    故默认 eps>0(略宽松地判遮挡); 求最短路时传 eps=0 并先收缩矩形(允许贴边绕行).
    """
    xmin, xmax, ymin, ymax, zmin, zmax = box
    d = (p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2])
    tmin, tmax = 0.0, 1.0
    for a, (lo, hi) in enumerate(((xmin, xmax), (ymin, ymax), (zmin, zmax))):
        if abs(d[a]) < 1e-15:
            if p0[a] < lo - eps or p0[a] > hi + eps:
                return False
        else:
            t1 = (lo - p0[a]) / d[a]
            t2 = (hi - p0[a]) / d[a]
            if t1 > t2:
                t1, t2 = t2, t1
            tmin = max(tmin, t1)
            tmax = min(tmax, t2)
            if tmin > tmax + eps:
                return False
    return True


def los_ok(p, target) -> bool:
    """相机(盘心, 高2m)到目标的视线不被任一建筑遮挡(切触算遮挡)."""
    p0 = (p[0], p[1], CAMERA_HEIGHT)
    p1 = (target.x, target.y, target.H)
    for b in BUILDINGS:
        if segment_hits_box(p0, p1, (b.x0, b.x1, b.y0, b.y1, 0.0, b.h)):
            return False
    return True


# ---------------------------------------------------------------- 约束3: 云台俯仰

def elevation_angle(p, target) -> float:
    r = math.hypot(target.x - p[0], target.y - p[1])
    return math.degrees(math.atan2(target.H - CAMERA_HEIGHT, r))


def elevation_ok(p, target) -> bool:
    e = elevation_angle(p, target)
    return ELEV_MIN - EPS <= e <= ELEV_MAX + EPS


# ---------------------------------------------------------------- 三约束合成

def is_station_feasible(p, target) -> bool:
    """三约束同时满足(可观测域判定核心谓词)."""
    return (clearance_ok(p)
            and elevation_ok(p, target)
            and los_ok(p, target))


# ---------------------------------------------------------------- 绕障最短路

def inflate_rect(x0, x1, y0, y1, m: float = VEHICLE_RADIUS):
    """建筑外扩 m 得盘心禁区(锐角版; 贴边允许)."""
    return (x0 - m, x1 + m, y0 - m, y1 + m)


def _seg_crosses_rect_interior(p, q, rect) -> bool:
    """2D 线段与(略收缩的)矩形内域相交则 True; 贴边/擦角允许."""
    xmin, xmax, ymin, ymax = rect
    shrink = 1e-9
    return segment_hits_box(
        (p[0], p[1], 0.0), (q[0], q[1], 0.0),
        (xmin + shrink, xmax - shrink, ymin + shrink, ymax - shrink, 0.0, 0.0),
        eps=0.0)


def shortest_path(start, goal, margin: float = VEHICLE_RADIUS):
    """膨胀矩形角点可见性图 + Dijkstra.

    返回 (路径长度m, 路径点列); 不可达返回 None.
    注: 采用锐角膨胀矩形, 路径长度对圆盘真实最优最多保守 0.643m/拐角
    (圆角弧段可再省), 行程时间估计略偏保守, 不会违反约束.
    """
    start = (float(start[0]), float(start[1]))
    goal = (float(goal[0]), float(goal[1]))
    if not clearance_ok(start, margin):
        raise ValueError(f"起点 {start} 违反车辆几何约束")
    if not clearance_ok(goal, margin):
        raise ValueError(f"终点 {goal} 违反车辆几何约束")
    if math.dist(start, goal) < EPS:
        return 0.0, [start]
    rects = [inflate_rect(b.x0, b.x1, b.y0, b.y1, margin) for b in BUILDINGS]
    nodes = [start, goal]
    for r in rects:
        nodes += [(r[0], r[2]), (r[1], r[2]), (r[1], r[3]), (r[0], r[3])]
    n = len(nodes)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if all(not _seg_crosses_rect_interior(nodes[i], nodes[j], r)
                   for r in rects):
                d = math.dist(nodes[i], nodes[j])
                adj[i].append((j, d))
                adj[j].append((i, d))
    INF = float("inf")
    dist = [INF] * n
    prev = [None] * n
    dist[0] = 0.0
    heap = [(0.0, 0)]
    while heap:
        d, u = heapq.heappop(heap)
        if d > dist[u] + 1e-12:
            continue
        if u == 1:
            break
        for v, w in adj[u]:
            nd = d + w
            if nd < dist[v] - 1e-12:
                dist[v] = nd
                prev[v] = u
                heapq.heappush(heap, (nd, v))
    if math.isinf(dist[1]):
        return None
    path = []
    u = 1
    while u is not None:
        path.append(nodes[u])
        u = prev[u]
    path.reverse()
    return dist[1], path


def travel_time(p, q, margin: float = VEHICLE_RADIUS):
    """行程时间 s = 绕障最短路长 / 3; 不可达返回 None."""
    res = shortest_path(p, q, margin)
    return None if res is None else res[0] / SPEED
