# -*- coding: utf-8 -*-
"""
swb_common.py — Shared geometry / kinematics / scheduling core for
2026 "ShuWei Cup" Problem B: Rapid UAV Localization via Coordinated Mobile Vehicles.

Field: 100 m x 100 m, origin SW corner, x->east, y->north.
Vehicle: disk r=1.5 m, speed 3 m/s, camera at disk center, height 2 m.
Gimbal: elevation [-10, 80] deg, pan/tilt 30 deg/s simultaneous (duration = slower axis).
Travel attitude: due east / 0 tilt; slew only after full stop; after a 5 s joint
observation the gimbal resets to due-east/0, then the vehicle departs.
"""
import math
import numpy as np

# ----------------------------------------------------------------------
# 1. Nominal model parameters
# ----------------------------------------------------------------------
FIELD = 100.0
R_V = 1.5                      # vehicle disk radius (m)
CAM_H = 2.0                    # camera height (m)
V_SPEED = 3.0                  # travel speed (m/s)
EL_MIN, EL_MAX = -10.0, 80.0   # gimbal elevation range (deg)
SLEW_RATE = 30.0               # deg/s
OBS_DUR = 5.0                  # joint observation duration (s)
CENTER = (50.0, 50.0)
SIG_ANGLE = 1.0                # azimuth measurement sigma (deg), Problem 2

# Buildings B1..B5: (x0, x1, y0, y1, height)
BUILDINGS = [
    (22.0, 34.0, 16.0, 26.0, 20.0),
    (63.0, 79.0, 59.0, 71.0, 28.0),
    (10.0, 24.0, 65.0, 76.0, 12.0),
    (72.0, 86.0, 22.0, 34.0, 24.0),
    (41.0, 51.0, 38.0, 46.0, 32.0),
]
BUILD_NAMES = ["B1", "B2", "B3", "B4", "B5"]

V_INIT = [(15.0, 14.0), (88.0, 10.0), (14.0, 86.0), (88.0, 68.0), (50.0, 93.0)]
V_NAMES = ["V1", "V2", "V3", "V4", "V5"]

ALERTS = {"T1": (30.0, 40.0, 45.0), "T2": (150.0, 45.0, 35.0), "T3": (285.0, 38.0, 55.0)}

D2R = math.pi / 180.0
R2D = 180.0 / math.pi


def target_xyz(alert):
    th, r, h = alert
    return (CENTER[0] + r * math.cos(th * D2R), CENTER[1] + r * math.sin(th * D2R), h)


TARGETS = {k: target_xyz(v) for k, v in ALERTS.items()}


# ----------------------------------------------------------------------
# 2. Angle / distance helpers
# ----------------------------------------------------------------------
def wrap180(a):
    return (a + 180.0) % 360.0 - 180.0


def wrap360(a):
    return a % 360.0


def point_rect_dist(px, py, rect):
    x0, x1, y0, y1 = rect
    dx = max(x0 - px, 0.0, px - x1)
    dy = max(y0 - py, 0.0, py - y1)
    return math.hypot(dx, dy)


def point_rect_dist_vec(px, py, rect):
    x0, x1, y0, y1 = rect
    dx = np.maximum.reduce([x0 - px, np.zeros_like(px), px - x1])
    dy = np.maximum.reduce([y0 - py, np.zeros_like(px), py - y1])
    return np.hypot(dx, dy)


def seg_min_rect_dist(P, Q, rect):
    """Exact minimum point_rect_dist over the segment P->Q.
    g(t) = point_rect_dist(P + t(Q-P)) is piecewise quadratic with breakpoints
    where a coordinate crosses a slab boundary; the minimum over each piece is
    either at a breakpoint or at the quadratic vertex -> evaluated exactly."""
    x0, x1, y0, y1 = rect
    best = min(point_rect_dist(P[0], P[1], rect), point_rect_dist(Q[0], Q[1], rect))
    ts = {0.0, 1.0}
    for (p, q, lo, hi) in ((P[0], Q[0], x0, x1), (P[1], Q[1], y0, y1)):
        if abs(q - p) > 1e-15:
            for bound in (lo, hi):
                t = (bound - p) / (q - p)
                if 0.0 < t < 1.0:
                    ts.add(t)
    ts = sorted(ts)
    bx, by = Q[0] - P[0], Q[1] - P[1]
    for ta, tb in zip(ts, ts[1:]):
        tm = 0.5 * (ta + tb)
        ax, ay = P[0] + tm * bx, P[1] + tm * by
        dxm = max(x0 - ax, 0.0, ax - x1)
        dym = max(y0 - ay, 0.0, ay - y1)
        # vertex of (dxm + bx*t')^2 + (dym + by*t')^2 on this piece
        den = bx * bx + by * by
        cand = [ta, tb]
        if den > 1e-18:
            tv = -(dxm * bx + dym * by) / den
            if ta < tv < tb:
                cand.append(tv)
        for t in cand:
            px, py = P[0] + t * bx, P[1] + t * by
            best = min(best, point_rect_dist(px, py, rect))
    return best


# ----------------------------------------------------------------------
# 3. Segment vs box / rect tests
# ----------------------------------------------------------------------
def seg_hits_box(A, B, box):
    """3-D slab test segment A->B vs box=(x0,x1,y0,y1,z0,z1); tangency counts as hit."""
    tmin, tmax = 0.0, 1.0
    eps = 1e-12
    for i in range(3):
        a, b = A[i], B[i]
        lo, hi = box[2 * i], box[2 * i + 1]
        if abs(b - a) < eps:
            if a < lo - eps or a > hi + eps:
                return False
        else:
            t1, t2 = (lo - a) / (b - a), (hi - a) / (b - a)
            if t1 > t2:
                t1, t2 = t2, t1
            tmin = max(tmin, t1)
            tmax = min(tmax, t2)
            if tmin > tmax + eps:
                return False
    return tmax >= -eps and tmin <= 1.0 + eps


def seg_hits_rect2d(P, Q, rect, strict=False):
    """2-D slab test segment P->Q vs rect=(x0,x1,y0,y1).
    strict=False: touching counts as a hit (LOS rule).
    strict=True : only interior overlap blocks (path may hug an inflated wall)."""
    tmin, tmax = 0.0, 1.0
    eps = 1e-12
    for i in range(2):
        p, q = P[i], Q[i]
        lo, hi = rect[2 * i], rect[2 * i + 1]
        if abs(q - p) < eps:
            if p < lo - eps or p > hi + eps:
                return False
        else:
            t1, t2 = (lo - p) / (q - p), (hi - p) / (q - p)
            if t1 > t2:
                t1, t2 = t2, t1
            tmin = max(tmin, t1)
            tmax = min(tmax, t2)
            if tmin > tmax + eps:
                return False
    if strict:
        return (tmax - tmin) > 1e-9 and tmax > 1e-9 and tmin < 1.0 - 1e-9
    return tmax >= -eps and tmin <= 1.0 + eps


# ----------------------------------------------------------------------
# 4. Observability (vectorized)
# ----------------------------------------------------------------------
def elev_and_range(sx, sy, tgt):
    """LOS elevation (deg) and horizontal camera-target range for arrays or scalars."""
    tx, ty, H = tgt
    d = np.hypot(tx - sx, ty - sy)
    el = np.degrees(np.arctan2(H - CAM_H, d))
    return el, d


def feasibility_mask(sx, sy, tgt, margin=0.0):
    """Vectorized three-constraint station feasibility.
    sx, sy : numpy arrays of candidate station coordinates (disk center)
    tgt    : (tx, ty, H)
    Returns (ok, el, d)."""
    tx, ty, H = tgt
    sx = np.asarray(sx, float)
    sy = np.asarray(sy, float)
    ok = ((sx >= R_V + margin) & (sx <= FIELD - R_V - margin) &
          (sy >= R_V + margin) & (sy <= FIELD - R_V - margin))
    for b in BUILDINGS:
        ok &= point_rect_dist_vec(sx, sy, b[:4]) >= R_V - 1e-9
    el, d = elev_and_range(sx, sy, tgt)
    ok &= (el >= EL_MIN - 1e-9) & (el <= EL_MAX + 1e-9)
    # LOS occlusion: slab test per building, vectorized
    A = np.stack([sx, sy, np.full_like(sx, CAM_H)], axis=-1)
    B = np.array([tx, ty, H])
    D = B - A                                   # (N,3)
    for (x0, x1, y0, y1, h) in BUILDINGS:
        box = np.array([x0, x1, y0, y1, 0.0, h])
        tmin = np.zeros_like(sx)
        tmax = np.ones_like(sx)
        inside = np.ones_like(sx, dtype=bool)
        for i in range(3):
            a = A[:, i]
            b = B[i]
            lo, hi = box[2 * i], box[2 * i + 1]
            db = b - a
            parallel = np.abs(db) < 1e-12
            lo_i, hi_i = (a < lo - 1e-12) | (a > hi + 1e-12), None
            # parallel: must lie inside slab
            par_in = (a >= lo - 1e-12) & (a <= hi + 1e-12)
            inside &= np.where(parallel, par_in, True)
            with np.errstate(divide="ignore", invalid="ignore"):
                t1 = (lo - a) / db
                t2 = (hi - a) / db
            swap = t1 > t2
            t1s = np.where(swap, t2, t1)
            t2s = np.where(swap, t1, t2)
            tmin = np.maximum(tmin, np.where(parallel, -np.inf, t1s))
            tmax = np.minimum(tmax, np.where(parallel, np.inf, t2s))
        hit = inside & (tmax >= -1e-12) & (tmin <= 1.0 + 1e-12) & (tmin <= tmax + 1e-12)
        ok &= ~hit
    return ok, el, d


def min_range_for_target(H):
    """Minimum horizontal range imposed by the elevation cap: (H-CAM_H)/tan(80deg)."""
    return (H - CAM_H) / math.tan(EL_MAX * D2R)


def station_feasible(sx, sy, tgt):
    ok, el, d = feasibility_mask(np.array([sx]), np.array([sy]), tgt)
    return bool(ok[0]), float(el[0]), float(d[0])


def los_blocked(sx, sy, tx, ty, H):
    A = (sx, sy, CAM_H)
    B = (tx, ty, H)
    for (x0, x1, y0, y1, h) in BUILDINGS:
        if seg_hits_box(A, B, (x0, x1, y0, y1, 0.0, h)):
            return True
    return False


# ----------------------------------------------------------------------
# 5. Shortest collision-free paths (visibility graph, Floyd-Warshall precompute)
# ----------------------------------------------------------------------
class PathPlanner:
    """Shortest path LENGTH for the disk CENTER; buildings inflated by R_V
    (touching the inflated boundary is allowed)."""

    def __init__(self):
        self.rects = []
        for (x0, x1, y0, y1, h) in BUILDINGS:
            rx0, rx1 = max(x0 - R_V, R_V), min(x1 + R_V, FIELD - R_V)
            ry0, ry1 = max(y0 - R_V, R_V), min(y1 + R_V, FIELD - R_V)
            self.rects.append((rx0, rx1, ry0, ry1))
        nodes = []
        for (x0, x1, y0, y1) in self.rects:
            for cx in (x0, x1):
                for cy in (y0, y1):
                    nodes.append((cx, cy))
        self.nodes = nodes
        n = len(nodes)
        INF = float("inf")
        D = [[INF] * n for _ in range(n)]
        for i in range(n):
            D[i][i] = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if self._clear(self.nodes[i], self.nodes[j]):
                    d = math.dist(self.nodes[i], self.nodes[j])
                    D[i][j] = D[j][i] = d
        for k in range(n):
            Dk, dk = D, D[k]
            for i in range(n):
                Di = D[i]
                dik = Di[k]
                if dik == INF:
                    continue
                for j in range(n):
                    alt = dik + dk[j]
                    if alt < Di[j]:
                        Di[j] = alt
        self.D = D
        self.INF = INF
        self._cache = {}

    def _clear(self, P, Q):
        """Edge test: the segment must stay outside every building's R_V-
        NEIGHBOURHOOD (rounded/Minkowski rectangle = disk swept along the rect).
        Tangency to that neighbourhood boundary is allowed (the disk may touch
        the wall but not enter its interior).  Using the exact min-distance
        avoids the classic false-unreachable bug for legal points sitting in
        the corner quadrant of the sharp inflated rectangle."""
        for (px, py) in (P, Q):
            if not (R_V - 1e-9 <= px <= FIELD - R_V + 1e-9 and
                    R_V - 1e-9 <= py <= FIELD - R_V + 1e-9):
                return False
        for b in BUILDINGS:
            if seg_min_rect_dist(P, Q, b[:4]) < R_V - 1e-9:
                return False
        return True

    def dist(self, A, B):
        key = (round(A[0], 6), round(A[1], 6), round(B[0], 6), round(B[1], 6))
        v = self._cache.get(key)
        if v is not None:
            return v
        if A == B:
            self._cache[key] = 0.0
            return 0.0
        best = math.dist(A, B) if self._clear(A, B) else self.INF
        visA = [(i, math.dist(A, self.nodes[i])) for i in range(len(self.nodes))
                if self._clear(A, self.nodes[i])]
        visB = [(j, math.dist(B, self.nodes[j])) for j in range(len(self.nodes))
                if self._clear(B, self.nodes[j])]
        for i, da in visA:
            Di = self.D[i]
            for j, db in visB:
                alt = da + Di[j] + db
                if alt < best:
                    best = alt
        self._cache[key] = best
        return best

    def time(self, A, B):
        return self.dist(A, B) / V_SPEED


# ----------------------------------------------------------------------
# 6. Gimbal pointing / slewing
# ----------------------------------------------------------------------
def pointing_at(sx, sy, tgt):
    tx, ty, H = tgt
    pan = wrap360(math.degrees(math.atan2(ty - sy, tx - sx)))
    d = math.hypot(tx - sx, ty - sy)
    tilt = math.degrees(math.atan2(H - CAM_H, d))
    return pan, tilt


def slew_time(cur_pan, cur_tilt, pan, tilt):
    dp = abs(wrap180(pan - cur_pan))
    dt = abs(tilt - cur_tilt)
    return max(dp, dt) / SLEW_RATE


# ----------------------------------------------------------------------
# 7. Candidate stations
# ----------------------------------------------------------------------
def candidate_stations(tgt, ang_step=2.0, rad_step=1.0, r_max=85.0):
    """Polar sweep of feasible stations around the target's horizontal projection.
    Returns dict(pts=[(x,y)...], el=ndarray, r=ndarray)."""
    tx, ty, H = tgt
    r_min = min_range_for_target(H)
    angs = np.arange(0.0, 360.0, ang_step)
    rads = np.arange(r_min, r_max + 1e-9, rad_step)
    A, R = np.meshgrid(np.deg2rad(angs), rads, indexing="ij")
    X = tx + R * np.cos(A)
    Y = ty + R * np.sin(A)
    ok, el, _ = feasibility_mask(X.ravel(), Y.ravel(), tgt)
    pts = np.stack([X.ravel()[ok], Y.ravel()[ok]], axis=1)
    els = el.ravel()[ok]
    rrs = R.ravel()[ok]
    return dict(pts=[(float(x), float(y)) for x, y in pts],
                el=els, r=rrs)


# ----------------------------------------------------------------------
# 8. Timeline simulation
# ----------------------------------------------------------------------
def simulate(assign, stations, planner, init_pos=None, init_state=None):
    """Exact timeline for a fixed assignment + station choice (fixed point on
    joint observation starts, because a vehicle's departure after its k-th task
    depends on the partner's ready time of that alert).
    assign   : {veh: [alert, ...]} ordered tasks
    stations : {(veh, alert): (x, y)}
    init_pos : {veh: (x, y)} default V_INIT
    init_state: {veh: dict(t, pos, pan, tilt)} overrides init_pos (for P4)
    Returns dict(timeline, alert_start, max_start, makespan, ready, ok)."""
    vehs = [v for v in assign if assign[v]]
    base = {}
    for v in vehs:
        if init_state and v in init_state:
            s0 = init_state[v]
            base[v] = dict(t=float(s0["t"]), pos=tuple(s0["pos"]),
                           pan=float(s0["pan"]), tilt=float(s0["tilt"]))
        else:
            p0 = init_pos[v] if init_pos else V_INIT[int(v[1]) - 1]
            base[v] = dict(t=0.0, pos=tuple(p0), pan=0.0, tilt=0.0)

    ready = {}
    cur = {v: base[v] for v in vehs}
    for _ in range(15):                       # fixed point (monotone increasing)
        changed = False
        for v in vehs:
            c = dict(base[v])
            for a in assign[v]:
                s = stations[(v, a)]
                t_arr = c["t"] + planner.time(c["pos"], s)
                pan, tilt = pointing_at(s[0], s[1], TARGETS[a])
                rdy = t_arr + slew_time(c["pan"], c["tilt"], pan, tilt)
                if ready.get((v, a)) != rdy:
                    ready[(v, a)] = rdy
                    changed = True
                partners = [u for u in vehs if u != v and a in assign[u]]
                start = rdy
                for u in partners:
                    ru = ready.get((u, a))
                    if ru is not None:
                        start = max(start, ru)
                reset = slew_time(pan, tilt, 0.0, 0.0)
                c = dict(t=start + OBS_DUR + reset, pos=(s[0], s[1]),
                         pan=0.0, tilt=0.0)
            cur[v] = c
        if not changed:
            break

    alerts = sorted({a for v in vehs for a in assign[v]})
    alert_start = {}
    for a in alerts:
        vals = [ready[(v, a)] for v in vehs if a in assign[v]]
        alert_start[a] = max(vals)

    timeline = {}
    for v in vehs:
        c = dict(base[v])
        ev = [dict(t=c["t"], kind="available", pos=c["pos"])]
        for a in assign[v]:
            s = stations[(v, a)]
            t_arr = c["t"] + planner.time(c["pos"], s)
            pan, tilt = pointing_at(s[0], s[1], TARGETS[a])
            rdy = t_arr + slew_time(c["pan"], c["tilt"], pan, tilt)
            start = alert_start[a]
            end = start + OBS_DUR
            reset = slew_time(pan, tilt, 0.0, 0.0)
            ev += [
                dict(t=c["t"], kind="depart", pos=c["pos"], alert=a),
                dict(t=t_arr, kind="arrive", pos=s, alert=a),
                dict(t=rdy, kind="ready", pan=pan, tilt=tilt, alert=a),
                dict(t=start, kind="obs_start", alert=a),
                dict(t=end, kind="obs_end", alert=a),
                dict(t=end + reset, kind="reset_done", alert=a),
            ]
            c = dict(t=end + reset, pos=(s[0], s[1]), pan=0.0, tilt=0.0)
        timeline[v] = ev
    makespan = max((ev[-1]["t"] for ev in timeline.values()), default=0.0)
    return dict(timeline=timeline, alert_start=alert_start,
                max_start=max(alert_start.values()) if alert_start else 0.0,
                makespan=makespan, ready=ready)


def approx_objective(assign, stations, planner, init_pos=None, init_state=None):
    """Cheap single-pass optimistic estimate of max alert start (no partner waits):
    used only to rank assignments before exact simulation."""
    vehs = [v for v in assign if assign[v]]
    starts = []
    for v in vehs:
        if init_state and v in init_state:
            s0 = init_state[v]
            c = dict(t=float(s0["t"]), pos=tuple(s0["pos"]),
                     pan=float(s0["pan"]), tilt=float(s0["tilt"]))
        else:
            p0 = init_pos[v] if init_pos else V_INIT[int(v[1]) - 1]
            c = dict(t=0.0, pos=tuple(p0), pan=0.0, tilt=0.0)
        for a in assign[v]:
            s = stations[(v, a)]
            t_arr = c["t"] + planner.time(c["pos"], s)
            pan, tilt = pointing_at(s[0], s[1], TARGETS[a])
            rdy = t_arr + slew_time(c["pan"], c["tilt"], pan, tilt)
            starts.append(rdy)
            reset = slew_time(pan, tilt, 0.0, 0.0)
            c = dict(t=rdy + OBS_DUR + reset, pos=(s[0], s[1]), pan=0.0, tilt=0.0)
    return max(starts) if starts else 0.0


# ----------------------------------------------------------------------
# 9. Assignment enumeration and station search
# ----------------------------------------------------------------------
def enumerate_assignments(alerts, vehs, max_tasks=3):
    """All {veh: ordered [alerts...]} with every alert given to two DISTINCT vehicles.
    max_tasks caps per-vehicle load (pruning)."""
    from itertools import combinations, permutations
    results = {}
    pairs = list(combinations(range(len(vehs)), 2))
    cur = {v: [] for v in vehs}

    def expand(i):
        if i == len(alerts):
            keys = [v for v in vehs if cur[v]]
            if any(len(cur[v]) > max_tasks for v in keys):
                return
            opts = []
            for v in keys:
                seen = sorted(set(permutations(cur[v])))
                opts.append((v, seen))

            def rec2(j, acc):
                if j == len(opts):
                    key = tuple((v, tuple(acc[v])) for v in vehs)
                    if key not in results:
                        results[key] = {v: list(acc[v]) for v in vehs}
                    return
                v, lst = opts[j]
                for order in lst:
                    acc[v] = list(order)
                    rec2(j + 1, acc)
                    acc[v] = []
            rec2(0, {v: [] for v in vehs})
            return
        a = alerts[i]
        for ia, ib in pairs:
            va, vb = vehs[ia], vehs[ib]
            if len(cur[va]) + 1 > max_tasks or len(cur[vb]) + 1 > max_tasks:
                continue
            cur[va].append(a)
            cur[vb].append(a)
            expand(i + 1)
            cur[va].pop()
            cur[vb].pop()

    expand(0)
    return list(results.values())


def _station_candidates_for(state, alert, cand, planner, topk=12):
    """Top-k candidate stations for a vehicle state by (exact path time + slew)."""
    tgt = TARGETS[alert]
    pts = cand[alert]["pts"]
    # pre-rank by Euclidean lower bound, then exact-path the best few
    eu = sorted(range(len(pts)), key=lambda i: math.dist(state["pos"], pts[i]))
    scored = []
    for i in eu[:40]:
        s = pts[i]
        tt = planner.time(state["pos"], s)
        pan, tilt = pointing_at(s[0], s[1], tgt)
        sl = slew_time(state["pan"], state["tilt"], pan, tilt)
        scored.append((tt + sl, s, pan, tilt))
    scored.sort(key=lambda x: x[0])
    return scored[:topk]


def greedy_stations(assign, planner, cand, init_pos=None, init_state=None):
    """Myopic initial station choice (min own travel+slew per task, round-robin)."""
    stations = {}
    vehs = [v for v in assign if assign[v]]
    state = {}
    for v in vehs:
        if init_state and v in init_state:
            s0 = init_state[v]
            state[v] = dict(t=float(s0["t"]), pos=tuple(s0["pos"]),
                            pan=float(s0["pan"]), tilt=float(s0["tilt"]))
        else:
            p0 = init_pos[v] if init_pos else V_INIT[int(v[1]) - 1]
            state[v] = dict(t=0.0, pos=tuple(p0), pan=0.0, tilt=0.0)
    pending = {v: list(assign[v]) for v in vehs}
    progressed = True
    while progressed:
        progressed = False
        for v in vehs:
            while pending[v]:
                a = pending[v][0]
                key = (v, a)
                if key not in stations:
                    best = _station_candidates_for(state[v], a, cand, planner, topk=1)[0]
                    stations[key] = best[1]
                s = stations[key]
                t_arr = state[v]["t"] + planner.time(state[v]["pos"], s)
                pan, tilt = pointing_at(s[0], s[1], TARGETS[a])
                rdy = t_arr + slew_time(state[v]["pan"], state[v]["tilt"], pan, tilt)
                partner = [u for u in vehs if u != v and a in assign[u]]
                start = rdy
                if partner:
                    pk = (partner[0], a)
                    if pk in stations and not pending[partner[0]]:
                        # partner fully committed -> its ready time known from simulate later
                        start = rdy
                state[v] = dict(t=start + OBS_DUR + slew_time(pan, tilt, 0.0, 0.0),
                                pos=(s[0], s[1]), pan=0.0, tilt=0.0)
                pending[v].pop(0)
                progressed = True
    return stations


def local_search(assign, stations, planner, cand, objective, init_pos=None,
                 init_state=None, passes=3, topk=10):
    """Coordinate descent over station slots using the exact simulate() objective.
    Alternatives for a slot are the candidates nearest to the slot's current station."""
    best_st = dict(stations)
    best_val = objective(simulate(assign, best_st, planner, init_pos, init_state))
    slots = [(v, a) for v in assign if assign[v] for a in assign[v]]
    for _ in range(passes):
        improved = False
        for (v, a) in slots:
            cur_s = best_st[(v, a)]
            pts = cand[a]["pts"]
            near = sorted(pts, key=lambda p: math.dist(p, cur_s))[:topk]
            for s in near:
                if s == cur_s:
                    continue
                trial = dict(best_st)
                trial[(v, a)] = s
                val = objective(simulate(assign, trial, planner, init_pos, init_state))
                if val < best_val - 1e-9:
                    best_val, best_st, improved = val, trial, True
                    cur_s = s
        if not improved:
            break
    return best_st, best_val


def sched_optimize(alerts, vehs, planner, cand, init_pos=None, init_state=None,
                   objective=None, max_tasks=3, prefilter=120, ls_top=6):
    """Enumerate assignments -> optimistic pre-rank -> exact simulate on the best
    candidates -> local station search. Returns (best_assign, best_stations, best_sim)."""
    if objective is None:
        objective = lambda sim: sim["max_start"]
    assigns = enumerate_assignments(alerts, vehs, max_tasks=max_tasks)
    approx = []
    veh_init = {}
    for v in vehs:
        if init_state and v in init_state:
            veh_init[v] = tuple(init_state[v]["pos"])
        else:
            veh_init[v] = tuple(init_pos[v]) if init_pos else V_INIT[int(v[1]) - 1]
    # pre-rank by a travel lower bound: sum over tasks of the nearest-candidate
    # Euclidean distance / speed (cheaper than path planning for every assignment)
    for asg in assigns:
        lb = 0.0
        for v in vehs:
            c = veh_init[v]
            for a in asg[v]:
                pts = cand[a]["pts"]
                lb += (min(math.dist(c, p) for p in pts) / V_SPEED if pts else 0.0)
                c = (TARGETS[a][0], TARGETS[a][1])
        approx.append((lb, asg))
    approx.sort(key=lambda x: x[0])
    pool = []
    for _, asg in approx[:prefilter]:
        st = greedy_stations(asg, planner, cand, init_pos, init_state)
        sim = simulate(asg, st, planner, init_pos, init_state)
        pool.append((objective(sim), asg, st))
    pool.sort(key=lambda x: x[0])
    best = None
    for val0, asg, st in pool[:8]:
        st2, val2 = local_search(asg, st, planner, cand, objective, init_pos,
                                 init_state, passes=3, topk=ls_top)
        if best is None or val2 < best[0]:
            sim2 = simulate(asg, st2, planner, init_pos, init_state)
            best = (val2, asg, st2, sim2)
    return best[1], best[2], best[3]


def nearest_heuristic(alerts, vehs, planner, cand, alert_order=None):
    """Baseline 'nearest-assignment heuristic': alerts are processed in a fixed
    order; each alert takes the two vehicles with the smallest optimistic
    earliest-ready time (travel from the vehicle's committed position to its
    nearest feasible station + slew); the vehicle's task is appended to its chain."""
    order = alert_order or list(alerts)
    chains = {v: [] for v in vehs}
    stations = {}
    for a in order:
        score = []
        for v in vehs:
            # vehicle current committed end state (optimistic: own timeline only)
            c = dict(t=0.0, pos=V_INIT[int(v[1]) - 1], pan=0.0, tilt=0.0)
            for a2 in chains[v]:
                s = stations[(v, a2)]
                t_arr = c["t"] + planner.time(c["pos"], s)
                pan, tilt = pointing_at(s[0], s[1], TARGETS[a2])
                rdy = t_arr + slew_time(c["pan"], c["tilt"], pan, tilt)
                reset = slew_time(pan, tilt, 0.0, 0.0)
                c = dict(t=rdy + OBS_DUR + reset, pos=(s[0], s[1]), pan=0.0, tilt=0.0)
            best = _station_candidates_for(c, a, cand, planner, topk=1)[0]
            score.append((best[0], v, best[1]))
        score.sort()
        chosen = [score[0][1]]
        # second vehicle: next best that is not the same vehicle
        for item in score[1:]:
            if item[1] not in chosen:
                chosen.append(item[1])
                break
        for v in chosen:
            chains[v].append(a)
            stations[(v, a)] = [x for x in score if x[1] == v][0][2]
    return chains, stations
