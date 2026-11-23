# -*- coding: utf-8 -*-
"""
Problem 2(2) relaxed - time-first, sigma < 2.0 m, ordinary overlap allowed.

This version is standalone and does NOT use the hand-written frozen plans.
It starts from the time-optimal station assignment, checks a strict
sigma < 2.0 m gate, and repairs only the violating alerts (latest first).
It still enumerates all assignment/order combinations, so a vehicle may be
assigned to its second-nearest alert if that improves the global schedule.
The ordinary alert is allowed to overlap the reused vehicle's second
observation; the chain is reported as a diagnostic, not a hard constraint.
A reused vehicle may also stay at its previous station when that station is
observable for its next alert, which removes unnecessary transfers.
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import swb_common as S
import p22_visuals as V

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

ALERTS = ["T1", "T2", "T3"]
VEHS = ["V1", "V2", "V3", "V4", "V5"]
planner = S.PathPlanner()

CLR_MIN = 0.10
EL_MAX = 79.9
SIGMA_MAX = 2.0       # strict upper bound; accept only sigma < 2.0 m
GAMMA_MIN = 30.0
BEARING_SIG = math.radians(1.0)
N_REPAIR = 300        # time-ranked assignments repaired
PAIR_TOPK = 300
LINE_STEP = 1.0
LINE_RANGE = 90.0


# ----------------------------------------------------------------------
# 1. Stations and quality
# ----------------------------------------------------------------------
def clearance_at(p):
    x, y = p
    d = min(x - S.R_V, S.FIELD - S.R_V - x, y - S.R_V, S.FIELD - S.R_V - y)
    for b in S.BUILDINGS:
        d = min(d, float(S.point_rect_dist_vec(
            np.array([x]), np.array([y]), b[:4])[0]) - S.R_V)
    return d


def station_ok(p, a):
    if clearance_at(p) < CLR_MIN - 1e-9:
        return False
    ok, el, _ = S.station_feasible(p[0], p[1], S.TARGETS[a])
    return bool(ok) and el <= EL_MAX + 1e-9


def make_set(a, ang, rad, rmax):
    full = S.candidate_stations(S.TARGETS[a], ang_step=ang,
                                rad_step=rad, r_max=rmax)
    pts, rs = [], []
    for p, r in zip(full["pts"], full["r"]):
        if station_ok(p, a):
            pts.append(p)
            rs.append(float(r))
    return dict(pts=pts, arr=np.asarray(pts, float) if pts else np.zeros((0, 2)),
                r=np.asarray(rs, float))


print("stations ...")
C_TIME = {a: make_set(a, 8.0, 2.0, 80.0) for a in ALERTS}
C_FULL = {a: make_set(a, 2.0, 1.0, 80.0) for a in ALERTS}
C_QUAL = {a: make_set(a, 3.0, 1.5, 42.0) for a in ALERTS}

# add per-vehicle nearest candidates to the time set
for a in ALERTS:
    pts = list(C_TIME[a]["pts"])
    for i, v in enumerate(VEHS):
        d = np.hypot(C_FULL[a]["arr"][:, 0] - S.V_INIT[i][0],
                     C_FULL[a]["arr"][:, 1] - S.V_INIT[i][1])
        for j in np.argsort(d)[:50]:
            pts.append(tuple(C_FULL[a]["arr"][j]))
    seen, uniq = set(), []
    for p in pts:
        k = (round(p[0], 6), round(p[1], 6))
        if k not in seen:
            seen.add(k)
            uniq.append(p)
    C_TIME[a] = dict(pts=uniq, arr=np.asarray(uniq, float),
                     r=np.array([math.hypot(S.TARGETS[a][0] - p[0],
                                            S.TARGETS[a][1] - p[1])
                                 for p in uniq]))
for a in ALERTS:
    print(f"  {a}: time={len(C_TIME[a]['pts'])}, quality={len(C_QUAL[a]['pts'])}")


def pair_quality(s1, s2, a):
    tx, ty, H = S.TARGETS[a]
    r1 = math.hypot(tx - s1[0], ty - s1[1])
    r2 = math.hypot(tx - s2[0], ty - s2[1])
    u1 = ((tx - s1[0]) / r1, (ty - s1[1]) / r1)
    u2 = ((tx - s2[0]) / r2, (ty - s2[1]) / r2)
    c = max(-1.0, min(1.0, u1[0] * u2[0] + u1[1] * u2[1]))
    g = math.degrees(math.acos(c))
    sig = BEARING_SIG * math.sqrt(r1 * r1 + r2 * r2) / max(
        math.sin(math.radians(g)), 1e-9)
    return dict(r1=r1, r2=r2, gamma=g, sigma=sig)


def quality_info(st, assign):
    info, bad = {}, []
    for a in ALERTS:
        vs = [v for v in VEHS if a in assign[v]]
        if len(vs) != 2:
            continue
        q = pair_quality(st[(vs[0], a)], st[(vs[1], a)], a)
        info[a] = q
        if q["sigma"] >= SIGMA_MAX - 1e-9 or q["gamma"] < GAMMA_MIN - 1e-9:
            bad.append(a)
    return info, bad


# ----------------------------------------------------------------------
# 2. Time-optimal initial stations (old 2.2 style)
# ----------------------------------------------------------------------
def nearest_time(state, a, topn=25):
    cand = C_TIME[a]
    arr = cand["arr"]
    d = np.hypot(arr[:, 0] - state["pos"][0], arr[:, 1] - state["pos"][1])
    best = None
    # A reused vehicle may stay at its previous station if that station is
    # also observable for the next alert.  This is the cheapest transfer.
    if station_ok(state["pos"], a):
        pan, tilt = S.pointing_at(state["pos"][0], state["pos"][1], S.TARGETS[a])
        sl = S.slew_time(state["pan"], state["tilt"], pan, tilt)
        best = (sl, (state["pos"][0], state["pos"][1]), pan, tilt)
    for i in np.argsort(d)[:min(topn, len(d))]:
        p = (float(arr[i, 0]), float(arr[i, 1]))
        tt = planner.time(state["pos"], p)
        pan, tilt = S.pointing_at(p[0], p[1], S.TARGETS[a])
        sl = S.slew_time(state["pan"], state["tilt"], pan, tilt)
        val = tt + sl
        if best is None or val < best[0]:
            best = (val, p, pan, tilt)
    return best


def time_stations(assign):
    st, state = {}, {}
    for v in VEHS:
        state[v] = dict(t=0.0, pos=S.V_INIT[VEHS.index(v)], pan=0.0, tilt=0.0)
    for v in VEHS:
        for a in assign[v]:
            _, p, pan, tilt = nearest_time(state[v], a)
            st[(v, a)] = p
            t_arr = state[v]["t"] + planner.time(state[v]["pos"], p)
            rdy = t_arr + S.slew_time(state[v]["pan"], state[v]["tilt"], pan, tilt)
            state[v] = dict(t=rdy + S.OBS_DUR + S.slew_time(pan, tilt, 0.0, 0.0),
                            pos=p, pan=0.0, tilt=0.0)
    return st


# ----------------------------------------------------------------------
# 3. Quality pair pool (used only for repairing violating alerts)
# ----------------------------------------------------------------------
_CACHE = {}


def perp_points(a, sL):
    key = (a, round(sL[0], 3), round(sL[1], 3))
    if key in _CACHE:
        return _CACHE[key]
    tx, ty, H = S.TARGETS[a]
    ux, uy = sL[0] - tx, sL[1] - ty
    n = math.hypot(ux, uy)
    if n < 1e-9:
        return []
    dx, dy = -uy / n, ux / n
    out = []
    for t in np.arange(-LINE_RANGE, LINE_RANGE + 1e-9, LINE_STEP):
        p = (tx + t * dx, ty + t * dy)
        if station_ok(p, a):
            out.append(p)
    _CACHE[key] = out
    return out


def build_pool(a):
    pts = [p for p, r in zip(C_QUAL[a]["pts"], C_QUAL[a]["r"]) if r <= 38.0]
    if len(pts) > 55:
        pts = [pts[i] for i in np.linspace(0, len(pts) - 1, 55, dtype=int)]
    pairs = {}
    for pL in pts:
        for pF in perp_points(a, pL):
            q = pair_quality(pL, pF, a)
            if q["sigma"] >= SIGMA_MAX - 1e-9 or q["gamma"] < GAMMA_MIN - 1e-9:
                continue
            key = tuple(sorted([(round(pL[0], 3), round(pL[1], 3)),
                                (round(pF[0], 3), round(pF[1], 3))]))
            rec = (max(q["r1"], q["r2"]), q["sigma"], pL, pF, q)
            if key not in pairs or rec[:2] < pairs[key][:2]:
                pairs[key] = rec
    return sorted(pairs.values(), key=lambda x: (x[0], x[1]))[:PAIR_TOPK]


PAIR_POOL = {a: build_pool(a) for a in ALERTS}
for a in ALERTS:
    print(f"  {a}: repair pairs={len(PAIR_POOL[a])}")


def depart_event(sim, v, a):
    for e in sim["timeline"][v]:
        if e["kind"] == "depart" and e["alert"] == a:
            return e
    return None


def candidate_list(state, a, topn=50):
    """Combined time-near and target-near candidates for one vehicle."""
    arr = C_FULL[a]["arr"]
    d = np.hypot(arr[:, 0] - state["pos"][0], arr[:, 1] - state["pos"][1])
    idx = np.argsort(d)[:min(topn, len(d))]
    return [(float(arr[i, 0]), float(arr[i, 1])) for i in idx]


def best_repair(assign, st, a, sim):
    """Search all candidate pairs for a violating alert.  No 90 deg requirement:
    any pair with sigma < 2 m and gamma >= 30 deg is allowed."""
    vs = [v for v in VEHS if a in assign[v]]
    if len(vs) != 2:
        return None
    e1, e2 = depart_event(sim, vs[0], a), depart_event(sim, vs[1], a)
    if e1 is None or e2 is None:
        return None
    other = {k: t for k, t in sim["alert_start"].items() if k != a}
    cand1 = candidate_list(e1, a, topn=50)
    cand2 = candidate_list(e2, a, topn=50)
    best = None
    for p1 in cand1:
        pan1, tilt1 = S.pointing_at(p1[0], p1[1], S.TARGETS[a])
        r1 = e1["t"] + planner.time(e1["pos"], p1) + S.slew_time(
            0.0, 0.0, pan1, tilt1)
        for p2 in cand2:
            q = pair_quality(p1, p2, a)
            if q["sigma"] >= SIGMA_MAX - 1e-9 or q["gamma"] < GAMMA_MIN - 1e-9:
                continue
            pan2, tilt2 = S.pointing_at(p2[0], p2[1], S.TARGETS[a])
            r2 = e2["t"] + planner.time(e2["pos"], p2) + S.slew_time(
                0.0, 0.0, pan2, tilt2)
            start = max(r1, r2)
            newmax = max(list(other.values()) + [start]) if other else start
            val = (newmax, start, q["sigma"])
            if best is None or val < best[0]:
                best = (val, p1, p2)
    return None if best is None else best[1:]


def chain_info(assign, st, sim):
    """Explicitly check the user's chain structure:
    one reused vehicle r serves A1 then A2; P1 helps A1, P2 pre-positions for
    A2, and the ordinary alert A3 must finish before A2 starts."""
    doubles = [v for v in VEHS if len(assign[v]) == 2]
    if len(doubles) != 1:
        return None
    rv = doubles[0]
    a1, a2 = assign[rv]
    rest = [a for a in ALERTS if a not in (a1, a2)]
    if len(rest) != 1:
        return None
    a3 = rest[0]
    p1 = [v for v in VEHS if v != rv and a1 in assign[v]]
    p2 = [v for v in VEHS if v != rv and a2 in assign[v]]
    if len(p1) != 1 or len(p2) != 1:
        return None
    p1, p2 = p1[0], p2[0]
    s1, s2 = st[(rv, a1)], st[(rv, a2)]
    rv_path = (planner.time(S.V_INIT[VEHS.index(rv)], s1)
               + planner.time(s1, s2))
    s2_start = sim["alert_start"][a2]
    ordinary_end = sim["alert_start"][a3] + S.OBS_DUR
    p2_ready = sim["ready"][(p2, a2)]
    p1_ready = sim["ready"][(p1, a1)]
    return dict(
        rv=rv, a1=a1, a2=a2, a3=a3, p1=p1, p2=p2,
        rv_path_m=rv_path * S.V_SPEED, rv_path_s=rv_path,
        s1=sim["alert_start"][a1], s2=s2_start,
        s3=sim["alert_start"][a3],
        ordinary_end=ordinary_end, ordinary_slack=s2_start - ordinary_end,
        p1_ready=p1_ready, p1_slack=sim["alert_start"][a1] - p1_ready,
        p2_ready=p2_ready, p2_slack=s2_start - p2_ready,
        chain_ok=(ordinary_end <= s2_start + 1e-9
                  and p2_ready <= s2_start + 1e-9),
    )


def repair_assignment(assign, st):
    cur = dict(st)
    sim = S.simulate(assign, cur, planner)
    info, bad = quality_info(cur, assign)
    if not bad:
        return cur, sim, info, 0.0
    for _ in range(3):
        info, bad = quality_info(cur, assign)
        if not bad:
            break
        # Repair the latest violating alert first (the likely bottleneck).
        a = max(bad, key=lambda x: sim["alert_start"].get(x, 0.0))
        res = best_repair(assign, cur, a, sim)
        if res is None:
            break
        vs = [v for v in VEHS if a in assign[v]]
        cur[(vs[0], a)], cur[(vs[1], a)] = res
        sim = S.simulate(assign, cur, planner)
    info, bad = quality_info(cur, assign)
    pen = 1000.0 * len(bad) + sum(
        max(0.0, q["sigma"] - (SIGMA_MAX - 1e-9)) for q in info.values())
    return cur, sim, info, pen


# ----------------------------------------------------------------------
# 4. Enumerate and repair
# ----------------------------------------------------------------------
print("enumerating assignments ...")
assigns = S.enumerate_assignments(ALERTS, VEHS, max_tasks=3)
print("  combinations:", len(assigns))

rows = []
for asg in assigns:
    st = time_stations(asg)
    sim = S.simulate(asg, st, planner)
    rows.append((sim["max_start"], asg, st, sim))
rows.sort(key=lambda x: x[0])
print("  time-ranked top:", rows[0][0], "s")

best = None
best_chain = None
for _, asg, st0, _ in rows[:N_REPAIR]:
    st, sim, info, pen = repair_assignment(asg, st0)
    feasible = pen <= 1e-9
    ci = chain_info(asg, st, sim)
    chain_ok = bool(ci and ci["chain_ok"])
    cand = (not feasible, sim["max_start"],
            sum(q["sigma"] for q in info.values()), asg, st, sim, info, pen)
    if best is None or cand[:3] < best[:3]:
        best = cand
        best_chain = ci

_, max_start, _, best_asg, best_st, best_sim, best_info, best_pen = best
feasible = best_pen <= 1e-9
print("\n=== HYBRID (time-first, repair violators) ===")
print("assignment:", {v: a for v, a in best_asg.items() if a})
print("stations:", {f"{v}-{a}": tuple(round(float(x), 2) for x in p)
                   for (v, a), p in best_st.items()})
print("starts:", {k: round(t, 3) for k, t in best_sim["alert_start"].items()})
for a in ALERTS:
    q = best_info[a]
    print(f"  {a}: gamma={q['gamma']:.1f}, sigma={q['sigma']:.3f} m")
max_sigma = max(q["sigma"] for q in best_info.values())
strict_sigma_ok = all(q["sigma"] < SIGMA_MAX - 1e-9
                      for q in best_info.values())
print(f"max start={max_start:.3f} s, makespan={best_sim['makespan']:.3f} s, "
      f"feasible={feasible}")
print(f"strict sigma < {SIGMA_MAX:.2f} m: {strict_sigma_ok} "
      f"(max {max_sigma:.3f} m)")
chain_ok = bool(best_chain and best_chain["chain_ok"])
accepted = feasible and strict_sigma_ok
print(f"chain diagnostic: chain_ok={chain_ok} (not required), accepted={accepted}")
if best_chain:
    c = best_chain
    print(f"chain: RV={c['rv']}, {c['a1']} -> {c['a2']}, "
          f"P1={c['p1']}, P2={c['p2']}, ordinary={c['a3']}")
    print(f"  RV path={c['rv_path_m']:.3f} m ({c['rv_path_s']:.3f} s), "
          f"ordinary ends {c['ordinary_end']:.3f} <= {c['s2']:.3f} "
          f"(slack {c['ordinary_slack']:.3f} s)")
    print(f"  P1 ready={c['p1_ready']:.3f} (waits {c['p1_slack']:.3f} s), "
          f"P2 ready={c['p2_ready']:.3f} (waits {c['p2_slack']:.3f} s)")

# ----------------------------------------------------------------------
# 5. Figures
# ----------------------------------------------------------------------
V.draw_gantt(
    best_sim,
    f"2 m schedule: max start = {max_start:.2f} s, sigma < 2.0 m",
    OUT / "fig_p2_2_2m_gantt.png",
)
V.draw_trajectory(
    best_st, best_asg, best_info, max_start,
    OUT / "fig_p2_2_2m_trajectory.png",
    f"2 m schedule: routes, stations and observation rays (max start = {max_start:.2f} s)",
)

# ----------------------------------------------------------------------
# 6. JSON
# ----------------------------------------------------------------------
def timeline_json(sim):
    out = {}
    for v, evs in sim["timeline"].items():
        rows2 = []
        for e in evs:
            row = dict(t=round(e["t"], 3), kind=e["kind"])
            if "alert" in e:
                row["alert"] = e["alert"]
            if "pos" in e:
                row["pos"] = [round(e["pos"][0], 3), round(e["pos"][1], 3)]
            if "pan" in e:
                row["pan_deg"] = round(e["pan"], 3)
                row["tilt_deg"] = round(e["tilt"], 3)
            rows2.append(row)
        out[v] = rows2
    return out


res = {
    "problem": "2(2)",
    "mode": "time-first; repair only quality-violating alerts",
    "constraints": {
        "clearance_min_m": CLR_MIN,
        "elevation_max_deg": EL_MAX,
        "sigma_max_m": SIGMA_MAX,
        "gamma_min_deg": GAMMA_MIN,
        "travel_cap_s": None,
    },
    "assignment": {v: a for v, a in best_asg.items() if a},
    "stations": {f"{v}-{a}": [round(float(p[0]), 4), round(float(p[1]), 4)]
                 for (v, a), p in best_st.items()},
    "alert_quality": {a: {k: round(float(v), 4) for k, v in q.items()}
                      for a, q in best_info.items()},
    "alert_start_s": {k: round(t, 3) for k, t in best_sim["alert_start"].items()},
    "max_start_s": round(max_start, 3),
    "makespan_s": round(best_sim["makespan"], 3),
    "max_sigma_m": round(max_sigma, 4),
    "strict_sigma_ok": bool(strict_sigma_ok),
    "chain_ok": bool(chain_ok),
    "accepted": bool(accepted),
    "quality_penalty": round(best_pen, 4),
    "chain": best_chain,
    "timeline": timeline_json(best_sim),
}
with open(OUT / "problem2_2_results_2m.json", "w", encoding="utf-8") as fp:
    json.dump(res, fp, ensure_ascii=False, indent=2)
print("results saved:", OUT / "problem2_2_results_2m.json")
