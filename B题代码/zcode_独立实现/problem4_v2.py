# -*- coding: utf-8 -*-
"""Problem 4 (v2) — complete per-question solver (unified v2 conventions).

Reuses the verified machinery of problem4.py (state capture, commitment-aware
rescheduling, failure-modified stubs) as a library, and adds:
  1. Assumption-7 margins enforced on ALL candidate stations
     (clearance >= 0.10 m beyond the buffer, elevation <= 79.9 deg);
  2. the P2 rule-based greedy strategy (paper's third strategy);
  3. the full 86-scenario family of the paper's §Q4 Table;
  4. +/-10% kinematic sensitivity of the P1 response (3 representative scenarios);
  5. outputs: q4_v2_results.json + four 300-dpi figures.

Baseline: the margin-compliant min-max schedule of Q2(2) (re-derived here).
All scenarios run under three strategies: P0 retain / P1 event-triggered
re-optimization / P2 rule-based greedy; metrics (eta, t_last), window [0,60] s.
"""
import json
import math
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
import problem4 as P4   # library import (machinery; its own suite is __main__-guarded)
import swb_common as S  # noqa: E402

OUT = BASE / "results"
OUT.mkdir(exist_ok=True)

WINDOW, T_D = 60.0, 5.0
MARGIN, EL_MAX = 0.10, 79.9   # Assumption-7 operating margins
VEHS = P4.VEHS
PLANNER = P4.planner


# ----------------------------------------------------------------------
# Assumption-7 margin filter on candidate stations
# ----------------------------------------------------------------------
def clearance_beyond_buffer(x, y):
    c = min(x - S.R_V, S.FIELD - S.R_V - x, y - S.R_V, S.FIELD - S.R_V - y)
    for (x0, x1, y0, y1, _h) in S.BUILDINGS:
        dx = max(x0 - x, 0.0, x - x1)
        dy = max(y0 - y, 0.0, y - y1)
        c = min(c, math.hypot(dx, dy) - S.R_V)
    return c


def margin_full_candidates(target, vehicle_positions):
    """Full sweep, margin+elevation filtered, per-vehicle-near + close-in subsets."""
    tx, ty, H = target
    r_min = (H - S.CAM_H) / math.tan(math.radians(S.EL_MAX))
    pts_all, r_all = [], []
    angs = np.arange(0.0, 360.0, 2.0)
    ca, sa = np.cos(np.deg2rad(angs)), np.sin(np.deg2rad(angs))
    for rr in np.arange(r_min, 112.0 + 1e-9, 1.0):
        pts_all.extend(zip(tx + rr * ca, ty + rr * sa))
        r_all.extend([rr] * len(angs))
    ok, el, _ = S.feasibility_mask(np.array([p[0] for p in pts_all]),
                                   np.array([p[1] for p in pts_all]), target)
    out = []
    for i, (px, py) in enumerate(pts_all):
        if not ok[i]:
            continue
        if clearance_beyond_buffer(px, py) < MARGIN - 1e-9:
            continue
        if math.degrees(math.atan2(H - S.CAM_H, math.hypot(tx - px, ty - py))) > EL_MAX + 1e-9:
            continue
        out.append((px, py))
    pts = out
    keep_ids = set()
    order_r = sorted(range(len(pts)), key=lambda i: math.hypot(tx - pts[i][0], ty - pts[i][1]))
    keep_ids.update(order_r[:30])
    for p0 in vehicle_positions:
        order = sorted(range(len(pts)), key=lambda i: math.dist(p0, pts[i]))
        keep_ids.update(order[:25])
    for i, (px, py) in enumerate(pts):
        pan = S.wrap360(math.degrees(math.atan2(ty - py, tx - px)))
        if min(abs(S.wrap180(pan)), pan) <= 12.0 and \
                math.hypot(tx - px, ty - py) <= 45.0:
            keep_ids.add(i)
    cl = sorted(keep_ids)
    return dict(pts=[pts[i] for i in cl])


# ----------------------------------------------------------------------
# margin-compliant candidate sets (Q1-Q3 base + per-scenario T4)
# ----------------------------------------------------------------------
CAND_V2 = {a: margin_full_candidates(S.TARGETS[a], S.V_INIT)
           for a in ("T1", "T2", "T3")}
P4.cand = CAND_V2            # reschedule()/greedy read the filtered sets

_orig_build = P4.build_sched_candidates


def build_sched_candidates_m(target, vehicle_positions):
    """problem4.build_sched_candidates with Assumption-7 margins applied."""
    return margin_full_candidates(target, vehicle_positions)


P4.build_sched_candidates = build_sched_candidates_m


# ----------------------------------------------------------------------
# best margin-compliant station per (vehicle state, alert), cached
# ----------------------------------------------------------------------
_station_cache = {}


def best_station(state, alert):
    """(ready_delta, station, pan, tilt): min over margin-compliant candidates
    of travel(state.pos -> station) + slew(state attitude -> target attitude)."""
    key = (round(state["pos"][0], 2), round(state["pos"][1], 2),
           alert, round(state["pan"], 1))
    hit = _station_cache.get(key)
    if hit is not None:
        return hit
    tgt = S.TARGETS[alert]
    bb = None
    for s in CAND_V2[alert]["pts"]:
        tt = PLANNER.time(state["pos"], s)
        pan, tilt = S.pointing_at(s[0], s[1], tgt)
        tau = S.slew_time(state["pan"], state["tilt"], pan, tilt)
        if bb is None or tt + tau < bb[0]:
            bb = (tt + tau, s, pan, tilt)
    _station_cache[key] = bb
    return bb


# ----------------------------------------------------------------------
# P2: rule-based greedy (paper's definition, margin-compliant candidates)
# ----------------------------------------------------------------------
def strategy_P2(states, alive, open_alerts):
    """Earliest-achievable-completion first: each round assigns the remaining
    alert whose best free pair gives the earliest completion (start+5+reset),
    by minimum joint readiness; iterates until no assignment is possible."""
    free = {v: dict(states[v]) for v in alive}
    starts = {}
    remaining = list(open_alerts)
    while remaining:
        best = None
        for a in remaining:
            for (v1, v2) in combinations(sorted(free), 2):
                b1 = best_station(free[v1], a)
                b2 = best_station(free[v2], a)
                start = max(free[v1]["t"] + b1[0], free[v2]["t"] + b2[0])
                reset1 = S.slew_time(b1[2], b1[3], 0.0, 0.0)
                reset2 = S.slew_time(b2[2], b2[3], 0.0, 0.0)
                completion = start + S.OBS_DUR + max(reset1, reset2)
                if best is None or (completion, a) < (best[0], best[1]):
                    best = (completion, a, start, v1, v2, b1, b2, reset1, reset2)
        if best is None:
            break
        completion, a, start, v1, v2, b1, b2, reset1, reset2 = best
        starts[a] = start
        for v, b, rst in ((v1, b1, reset1), (v2, b2, reset2)):
            free[v] = dict(t=start + S.OBS_DUR + rst, pos=b[1], pan=0.0, tilt=0.0)
        remaining.remove(a)
    return starts


# ----------------------------------------------------------------------
# scenario runner: P0 / P1 / P2
# ----------------------------------------------------------------------
def run_v2(name, failed=None, new_alert=None, t_rec=None):
    decision_t = t_rec if t_rec is not None else T_D
    all_alerts = ["T1", "T2", "T3"] + (["T4"] if new_alert is not None else [])
    row = dict(scenario=name, failed=failed, new_alert=new_alert, t_rec=t_rec,
               N=len(all_alerts))

    # ---- P0 (retain; failure-modified stub from problem4.py) ----
    starts_k, stub = P4.simulate_keep(failed=failed)
    if new_alert is not None:
        starts_k = dict(starts_k); starts_k["T4"] = float("inf")
    eta_k, tl_k = P4.evaluate(starts_k, all_alerts)
    row["P0"] = dict(eta=eta_k, t_last=tl_k)

    # ---- shared state capture ----
    alive = [v for v in VEHS if v != failed]
    states, completed, pending = P4.extract_states(stub, decision_t, set(VEHS))
    open_alerts = [a for a in pending if a not in completed]
    if new_alert is not None:
        open_alerts = open_alerts + ["T4"]

    # ---- P1 (event-triggered re-optimization, margin-compliant candidates) ----
    tgt_map = {}
    if new_alert is not None:
        tgt_map["T4"] = S.target_xyz(new_alert)
    saved = dict(S.TARGETS)
    S.TARGETS.update(tgt_map)
    try:
        sim_r, alerts_r, _completed = P4.reschedule(
            decision_t, failed=failed, new_alert=new_alert,
            keep_stub=(starts_k if new_alert is None
                       else {k: v for k, v in starts_k.items() if k != "T4"}, stub))
        starts_r = {a: stub["alert_start"][a] for a in completed}
        starts_r.update({a: sim_r["alert_start"][a] for a in sim_r["alert_start"]})
    finally:
        S.TARGETS.clear(); S.TARGETS.update(saved)
    eta_r, tl_r = P4.evaluate(starts_r, all_alerts)
    row["P1"] = dict(eta=eta_r, t_last=tl_r,
                     starts={a: (round(t, 2) if t < 1e9 else None)
                             for a, t in starts_r.items()})

    # ---- P2 (rule-based greedy on the same captured states) ----
    tgt4 = S.target_xyz(new_alert) if new_alert else None
    saved2 = dict(S.TARGETS)
    if new_alert is not None:
        S.TARGETS["T4"] = tgt4
        CAND_V2["T4"] = margin_full_candidates(
            tgt4, [states[v]["pos"] for v in alive])
    try:
        g = strategy_P2(states, alive, [a for a in open_alerts if a != "T4"]
                        + (["T4"] if new_alert else []))
        starts_2 = {a: stub["alert_start"][a] for a in completed}
        starts_2.update(g)
    except Exception as exc:
        import traceback, json as _json
        with open(OUT / "debug_P2_failure.json", "w", encoding="utf-8") as df:
            _json.dump({"scenario": name, "error": repr(exc),
                        "trace": traceback.format_exc(),
                        "states": {v: {k: (list(x) if isinstance(x, tuple) else x)
                                       for k, x in st.items()}
                                   for v, st in states.items()},
                        "alive": alive, "open": open_alerts,
                        "completed": completed}, df, ensure_ascii=False, indent=2)
        raise
    finally:
        S.TARGETS.clear(); S.TARGETS.update(saved2)
    eta_2, tl_2 = P4.evaluate(starts_2, all_alerts)
    row["P2"] = dict(eta=eta_2, t_last=tl_2,
                     starts={a: (round(t, 2) if t < 1e9 else None)
                             for a, t in starts_2.items()})
    print(f"{name:36s} N={row['N']}  P0 ({eta_k:.2f},{tl_k:6.2f})  "
          f"P1 ({eta_r:.2f},{tl_r:6.2f})  P2 ({eta_2:.2f},{tl_2:6.2f})")
    return row


# ----------------------------------------------------------------------
# 86-scenario family (paper §Q4 Table design)
# ----------------------------------------------------------------------
scenarios = []
for th in range(15, 360, 30):
    for R in (30.0, 40.0, 50.0):
        scenarios.append(dict(name=f"(a) NA θ={th}° R={R:.0f}",
                              new_alert=(float(th), R, 50.0)))
for v in VEHS:
    scenarios.append(dict(name=f"(b) FA {v}", failed=v))
for tr in (10.0, 15.0, 20.0, 25.0, 30.0):
    scenarios.append(dict(name=f"(c) CO t_rec={tr:.0f}", t_rec=tr))
for th in range(15, 360, 30):
    for tr in (10.0, 30.0):
        scenarios.append(dict(name=f"(a+c) NA θ={th}° R=40 t_rec={tr:.0f}",
                              new_alert=(float(th), 40.0, 50.0), t_rec=tr))
for tr in (7.5, 12.5, 15.0, 25.0):
    scenarios.append(dict(name=f"(a+c sweep) NA θ=15° R=40 t_rec={tr:.1f}",
                          new_alert=(15.0, 40.0, 50.0), t_rec=tr))
for v in VEHS:
    scenarios.append(dict(name=f"(b+c) FA {v}+CO(10)", failed=v, t_rec=10.0))
worst_theta = {"V1": 225.0, "V2": 225.0, "V3": 225.0, "V4": 135.0, "V5": 225.0}
for v in VEHS:
    scenarios.append(dict(name=f"(a+b) NA θ={worst_theta[v]:.0f}°+FA {v}",
                          new_alert=(worst_theta[v], 40.0, 50.0), failed=v))
for v in ("V2", "V3"):
    scenarios.append(dict(name=f"(a+b+c) NA θ=225°+FA {v}+CO(10)",
                          new_alert=(225.0, 40.0, 50.0), failed=v, t_rec=10.0))
assert len(scenarios) == 86, f"scenario count {len(scenarios)} != 86"
print(f"=== 86-scenario suite (P0 / P1 / P2) ===")
rows = [run_v2(**sc) for sc in scenarios]

# ----------------------------------------------------------------------
# kinematic sensitivity of the P1 response (3 representative scenarios)
# ----------------------------------------------------------------------
BASE_SPEED, BASE_SLEW = S.V_SPEED, S.SLEW_RATE
SENS = []
for name, kw in (("NA worst (θ=225,R=50)", dict(new_alert=(225.0, 50.0, 50.0))),
                 ("NA easy  (θ=315,R=30)", dict(new_alert=(315.0, 30.0, 50.0))),
                 ("FA V3", dict(failed="V3"))):
    cfgs = {}
    for cname, ss, sl in (("nominal", 1.0, 1.0), ("speed−10%", 0.9, 1.0),
                          ("speed+10%", 1.1, 1.0), ("slew−10%", 1.0, 0.9),
                          ("slew+10%", 1.0, 1.1)):
        S.V_SPEED, S.SLEW_RATE = BASE_SPEED * ss, BASE_SLEW * sl
        starts_s, stub_s = P4.simulate_keep(failed=kw.get("failed"))
        sim_r, _, _ = P4.reschedule(T_D, failed=kw.get("failed"),
                                    new_alert=kw.get("new_alert"),
                                    keep_stub=(starts_s, stub_s))
        cfgs[cname] = round(max(sim_r["alert_start"].values()), 3)
    S.V_SPEED, S.SLEW_RATE = BASE_SPEED, BASE_SLEW
    SENS.append(dict(scenario=name, **cfgs))
    print(f"sens {name}: {cfgs}")

# ----------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------
na_rows = [r for r in rows if r["scenario"].startswith("(a) NA ")]
thetas = sorted({float(r["new_alert"][0]) for r in na_rows})
radii = sorted({float(r["new_alert"][1]) for r in na_rows})

# (1) new-alert grid heat-map
grid = np.full((len(radii), len(thetas)), np.nan)
for r_ in na_rows:
    th, R = float(r_["new_alert"][0]), float(r_["new_alert"][1])
    grid[radii.index(R), thetas.index(th)] = r_["P1"]["t_last"]
fig, ax = plt.subplots(figsize=(8.6, 3.6), dpi=300)
im = ax.imshow(grid, cmap="viridis", aspect="auto",
               extent=[min(thetas) - 15, max(thetas) + 15,
                       min(radii) - 5, max(radii) + 5], origin="lower")
for ri, R in enumerate(radii):
    for gi, th in enumerate(thetas):
        ax.text(th, R, f"{grid[ri, gi]:.1f}", ha="center", va="center",
                fontsize=7.5, color="w")
ax.set_xticks(thetas); ax.set_yticks(radii)
ax.set_xlabel("new-alert azimuth θ (°)"); ax.set_ylabel("new-alert radius R (m)")
ax.set_title("P1 response time over the new-alert grid (η = 1 in all 36 cells)",
             fontsize=10)
fig.colorbar(im, ax=ax, label="t_last (s)")
fig.tight_layout()
fig.savefig(OUT / "fig_q4_v2_grid.png", dpi=300)

# (2) failure bars
fb = [r for r in rows if r["scenario"].startswith("(b) FA ") and "CO" not in r["scenario"]]
names = [r["scenario"].split()[-1] for r in fb]
x = np.arange(len(fb))
fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=300)
ax.bar(x - 0.27, [r["P0"]["t_last"] for r in fb], width=0.27,
       label="P0 retain", color="#c44")
ax.bar(x, [r["P1"]["t_last"] for r in fb], width=0.27,
       label="P1 re-optimization", color="#36c")
ax.bar(x + 0.27, [r["P2"]["t_last"] for r in fb], width=0.27,
       label="P2 rule-based greedy", color="#3a3")
for xi, r in zip(x, fb):
    ax.text(xi - 0.27, r["P0"]["t_last"] + 0.5, f"η={r['P0']['eta']:.2f}",
            ha="center", fontsize=7.5)
    ax.text(xi + 0.27, r["P2"]["t_last"] + 0.5, f"η={r['P2']['eta']:.2f}",
            ha="center", fontsize=7.5)
ax.set_xticks(x); ax.set_xticklabels(names)
ax.set_ylabel("t_last (s)"); ax.set_ylim(0, 27)
ax.set_title("Failure family: three strategies (window [0,60] s)", fontsize=10)
ax.legend(fontsize=8.5); ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "fig_q4_v2_failure.png", dpi=300)

# (3) outage sweep
co = [r for r in rows if r["scenario"].startswith("(c) CO ")]
trs = [r["t_rec"] for r in co]
fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=300)
ax.plot(trs, [r["P0"]["t_last"] for r in co], "s-", color="#c44", label="P0 retain")
ax.plot(trs, [r["P1"]["t_last"] for r in co], "o-", color="#36c",
        label="P1 re-optimization")
ax.plot(trs, [r["P2"]["t_last"] for r in co], "^-", color="#3a3",
        label="P2 rule-based greedy")
ax.set_xlabel("outage recovery time t_rec (s)"); ax.set_ylabel("t_last (s)")
ax.set_title("Pure-outage family: all strategies coincide (zero value boundary)",
             fontsize=10)
ax.legend(fontsize=8.5); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "fig_q4_v2_trec.png", dpi=300)

# (4) value-boundary scatter
fig, ax = plt.subplots(figsize=(7.4, 5.4), dpi=300)
for r in rows:
    d_eta = r["P1"]["eta"] - r["P0"]["eta"]
    d_tl = r["P1"]["t_last"] - r["P0"]["t_last"]
    ax.scatter(d_tl, d_eta, s=26,
               color="#36c" if d_eta > 0 else "#888", alpha=0.85)
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("Δ t_last = t_last(P1) − t_last(P0) (s)")
ax.set_ylabel("Δ η = η(P1) − η(P0)")
ax.set_title("Value of re-optimization across all 86 scenarios:\n"
             "coverage is gained exactly when P0 loses an alert", fontsize=10)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "fig_q4_v2_value.png", dpi=300)

# ----------------------------------------------------------------------
# JSON
# ----------------------------------------------------------------------
res = {
    "baseline": {
        "description": "margin-compliant min-max schedule of Q2(2) (v2 conventions)",
        "assignment": {v: a for v, a in P4.best_asg.items() if a},
        "stations": {f"{v}-{a}": [round(x, 3) for x in s]
                     for (v, a), s in P4.best_st.items()},
        "alert_start_s": {k: round(t, 3) for k, t in
                          P4.base_sim["alert_start"].items()},
    },
    "margins": {"clearance_beyond_buffer_m": MARGIN, "elevation_max_deg": EL_MAX},
    "strategies": ["P0 retain", "P1 event-triggered re-optimization",
                   "P2 rule-based greedy"],
    "scenario_design": "86 = (a) 12x3 new-alert grid + (b) 5 failures + "
                       "(c) 5 outage lengths + (a+c) 28 + (b+c) 5 + (a+b) 5 "
                       "+ (a+b+c) 2",
    "kinematic_sensitivity": SENS,
    "scenarios": rows,
}
with open(OUT / "q4_v2_results.json", "w", encoding="utf-8") as fp:
    json.dump(res, fp, ensure_ascii=False, indent=2)
print("results saved:", OUT / "q4_v2_results.json")
