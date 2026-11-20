# -*- coding: utf-8 -*-
"""
Problem 2(2) — Vehicle reuse and routing under three concurrent alerts.

Model
-----
* 6 vehicle-task slots (2 vehicles x 3 alerts), 5 vehicles -> exactly one vehicle
  serves two alerts (chain), the others one each (we also allow a 3-task chain
  and an idle vehicle in the enumeration for completeness).
* Each (vehicle, alert) task needs a station inside the alert's observable region;
  the vehicle timeline obeys the full timing sequence:
  travel -> slew (30 deg/s, slower axis governs) -> wait -> 5 s joint observation
  -> reset to due-east -> depart.
* Objective: minimize the maximum joint-observation start time over the three
  alerts.  Solution method: exhaustive assignment enumeration (pruned) +
  travel-lower-bound pre-ranking + exact fixed-point timeline simulation +
  coordinate-descent local search on station choices.
* Baseline: a "nearest-assignment heuristic" processed in alert label order
  (each alert takes the two vehicles with the smallest optimistic earliest-ready
  time from their committed position). Same constraints, no lookahead.

Outputs: results/problem2_2_results.json, fig_p2_2_gantt.png
"""
import json
import math
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import swb_common as S

OUT = "results"
os.makedirs(OUT, exist_ok=True)

ALERT_NAMES = ["T1", "T2", "T3"]
VEHS = ["V1", "V2", "V3", "V4", "V5"]
planner = S.PathPlanner()

# ----------------------------------------------------------------------
# 1. Candidate stations per alert (scheduler subset)
# ----------------------------------------------------------------------
def build_sched_candidates(alert_names):
    cand = {}
    for a in alert_names:
        # r_max must cover the whole field: some vehicles can observe from their
        # standing positions far from the target (e.g. V3 sees T3 from r~86 m)
        full = S.candidate_stations(S.TARGETS[a], ang_step=2.0, rad_step=1.0, r_max=112.0)
        pts = full["pts"]
        rr = full["r"]
        keep = set()
        # (a) all close-in candidates r <= 14
        for i, p in enumerate(pts):
            if rr[i] <= 14.0:
                keep.add(i)
        # (b) per-vehicle 25 nearest feasible candidates (approach directions)
        for vi, v in enumerate(VEHS):
            order = sorted(range(len(pts)), key=lambda i: math.dist(S.V_INIT[vi], pts[i]))
            keep.update(order[:25])
        # (c) 30 closest-to-target overall
        order = sorted(range(len(pts)), key=lambda i: rr[i])
        keep.update(order[:30])
        # (d) due-east-bearing options (pan ~ 0 -> minimal slew), r <= 45
        for i, p in enumerate(pts):
            tx, ty, H = S.TARGETS[a]
            pan = S.wrap360(math.degrees(math.atan2(ty - p[1], tx - p[0])))
            if min(S.wrap180(pan), abs(S.wrap180(pan))) <= 12.0 and rr[i] <= 45.0:
                keep.add(i)
        cl = sorted(keep)
        cand[a] = dict(pts=[pts[i] for i in cl],
                       el=full["el"][cl], r=rr[cl])
        print(f"{a}: {len(pts)} feasible -> {len(cand[a]['pts'])} scheduler candidates")
    return cand


cand = build_sched_candidates(ALERT_NAMES)

# ----------------------------------------------------------------------
# 2. Optimized scheme
# ----------------------------------------------------------------------
def objective(sim):
    return sim["max_start"] + 1e-4 * sum(sim["alert_start"].values())


best_asg, best_st, best_sim = S.sched_optimize(
    ALERT_NAMES, VEHS, planner, cand, objective=objective,
    max_tasks=2, prefilter=150, ls_top=10)

print("\n=== OPTIMIZED SCHEME ===")
print("assignment:", {v: a for v, a in best_asg.items() if a})
print("stations:", {f"{v}-{a}": tuple(round(x, 2) for x in s)
                   for (v, a), s in best_st.items()})
print("alert obs starts:", {k: round(t, 2) for k, t in best_sim["alert_start"].items()})
print("max obs start: %.2f s" % best_sim["max_start"])

# ----------------------------------------------------------------------
# 3. Nearest-assignment heuristic baseline
# ----------------------------------------------------------------------
h_asg, h_st = S.nearest_heuristic(ALERT_NAMES, VEHS, planner, cand,
                                  alert_order=["T1", "T2", "T3"])
h_sim = S.simulate(h_asg, h_st, planner)
print("\n=== NEAREST-ASSIGNMENT HEURISTIC ===")
print("assignment:", {v: a for v, a in h_asg.items() if a})
print("alert obs starts:", {k: round(t, 2) for k, t in h_sim["alert_start"].items()})
print("max obs start: %.2f s" % h_sim["max_start"])
print("improvement of optimized over heuristic: %.2f s (%.1f%%)" % (
    h_sim["max_start"] - best_sim["max_start"],
    100 * (h_sim["max_start"] - best_sim["max_start"]) / h_sim["max_start"]))

# ----------------------------------------------------------------------
# 4. Gantt chart
# ----------------------------------------------------------------------
COLORS = dict(travel="#9e9e9e", slew="#ff9800", wait="#fff3c4", obs="#4caf50",
              reset="#29b6f6")

def draw_gantt(sim, title, fname):
    fig, ax = plt.subplots(figsize=(11, 4.6), dpi=300)
    active = [v for v in VEHS if sim["timeline"].get(v)]
    ylabels = []
    for yi, v in enumerate(active):
        ylabels.append(v)
        t_prev = None
        for ev in sim["timeline"][v]:
            if ev["kind"] == "depart":
                t0 = ev["t"]; pos0 = ev["pos"]; alert = ev["alert"]
            elif ev["kind"] == "arrive":
                ax.add_patch(plt.Rectangle((t0, yi - 0.33), ev["t"] - t0, 0.66,
                                           color=COLORS["travel"], alpha=0.8))
                slew_end = next(e["t"] for e in sim["timeline"][v]
                                if e["kind"] == "ready" and e["alert"] == alert)
                ax.add_patch(plt.Rectangle((ev["t"], yi - 0.33), slew_end - ev["t"], 0.66,
                                           color=COLORS["slew"], alpha=0.9))
                obs_start = sim["alert_start"][alert]
                obs_end = obs_start + S.OBS_DUR
                reset_end = next(e["t"] for e in sim["timeline"][v]
                                 if e["kind"] == "reset_done" and e["alert"] == alert)
                ax.add_patch(plt.Rectangle((slew_end, yi - 0.33),
                                           obs_start - slew_end, 0.66,
                                           color=COLORS["wait"], alpha=0.9))
                ax.add_patch(plt.Rectangle((obs_start, yi - 0.33), S.OBS_DUR, 0.66,
                                           color=COLORS["obs"], alpha=0.95))
                ax.add_patch(plt.Rectangle((obs_end, yi - 0.33),
                                           reset_end - obs_end, 0.66,
                                           color=COLORS["reset"], alpha=0.9))
                ax.text((t0 + ev["t"]) / 2, yi, alert, ha="center", va="center",
                        fontsize=7, color="w", weight="bold")
    for a, t in sim["alert_start"].items():
        ax.axvline(t, color="red", ls=":", lw=1.0, alpha=0.8)
        ax.text(t, len(active) - 0.35, f"{a} @{t:.1f}s", color="red", fontsize=8)
    ax.set_yticks(range(len(active)))
    ax.set_yticklabels(ylabels)
    ax.invert_yaxis()
    ax.set_xlabel("time (s)")
    ax.set_xlim(0, max(60, sim["makespan"] + 2))
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.3)
    from matplotlib.patches import Patch
    handles = [Patch(color=COLORS["travel"], label="travel"),
               Patch(color=COLORS["slew"], label="gimbal slew"),
               Patch(color=COLORS["wait"], label="hold (wait partner)"),
               Patch(color=COLORS["obs"], label="joint observation (5 s)"),
               Patch(color=COLORS["reset"], label="gimbal reset")]
    ax.legend(handles=handles, loc="upper right", fontsize=8, ncols=5)
    fig.tight_layout()
    fig.savefig(fname, dpi=300)
    print("figure saved:", fname)


draw_gantt(best_sim, f"Optimized scheme — max observation start = "
                     f"{best_sim['max_start']:.1f} s",
           f"{OUT}/fig_p2_2_gantt.png")
draw_gantt(h_sim, f"Nearest-assignment heuristic — max observation start = "
                  f"{h_sim['max_start']:.1f} s",
           f"{OUT}/fig_p2_2_gantt_heuristic.png")

# ----------------------------------------------------------------------
# 5. Save results
# ----------------------------------------------------------------------
def timeline_json(sim):
    out = {}
    for v, evs in sim["timeline"].items():
        rows = []
        for ev in evs:
            row = dict(t=round(ev["t"], 3), kind=ev["kind"])
            if "alert" in ev:
                row["alert"] = ev["alert"]
            if "pos" in ev:
                row["pos"] = [round(ev["pos"][0], 2), round(ev["pos"][1], 2)]
            if "pan" in ev:
                row["pan_deg"] = round(ev["pan"], 2); row["tilt_deg"] = round(ev["tilt"], 2)
            rows.append(row)
        out[v] = rows
    return out


res = {
    "optimized": {
        "assignment": {v: a for v, a in best_asg.items() if a},
        "stations": {f"{v}-{a}": [round(x, 3) for x in s] for (v, a), s in best_st.items()},
        "alert_start_s": {k: round(t, 3) for k, t in best_sim["alert_start"].items()},
        "max_start_s": round(best_sim["max_start"], 3),
        "makespan_s": round(best_sim["makespan"], 3),
        "timeline": timeline_json(best_sim),
    },
    "heuristic": {
        "rule": "alerts in label order; each takes the two vehicles with smallest "
                "optimistic earliest-ready time; tasks appended to vehicle chains",
        "assignment": {v: a for v, a in h_asg.items() if a},
        "stations": {f"{v}-{a}": [round(x, 3) for x in s] for (v, a), s in h_st.items()},
        "alert_start_s": {k: round(t, 3) for k, t in h_sim["alert_start"].items()},
        "max_start_s": round(h_sim["max_start"], 3),
        "makespan_s": round(h_sim["makespan"], 3),
        "timeline": timeline_json(h_sim),
    },
    "improvement_s": round(h_sim["max_start"] - best_sim["max_start"], 3),
}
with open(f"{OUT}/problem2_2_results.json", "w", encoding="utf-8") as fp:
    json.dump(res, fp, ensure_ascii=False, indent=2)
print("results saved:", f"{OUT}/problem2_2_results.json")
