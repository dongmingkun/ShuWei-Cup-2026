# -*- coding: utf-8 -*-
"""
Problem 4 — Rescheduling under disturbances at t = 5 s, against the Problem 2(2)
optimized baseline.

Disturbances
  NA  new alert T4(theta, R, H=50), two vehicles, theta in [0,360), R in [30,50]
  FA  random failure of exactly one vehicle at t = 5 (equal probability)
  CO  communication outage [5, t_rec], t_rec in (5,30]: vehicles continue ongoing
      tasks; new commands only after recovery.

Rules implemented
  * an alert whose observation started before the disturbance is completed and
    cannot be rescheduled (its vehicles locked until obs_end + reset);
  * dispatched-but-not-observing vehicles may be reassigned (mid-travel: available
    at the interpolated position; mid-slew: available at the station with the
    commanded attitude; waiting: already pointed);
  * neither policy knows the disturbance at t = 0; rescheduling triggers only
    after the disturbance is known (t = 5, or t_rec under an outage).

Policies
  KEEP  execute the original t = 0 schedule regardless.
  RES   event-triggered re-optimization over pending alerts from frozen states,
        using the same enumerator/simulator as Problem 2(2), then an exact
        local search honoring in-progress commitments.

Evaluation window [0,60] s: eta = k/N, t_last = latest start (60 if some alert
cannot start).  Scenario design is justified in the README.

Outputs: results/problem4_results.json, fig_p4_comparison.png
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

VEHS = ["V1", "V2", "V3", "V4", "V5"]
planner = S.PathPlanner()
WINDOW = 60.0
T_D = 5.0


# ----------------------------------------------------------------------
# 0. Baseline schedule (Problem 2(2) optimum) + candidate sets
# ----------------------------------------------------------------------
def build_sched_candidates(target_xyz_, veh_positions):
    full = S.candidate_stations(target_xyz_, ang_step=2.0, rad_step=1.0, r_max=112.0)
    pts = full["pts"]; rr = full["r"]; keep = set()
    for i, p in enumerate(pts):
        if rr[i] <= 14.0:
            keep.add(i)
    for p0 in veh_positions:
        order = sorted(range(len(pts)), key=lambda i: math.dist(p0, pts[i]))
        keep.update(order[:25])
    order = sorted(range(len(pts)), key=lambda i: rr[i])
    keep.update(order[:30])
    for i, p in enumerate(pts):
        tx, ty, _ = target_xyz_
        pan = S.wrap360(math.degrees(math.atan2(ty - p[1], tx - p[0])))
        if min(S.wrap180(pan), abs(S.wrap180(pan))) <= 12.0 and rr[i] <= 45.0:
            keep.add(i)
    cl = sorted(keep)
    return dict(pts=[pts[i] for i in cl], el=full["el"][cl], r=rr[cl])


_c0 = build_sched_candidates(S.TARGETS["T1"], S.V_INIT)
cand = {"T1": _c0,
        "T2": build_sched_candidates(S.TARGETS["T2"], S.V_INIT),
        "T3": build_sched_candidates(S.TARGETS["T3"], S.V_INIT)}

best_asg, best_st, base_sim = S.sched_optimize(
    ["T1", "T2", "T3"], VEHS, planner, cand,
    objective=lambda s: s["max_start"] + 1e-4 * sum(s["alert_start"].values()),
    max_tasks=2, prefilter=150, ls_top=10)
best_asg = {v: a for v, a in best_asg.items() if a}
print("baseline assignment:", best_asg)
print("baseline starts:", {k: round(t, 2) for k, t in base_sim["alert_start"].items()})


# ----------------------------------------------------------------------
# 1. Freeze the planned timeline at a decision time
# ----------------------------------------------------------------------
def extract_states(sim, t_d, alive):
    """States of `alive` vehicles at t_d assuming the plan ran undisturbed to t_d.
    Returns (states, completed, pending).  state: t, pos, pan, tilt,
    commit=(alert, station, planned_ready) for the interrupted task."""
    alert_start = sim["alert_start"]
    completed = [a for a, t0 in alert_start.items() if t0 <= t_d]
    pending = [a for a, t0 in alert_start.items() if t0 > t_d]
    states = {}
    for v, evs in sim["timeline"].items():
        if v not in alive:
            continue
        state = None
        for idx, ev in enumerate(evs):
            k = ev["kind"]
            if k == "depart" and ev["t"] > t_d:
                prev = evs[idx - 1]
                a = ev["alert"]
                arr = next(e for e in evs if e["kind"] == "arrive" and e["alert"] == a)
                rdy = next(e for e in evs if e["kind"] == "ready" and e["alert"] == a)
                prev_pos = prev.get("pos") or next(
                    (e["pos"] for e in reversed(evs[:idx]) if "pos" in e),
                    S.V_INIT[0])
                if prev["t"] >= t_d:      # still at previous station (idle)
                    state = dict(t=t_d, pos=tuple(prev_pos), pan=0.0, tilt=0.0,
                                 commit=None)
                else:                     # mid-travel
                    frac = (t_d - prev["t"]) / max(arr["t"] - prev["t"], 1e-9)
                    pos = (prev_pos[0] + frac * (arr["pos"][0] - prev_pos[0]),
                           prev_pos[1] + frac * (arr["pos"][1] - prev_pos[1]))
                    state = dict(t=t_d, pos=pos, pan=0.0, tilt=0.0,
                                 commit=(a, tuple(arr["pos"]), rdy["t"]))
                break
            if k == "arrive" and ev["t"] <= t_d:
                rdy = next(e for e in evs if e["kind"] == "ready" and e["alert"] == ev["alert"])
                if rdy["t"] > t_d:        # mid-slew
                    state = dict(t=t_d, pos=tuple(ev["pos"]), pan=rdy["pan"],
                                 tilt=rdy["tilt"],
                                 commit=(ev["alert"], tuple(ev["pos"]), rdy["t"]))
                    break
            if k == "ready" and ev["t"] <= t_d:
                a = ev["alert"]
                arr = next(e for e in evs if e["kind"] == "arrive" and e["alert"] == a)
                if alert_start[a] <= t_d:  # observing / resetting (locked)
                    rst = next(e for e in evs if e["kind"] == "reset_done" and e["alert"] == a)
                    if rst["t"] > t_d:
                        state = dict(t=rst["t"], pos=tuple(arr["pos"]), pan=0.0,
                                     tilt=0.0, commit=None)
                        break
                    # reset finished -> free, keep scanning
                else:                      # pointed, waiting for partner
                    state = dict(t=t_d, pos=tuple(arr["pos"]), pan=ev["pan"],
                                 tilt=ev["tilt"],
                                 commit=(a, tuple(arr["pos"]), ev["t"]))
                    break
        if state is None:                 # all tasks done before t_d
            last_pos = evs[-1].get("pos") or next(
                (e["pos"] for e in reversed(evs) if "pos" in e), S.V_INIT[0])
            state = dict(t=max(t_d, evs[-1]["t"]), pos=tuple(last_pos),
                         pan=0.0, tilt=0.0, commit=None)
        states[v] = state
    return states, sorted(completed), sorted(pending)


# ----------------------------------------------------------------------
# 2. Exact simulator with in-progress commitments (for RES)
# ----------------------------------------------------------------------
def simulate_commit(assign, stations, init_state):
    """Like S.simulate but: a vehicle whose FIRST task equals its committed
    (alert, station) becomes ready exactly at the planned time."""
    vehs = [v for v in assign if assign[v]]
    base = {v: dict(t=init_state[v]["t"], pos=tuple(init_state[v]["pos"]),
                    pan=init_state[v]["pan"], tilt=init_state[v]["tilt"])
            for v in vehs}
    ready = {}
    for _ in range(15):
        changed = False
        for v in vehs:
            c = dict(base[v])
            for a in assign[v]:
                s = stations[(v, a)]
                cm = init_state[v].get("commit")
                if cm and cm[0] == a and tuple(cm[1]) == tuple(s):
                    rdy = max(cm[2], c["t"])
                else:
                    t_arr = c["t"] + planner.time(c["pos"], s)
                    pan, tilt = S.pointing_at(s[0], s[1], S.TARGETS[a])
                    rdy = t_arr + S.slew_time(c["pan"], c["tilt"], pan, tilt)
                if ready.get((v, a)) != rdy:
                    ready[(v, a)] = rdy
                    changed = True
                partners = [u for u in vehs if u != v and a in assign[u]]
                start = rdy
                for u in partners:
                    ru = ready.get((u, a))
                    if ru is not None:
                        start = max(start, ru)
                pan_, tilt_ = S.pointing_at(s[0], s[1], S.TARGETS[a])
                reset = S.slew_time(pan_, tilt_, 0.0, 0.0)
                c = dict(t=start + S.OBS_DUR + reset, pos=(s[0], s[1]),
                         pan=0.0, tilt=0.0)
        if not changed:
            break
    alerts = sorted({a for v in vehs for a in assign[v]})
    alert_start = {a: max(ready[(v, a)] for v in vehs if a in assign[v])
                   for a in alerts}
    timeline = {}
    for v in vehs:
        c = dict(base[v])
        ev = [dict(t=c["t"], kind="available", pos=c["pos"])]
        for a in assign[v]:
            s = stations[(v, a)]
            cm = init_state[v].get("commit")
            pan, tilt = S.pointing_at(s[0], s[1], S.TARGETS[a])
            if cm and cm[0] == a and tuple(cm[1]) == tuple(s):
                t_arr = min(c["t"], cm[2]); rdy = max(cm[2], c["t"])
            else:
                t_arr = c["t"] + planner.time(c["pos"], s)
                rdy = t_arr + S.slew_time(c["pan"], c["tilt"], pan, tilt)
            start = alert_start[a]
            end = start + S.OBS_DUR
            reset = S.slew_time(pan, tilt, 0.0, 0.0)
            ev += [dict(t=c["t"], kind="depart", pos=c["pos"], alert=a),
                   dict(t=t_arr, kind="arrive", pos=s, alert=a),
                   dict(t=rdy, kind="ready", pan=pan, tilt=tilt, alert=a),
                   dict(t=start, kind="obs_start", alert=a),
                   dict(t=end, kind="obs_end", alert=a),
                   dict(t=end + reset, kind="reset_done", alert=a)]
            c = dict(t=end + reset, pos=(s[0], s[1]), pan=0.0, tilt=0.0)
        timeline[v] = ev
    return dict(timeline=timeline, alert_start=alert_start, ready=ready,
                max_start=max(alert_start.values()) if alert_start else 0.0)


def exact_local_search(assign, stations, init_state, cand_r, passes=3, topk=8):
    best_st = dict(stations)
    best = simulate_commit(assign, best_st, init_state)
    best_val = best["max_start"]
    slots = [(v, a) for v in assign if assign[v] for a in assign[v]]
    for _ in range(passes):
        improved = False
        for (v, a) in slots:
            cur_s = best_st[(v, a)]
            near = sorted(cand_r[a]["pts"], key=lambda p: math.dist(p, cur_s))[:topk]
            for s in near:
                if s == cur_s:
                    continue
                trial = dict(best_st)
                trial[(v, a)] = s
                sim = simulate_commit(assign, trial, init_state)
                if sim["max_start"] < best_val - 1e-9:
                    best_val, best_st, best, improved = sim["max_start"], trial, sim, True
        if not improved:
            break
    return best_st, best


# ----------------------------------------------------------------------
# 3. Policy KEEP: original plan, optional vehicle loss at t = 5
#    Returns (starts, stub) where stub is the failure-modified event timeline:
#    alerts whose planned obs started <= T_D are completed; a dead vehicle
#    stalls its chains (partners wait forever, pointed, downstream tasks lost).
# ----------------------------------------------------------------------
def simulate_keep(failed=None):
    vehs = [v for v in best_asg if best_asg[v]]
    base_ready = base_sim["ready"]
    completed_early = {a for a, t0 in base_sim["alert_start"].items() if t0 <= T_D}
    ready = {}
    base = {v: dict(t=0.0, pos=tuple(S.V_INIT[int(v[1]) - 1]), pan=0.0, tilt=0.0)
            for v in vehs}
    for _ in range(15):
        changed = False
        for v in vehs:
            c = dict(base[v])
            for a in best_asg[v]:
                if a in completed_early:
                    ready[(v, a)] = base_ready[(v, a)]
                    rst = next(e for e in base_sim["timeline"][v]
                               if e["kind"] == "reset_done" and e["alert"] == a)
                    arr = next(e for e in base_sim["timeline"][v]
                               if e["kind"] == "arrive" and e["alert"] == a)
                    c = dict(t=rst["t"], pos=tuple(arr["pos"]), pan=0.0, tilt=0.0)
                    continue
                if v == failed:
                    if ready.get((v, a)) != float("inf"):
                        ready[(v, a)] = float("inf")
                        changed = True
                    break
                s = best_st[(v, a)]
                t_arr = c["t"] + planner.time(c["pos"], s)
                pan, tilt = S.pointing_at(s[0], s[1], S.TARGETS[a])
                rdy = t_arr + S.slew_time(c["pan"], c["tilt"], pan, tilt)
                if ready.get((v, a)) != rdy:
                    ready[(v, a)] = rdy
                    changed = True
                partners = [u for u in vehs if u != v and a in best_asg[u]]
                start = rdy
                for u in partners:
                    ru = ready.get((u, a))
                    if ru is not None:
                        start = max(start, ru)
                reset = S.slew_time(pan, tilt, 0.0, 0.0)
                c = dict(t=start + S.OBS_DUR + reset, pos=(s[0], s[1]),
                         pan=0.0, tilt=0.0)
        if not changed:
            break
    starts = {}
    for a in ["T1", "T2", "T3"]:
        vs = [v for v in vehs if a in best_asg[v]]
        if not vs:
            starts[a] = float("inf")
            continue
        vals = []
        for v in vs:
            if a in completed_early:
                vals.append(base_ready[(v, a)])
            elif v == failed:
                vals.append(float("inf"))
            else:
                vals.append(ready.get((v, a), float("inf")))
        starts[a] = max(vals)
    # pass 2: failure-modified event stubs
    stub = {}
    for v in vehs:
        ev = [dict(t=0.0, kind="available", pos=tuple(S.V_INIT[int(v[1]) - 1]))]
        c = dict(base[v])
        for a in best_asg[v]:
            if a in completed_early:
                arr = next(e for e in base_sim["timeline"][v]
                           if e["kind"] == "arrive" and e.get("alert") == a)
                for e in base_sim["timeline"][v]:
                    if e["kind"] == "available":
                        continue
                    if e.get("alert") == a:
                        ev.append(dict(e))
                    if e["kind"] == "reset_done" and e.get("alert") == a:
                        c = dict(t=e["t"], pos=tuple(arr["pos"]), pan=0.0, tilt=0.0)
                        break
            elif v == failed:
                break
            else:
                s = best_st[(v, a)]
                pan, tilt = S.pointing_at(s[0], s[1], S.TARGETS[a])
                t_arr = c["t"] + planner.time(c["pos"], s)
                rdy = ready[(v, a)]
                start = starts[a]
                ev += [dict(t=c["t"], kind="depart", pos=c["pos"], alert=a),
                       dict(t=t_arr, kind="arrive", pos=s, alert=a),
                       dict(t=rdy, kind="ready", pan=pan, tilt=tilt, alert=a)]
                if start == float("inf"):
                    c = dict(t=float("inf"), pos=(s[0], s[1]), pan=pan, tilt=tilt)
                    break
                ev += [dict(t=start, kind="obs_start", alert=a),
                       dict(t=start + S.OBS_DUR, kind="obs_end", alert=a),
                       dict(t=start + S.OBS_DUR + S.slew_time(pan, tilt, 0, 0),
                            kind="reset_done", alert=a)]
                c = dict(t=start + S.OBS_DUR + S.slew_time(pan, tilt, 0, 0),
                         pos=(s[0], s[1]), pan=0.0, tilt=0.0)
        stub[v] = ev
    return starts, dict(timeline=stub, alert_start=starts)


# ----------------------------------------------------------------------
# 4. Policy RES: event-triggered re-optimization
# ----------------------------------------------------------------------
def reschedule(decision_t, failed=None, new_alert=None, keep_stub=None):
    """Re-optimize pending alerts from the failure-modified state at decision_t.
    keep_stub: (starts, sim-like dict) from simulate_keep — vehicles executed the
    original plan during any outage; stalled alerts never started."""
    starts_stub, stub_sim = keep_stub
    states, completed, pending = extract_states(stub_sim, decision_t, set(VEHS))
    alive = [v for v in VEHS if v != failed]
    alerts = [a for a in pending if a not in completed]
    cand_r = {a: cand[a] for a in alerts}
    tgt_map = {}
    if new_alert is not None:
        theta, R, H = new_alert
        tgt_map["T4"] = S.target_xyz((theta, R, H))
        cand_r["T4"] = build_sched_candidates(tgt_map["T4"],
                                              [states[v]["pos"] for v in alive])
        alerts = alerts + ["T4"]

    saved = dict(S.TARGETS)
    S.TARGETS.update(tgt_map)
    try:
        init_state = {}
        for v in alive:
            st = states[v]
            d = dict(t=st["t"], pos=st["pos"], pan=st["pan"], tilt=st["tilt"])
            if st["commit"] is not None and st["commit"][0] in alerts:
                d["commit"] = st["commit"]
            init_state[v] = d
        if not alerts:
            return dict(alert_start={}, timeline={}), alerts, completed
        asg, stn, _ = S.sched_optimize(
            alerts, alive, planner, cand_r, init_state=init_state,
            objective=lambda s: s["max_start"] + 1e-4 * sum(s["alert_start"].values()),
            max_tasks=2, prefilter=120, ls_top=8)
        stn, sim = exact_local_search(asg, stn, init_state, cand_r)
        asg = {v: a for v, a in asg.items() if a}
        sim["assignment"] = asg
        sim["stations"] = stn
    finally:
        S.TARGETS.clear()
        S.TARGETS.update(saved)
    return sim, alerts, completed


# ----------------------------------------------------------------------
# 5. Scenario evaluation
# ----------------------------------------------------------------------
def evaluate(starts, all_alerts):
    k = sum(1 for a in all_alerts if starts.get(a, float("inf")) <= WINDOW)
    eta = k / len(all_alerts)
    if k < len(all_alerts):
        t_last = WINDOW
    else:
        t_last = max(starts[a] for a in all_alerts)
    return eta, t_last


def run_scenario(name, failed=None, new_alert=None, t_rec=None):
    row = dict(scenario=name, failed=failed, new_alert=new_alert, t_rec=t_rec)
    all_alerts = ["T1", "T2", "T3"] + (["T4"] if new_alert is not None else [])
    row["N"] = len(all_alerts)
    # KEEP (also produces the failure-modified execution stub)
    starts_k, stub = simulate_keep(failed=failed)
    if new_alert is not None:
        starts_k = dict(starts_k)
        starts_k["T4"] = float("inf")
    eta_k, tl_k = evaluate(starts_k, all_alerts)
    row["keep"] = dict(eta=round(eta_k, 4), t_last=round(tl_k, 2),
                       starts={a: (round(t, 2) if t < 1e9 else None)
                               for a, t in starts_k.items()})
    # RES
    decision_t = t_rec if t_rec is not None else T_D
    sim_r, alerts_r, completed = reschedule(decision_t, failed, new_alert,
                                            keep_stub=(starts_k if new_alert is None
                                                       else {k: v for k, v in starts_k.items()
                                                             if k != "T4"}, stub))
    starts_r = {a: stub["alert_start"][a] for a in completed}
    starts_r.update({a: sim_r["alert_start"][a] for a in sim_r["alert_start"]})
    eta_r, tl_r = evaluate(starts_r, all_alerts)
    row["res"] = dict(eta=round(eta_r, 4), t_last=round(tl_r, 2),
                      starts={a: (round(t, 2) if t < 1e9 else None)
                              for a, t in starts_r.items()},
                      assignment=sim_r.get("assignment", {}),
                      stations={f"{v}-{a}": [round(x, 2) for x in s]
                                for (v, a), s in (sim_r.get("stations") or {}).items()},
                      timeline={v: [dict(t=round(e["t"], 2), kind=e["kind"],
                                         alert=e.get("alert"))
                                    for e in evs]
                                for v, evs in (sim_r.get("timeline") or {}).items()})
    print(f"{name:30s} N={row['N']}  KEEP eta={eta_k:.3f} t={tl_k:6.2f} | "
          f"RES eta={eta_r:.3f} t={tl_r:6.2f}")
    return row


# ----------------------------------------------------------------------
# 6. Scenario set
# ----------------------------------------------------------------------
scenarios = []
for th in range(0, 360, 45):
    for R in (30.0, 50.0):
        scenarios.append(dict(name=f"NA θ={th}° R={R:.0f}", new_alert=(float(th), R, 50.0)))
for v in VEHS:
    scenarios.append(dict(name=f"FA {v}", failed=v))
for tr in (10.0, 20.0, 30.0):
    scenarios.append(dict(name=f"CO t_rec={tr:.0f}", t_rec=tr))
for th in (0.0, 180.0):
    for v in VEHS:
        scenarios.append(dict(name=f"NA θ={th:.0f}°+FA {v}",
                              new_alert=(th, 40.0, 50.0), failed=v))
for v in VEHS:
    scenarios.append(dict(name=f"FA {v}+CO(20)", failed=v, t_rec=20.0))
for th in (0.0, 90.0, 180.0, 270.0):
    scenarios.append(dict(name=f"NA θ={th:.0f}°+CO(20)",
                          new_alert=(th, 40.0, 50.0), t_rec=20.0))
scenarios.append(dict(name="NA θ=0°+FA V3+CO(20)",
                      new_alert=(0.0, 40.0, 50.0), failed="V3", t_rec=20.0))
scenarios.append(dict(name="NA θ=180°+FA V3+CO(20)",
                      new_alert=(180.0, 40.0, 50.0), failed="V3", t_rec=20.0))

print(f"\n=== {len(scenarios)} scenarios ===")
results = [run_scenario(**sc) for sc in scenarios]

# ----------------------------------------------------------------------
# 7. Figures + save
# ----------------------------------------------------------------------
fig, axs = plt.subplots(2, 1, figsize=(13.5, 8.5), dpi=300, sharex=True)
x = np.arange(len(results))
names = [r["scenario"] for r in results]
axs[0].bar(x - 0.2, [r["keep"]["eta"] for r in results], width=0.4,
           label="KEEP original schedule", color="#ef5350")
axs[0].bar(x + 0.2, [r["res"]["eta"] for r in results], width=0.4,
           label="RES event-triggered rescheduling", color="#42a5f5")
axs[0].set_ylabel("completion ratio η = k/N")
axs[0].set_ylim(0, 1.15)
axs[0].legend(fontsize=9, loc="lower left")
axs[0].set_title("Problem 4 — disturbance scenarios, evaluation window [0, 60] s")
axs[0].grid(axis="y", alpha=0.3)
axs[1].bar(x - 0.2, [r["keep"]["t_last"] for r in results], width=0.4, color="#ef5350")
axs[1].bar(x + 0.2, [r["res"]["t_last"] for r in results], width=0.4, color="#42a5f5")
axs[1].set_ylabel("t_last (s)")
axs[1].set_xticks(x)
axs[1].set_xticklabels(names, rotation=62, ha="right", fontsize=7)
axs[1].grid(axis="y", alpha=0.3)
axs[1].set_xlabel("scenario")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_p4_comparison.png", dpi=300)

with open(f"{OUT}/problem4_results.json", "w", encoding="utf-8") as fp:
    json.dump(dict(baseline=dict(assignment=best_asg,
                                 stations={f"{v}-{a}": [round(x, 3) for x in s]
                                          for (v, a), s in best_st.items()},
                                 starts={k: round(t, 3) for k, t in
                                         base_sim["alert_start"].items()}),
                     scenarios=results),
              fp, ensure_ascii=False, indent=2)
print("saved results + figure.")
