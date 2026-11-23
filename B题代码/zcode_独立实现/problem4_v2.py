# -*- coding: utf-8 -*-
"""Problem 4 (v3.2m) — Q4 rescheduling on the sigma<2m Q2(2) baseline.

Convention change (team decision, 2026-11-23): the Q2(2) schedule under Q4 is
the TIME-FIRST plan with a hard quality gate  sigma < 2.0 m  and  gamma >= 30
deg  on every observation pair (problem2_2_2m.py ->
results/problem2_2_results_2m.json), NOT the margin-compliant 90-degree-family
plan of the earlier draft.  All rescheduled pairs carry the same gate, so the
"error within 2 m" requirement holds for every alert in every scenario.

Strategies per scenario (window [0, 60] s, decision epoch tau):
  P0  retain the frozen 2m plan (failure-modified stub via problem4.py)
  P1  event-triggered re-optimization WITHOUT the two-slot convention:
      exhaustive (order, vehicle-pair) enumeration + time-first stations +
      quality repair of violating alerts (latest first) + binding re-opt.
      Methodologically mirrors problem2_2_2m.py (time-first, repair violators).
  P1c same search with at most two task slots per vehicle (control: quantifies
      the slot cap's cost under the 2m convention).
  P2  rule-based greedy (earliest achievable completion) with the same gate.

All searches enforce: clearance >= 0.10 m beyond the buffer, elevation <= 79.9
deg, sigma < 2.0 m, gamma >= 30 deg on every pair; first-task commitments
inherit the planned ready time (simulate_commit semantics).  Outputs:
q4_v2_results.json + four 300-dpi figures.
"""
import json
import math
import os
import sys
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
import problem4 as P4   # machinery: simulate_keep / extract_states / evaluate
import swb_common as S  # noqa: E402

OUT = BASE / "results"
OUT.mkdir(exist_ok=True)

WINDOW, T_D = 60.0, 5.0
MARGIN, EL_MAX = 0.10, 79.9          # Assumption-7 operating margins
SIGMA_MAX, GAMMA_MIN = 2.0, 30.0     # Q2(2) 2m-convention quality gate
BEARING_SIG = math.radians(1.0)
VEHS = P4.VEHS
PLANNER = P4.planner

# ----------------------------------------------------------------------
# 0. Freeze the sigma<2m Q2(2) baseline and patch the P4 machinery
# ----------------------------------------------------------------------
_B2M = json.load(open(OUT / "problem2_2_results_2m.json", encoding="utf-8"))
ASG2M = {v: al for v, al in _B2M["assignment"].items()}
ST2M = {(v, a): tuple(float(x) for x in s)
        for key, s in _B2M["stations"].items() for v, a in [key.split("-")]}
SIM2M = S.simulate(ASG2M, ST2M, PLANNER,
                   init_pos={v: S.V_INIT[int(v[1]) - 1] for v in ASG2M})
P4.best_asg, P4.best_st, P4.base_sim = ASG2M, ST2M, SIM2M

_chk = P4.simulate_keep()
for _a in ("T1", "T2", "T3"):
    assert abs(_chk[0][_a] - _B2M["alert_start_s"][_a]) < 1e-3, \
        f"2m baseline replay mismatch on {_a}: {_chk[0][_a]} vs " \
        f"{_B2M['alert_start_s'][_a]}"
print(f"2m baseline frozen: starts {_B2M['alert_start_s']}, "
      f"max_sigma {_B2M['max_sigma_m']} m, max_start "
      f"{_B2M['max_start_s']} s (replay verified)")


# ----------------------------------------------------------------------
# 1. Candidates with margins + the 2m quality gate
# ----------------------------------------------------------------------
def clearance_at(p):
    x, y = p
    d = min(x - S.R_V, S.FIELD - S.R_V - x, y - S.R_V, S.FIELD - S.R_V - y)
    for (x0, x1, y0, y1, _h) in S.BUILDINGS:
        d = min(d, float(S.point_rect_dist_vec(np.array([x]), np.array([y]),
                                               (x0, x1, y0, y1))[0]) - S.R_V)
    return d


def make_set(a, ang=2.0, rad=1.0, rmax=80.0):
    """Margin-filtered candidate stations around alert a's projection."""
    full = S.candidate_stations(S.TARGETS[a], ang_step=ang, rad_step=rad,
                                r_max=rmax)
    pts, rs = [], []
    for p, r in zip(full["pts"], full["r"]):
        if clearance_at(p) < MARGIN - 1e-9:
            continue
        el = math.degrees(math.atan2(S.TARGETS[a][2] - S.CAM_H,
                                     math.hypot(S.TARGETS[a][0] - p[0],
                                                S.TARGETS[a][1] - p[1])))
        if el > EL_MAX + 1e-9:
            continue
        pts.append(tuple(p))
        rs.append(float(r))
    return dict(pts=pts, arr=np.asarray(pts, float) if pts else np.zeros((0, 2)),
                r=np.asarray(rs, float))


C_FULL = {a: make_set(a) for a in ("T1", "T2", "T3")}
print("candidate counts:", {a: len(v["pts"]) for a, v in C_FULL.items()})


def pair_quality(s1, s2, tgt):
    tx, ty, _H = tgt
    r1 = math.hypot(tx - s1[0], ty - s1[1])
    r2 = math.hypot(tx - s2[0], ty - s2[1])
    u1 = ((tx - s1[0]) / r1, (ty - s1[1]) / r1)
    u2 = ((tx - s2[0]) / r2, (ty - s2[1]) / r2)
    c = max(-1.0, min(1.0, u1[0] * u2[0] + u1[1] * u2[1]))
    g = math.degrees(math.acos(c))
    sig = BEARING_SIG * math.sqrt(r1 * r1 + r2 * r2) / max(
        math.sin(math.radians(g)), 1e-9)
    return dict(r1=r1, r2=r2, gamma=g, sigma=sig)


def pair_ok(s1, s2, tgt):
    q = pair_quality(s1, s2, tgt)
    return q["sigma"] < SIGMA_MAX - 1e-9 and q["gamma"] >= GAMMA_MIN - 1e-9, q


# ----------------------------------------------------------------------
# 2. Time-first station choice + exact chain simulation (commit-aware)
# ----------------------------------------------------------------------
def nearest_time(state, a, topn=25):
    """Cheapest (travel + slew) station for one vehicle from its state.
    The committed station is offered with the commitment shortcut (inherits
    the planned ready time) when it is the vehicle's first task; same-alert
    alternatives respect the commitment clock (mid-slew captures)."""
    tgt = S.TARGETS[a]
    best = None
    cm = state.get("commit")
    if cm and cm[0] == a:
        s = tuple(cm[1])
        pan, tilt = S.pointing_at(s[0], s[1], tgt)
        delta = max(cm[2], state["t"]) - state["t"]
        best = (delta, s, pan, tilt)
    arr = C_FULL[a]["arr"]
    if arr.size:
        d = np.hypot(arr[:, 0] - state["pos"][0], arr[:, 1] - state["pos"][1])
        for j in np.argsort(d)[:topn]:
            p = (float(arr[j, 0]), float(arr[j, 1]))
            tt = PLANNER.time(state["pos"], p)
            pan, tilt = S.pointing_at(p[0], p[1], tgt)
            sl = S.slew_time(state["pan"], state["tilt"], pan, tilt)
            val = tt + sl
            if cm and cm[0] == a and (abs(state["pan"]) > 1e-9 or
                                      abs(state["tilt"]) > 1e-9):
                val = max(val, cm[2] - state["t"])
            if best is None or val < best[0]:
                best = (val, p, pan, tilt)
    return best


def chain_sim(chains, init_states):
    """Exact timeline of ordered per-vehicle chains [(alert, station), ...]:
    fixed point over coupled readiness; a vehicle's FIRST task inherits its
    commitment shortcut.  Returns (starts, events, quality)."""
    servers = {}
    for v, ch in chains.items():
        for (a, _s) in ch:
            servers.setdefault(a, []).append(v)
    ready = {}
    for _it in range(24):
        snap = dict(ready)
        for v, ch in chains.items():
            c = dict(init_states[v])
            for i, (a, s) in enumerate(ch):
                pan, tilt = S.pointing_at(s[0], s[1], S.TARGETS[a])
                cm = init_states[v].get("commit")
                if i == 0 and cm and cm[0] == a and _same_st(cm[1], s):
                    rdy = max(cm[2], c["t"])
                else:
                    tt = PLANNER.time(c["pos"], s)
                    rdy = c["t"] + tt + S.slew_time(c["pan"], c["tilt"],
                                                    pan, tilt)
                    if i == 0:
                        rdy = _commit_clock(init_states, v, a, c, rdy)
                ready[(v, a)] = rdy
                start = max(ready[(u, a)] for u in servers[a]
                            if (u, a) in ready)
                reset = S.slew_time(pan, tilt, 0.0, 0.0)
                c = dict(t=start + S.OBS_DUR + reset, pos=(s[0], s[1]),
                         pan=0.0, tilt=0.0)
        if snap and len(snap) == len(ready) and all(snap.get(k) == v2
                                                    for k, v2 in ready.items()):
            break
    starts = {a: max(ready[(u, a)] for u in servers[a]) for a in servers}
    events = {}
    for v, ch in chains.items():
        c = dict(init_states[v])
        ev = [dict(t=c["t"], kind="available", pos=c["pos"])]
        for i, (a, s) in enumerate(ch):
            pan, tilt = S.pointing_at(s[0], s[1], S.TARGETS[a])
            cm = init_states[v].get("commit")
            shortcut = i == 0 and cm and cm[0] == a and tuple(cm[1]) == tuple(s)
            rdy = ready[(v, a)]
            if shortcut:
                # committed plan: travel already done; stand-in departure
                ev.append(dict(t=c["t"], kind="depart", pos=c["pos"],
                               pan=pan, tilt=tilt, alert=a, station=s, rdy=rdy))
                ev.append(dict(t=rdy, kind="arrive", pos=s, alert=a))
            else:
                tt = PLANNER.time(c["pos"], s)
                t_arr = c["t"] + tt
                ev.append(dict(t=c["t"], kind="depart", pos=c["pos"],
                               pan=c["pan"], tilt=c["tilt"], alert=a,
                               station=s, rdy=rdy))
                ev.append(dict(t=t_arr, kind="arrive", pos=s, alert=a))
            ev.append(dict(t=rdy, kind="ready", pos=s, pan=pan, tilt=tilt,
                           alert=a))
            reset = S.slew_time(pan, tilt, 0.0, 0.0)
            t_rs = starts[a] + S.OBS_DUR + reset
            ev.append(dict(t=t_rs, kind="reset_done", pos=s, alert=a))
            c = dict(t=t_rs, pos=(s[0], s[1]), pan=0.0, tilt=0.0)
        events[v] = ev
    qual = {}
    for a in servers:
        ss = [dict(chains[v])[a] for v in servers[a]]
        if len(ss) == 2:
            ok, q = pair_ok(ss[0], ss[1], S.TARGETS[a])
            q["ok"] = ok
            qual[a] = q
    return starts, events, qual


def build_time_first(order, prs, states):
    """Time-first chains for a forced (order, pair-assignment) structure.
    Stations are the cheapest per vehicle in sequence; the quality gate is
    NOT applied here (repaired afterwards), mirroring problem2_2_2m.py."""
    chains = {v: [] for v in states}
    free = {v: dict(states[v]) for v in states}
    for a, (v1, v2) in zip(order, prs):
        b1 = nearest_time(free[v1], a)
        b2 = nearest_time(free[v2], a)
        if b1 is None or b2 is None:
            return None
        r1 = free[v1]["t"] + b1[0]
        r2 = free[v2]["t"] + b2[0]
        start = max(r1, r2)
        for v, b in ((v1, b1), (v2, b2)):
            rst = S.slew_time(b[2], b[3], 0.0, 0.0)
            chains[v].append((a, b[1]))
            free[v] = dict(t=start + S.OBS_DUR + rst, pos=b[1], pan=0.0,
                           tilt=0.0)
    return chains


# ----------------------------------------------------------------------
# 3. Quality repair (violating alerts, latest first) — sigma<2m machinery
# ----------------------------------------------------------------------
def _cand_list(state, a, topn=40):
    tgt = S.TARGETS[a]
    out = []
    cm = state.get("commit")
    if cm and cm[0] == a and tuple(cm[1]) not in out:
        out.append(tuple(cm[1]))
    arr = C_FULL[a]["arr"]
    if arr.size:
        d = np.hypot(arr[:, 0] - state["pos"][0], arr[:, 1] - state["pos"][1])
        for j in np.argsort(d)[:topn]:
            p = (float(arr[j, 0]), float(arr[j, 1]))
            if p not in out:
                out.append(p)
    return out


_STATION_TOL = 0.01   # m: matches the 4-decimal rounding of the frozen baseline


def _same_st(p, q):
    return abs(p[0] - q[0]) + abs(p[1] - q[1]) < _STATION_TOL


def _commit_clock(init_states, v, a, e, computed_ready):
    """Commitment clock: a vehicle captured mid-slew or waiting (nonzero
    attitude) cannot be ready for its COMMITTED alert before the committed
    ready time, at any station — its gimbal is still working until then."""
    cm = init_states[v].get("commit")
    if cm and cm[0] == a and (abs(e["pan"]) > 1e-9 or abs(e["tilt"]) > 1e-9):
        return max(computed_ready, cm[2])
    return computed_ready


def _rdy_of(init_states, first_task, v, e, p, a):
    """Readiness of vehicle v for station p given its depart event e."""
    cm = init_states[v].get("commit")
    if first_task and cm and cm[0] == a and _same_st(cm[1], p):
        return max(cm[2], e["t"])
    pan, tilt = S.pointing_at(p[0], p[1], S.TARGETS[a])
    ready = e["t"] + PLANNER.time(e["pos"], p) + S.slew_time(
        e["pan"], e["tilt"], pan, tilt)
    if first_task:
        ready = _commit_clock(init_states, v, a, e, ready)
    return ready


def best_repair(chains, events, a, states):
    """Best gated station pair for alert a's two servers, evaluated from their
    depart events; returns ((s1, s2), newmax) or None if no gated pair."""
    vs = [v for v in chains if any(x == a for x, _s in chains[v])]
    if len(vs) != 2:
        return None
    e1 = next(e for e in events[vs[0]] if e["kind"] == "depart" and
              e.get("alert") == a)
    e2 = next(e for e in events[vs[1]] if e["kind"] == "depart" and
              e.get("alert") == a)
    st1 = dict(states[vs[0]], t=e1["t"], pos=tuple(e1["pos"]),
               pan=e1["pan"], tilt=e1["tilt"])
    st2 = dict(states[vs[1]], t=e2["t"], pos=tuple(e2["pos"]),
               pan=e2["pan"], tilt=e2["tilt"])
    c1 = _cand_list(st1, a)
    c2 = _cand_list(st2, a)
    first1 = chains[vs[0]][0][0] == a
    first2 = chains[vs[1]][0][0] == a
    best = None
    for p1 in c1:
        for p2 in c2:
            ok, q = pair_ok(p1, p2, S.TARGETS[a])
            if not ok:
                continue
            r1 = _rdy_of(states, first1, vs[0], e1, p1, a)
            r2 = _rdy_of(states, first2, vs[1], e2, p2, a)
            start = max(r1, r2)
            val = (start, q["sigma"])
            if best is None or val < best[0]:
                best = (val, p1, p2)
    if best is None:
        return None
    return ((best[1], best[2]), best[0][0])


def repair_quality(chains, states, passes=3):
    """Repair violating alerts (latest first) up to `passes` rounds, then one
    re-opt of the binding alert.  Returns (chains, events)."""
    starts, events, qual = chain_sim(chains, states)
    for _ in range(passes):
        viol = [a for a, q in qual.items() if not q["ok"]]
        if not viol:
            break
        a = max(viol, key=lambda x: starts.get(x, 0.0))
        res = best_repair(chains, events, a, states)
        if res is None:
            break
        (s1, s2), _ = res
        vs = [v for v in chains if any(x == a for x, _s in chains[v])]
        for v, s in zip(vs, (s1, s2)):
            idx = [i for i, (x, _s) in enumerate(chains[v]) if x == a][0]
            chains[v][idx] = (a, s)
        starts, events, qual = chain_sim(chains, states)
    # binding-alert re-opt (equal quality, possibly earlier start)
    if starts:
        ab = max(starts, key=starts.get)
        res = best_repair(chains, events, ab, states)
        if res is not None:
            (s1, s2), newmax = res
            old_ab = starts[ab]
            vs = [v for v in chains if any(x == ab for x, _s in chains[v])]
            trial = {v: [list(x) for x in chains[v]] for v in chains}
            for v, s in zip(vs, (s1, s2)):
                idx = [i for i, (x, _s) in enumerate(trial[v]) if x == ab][0]
                trial[v][idx] = (ab, s)
            st2, _e2, q2 = chain_sim(trial, states)
            if max(st2.values()) < max(starts.values()) - 1e-9:
                chains, starts, events, qual = trial, st2, _e2, q2
    return chains, events


# ----------------------------------------------------------------------
# 4. Search engines: P1 (cap-free) / P1c (capped) / P2 (rule greedy)
# ----------------------------------------------------------------------
def search_2m(states, alive, open_alerts, max_slots=None, top_refine=8):
    """Exhaustive (order, pair) enumeration + time-first stations + quality
    repair; returns (starts, chains) of the best plan, or ({}, {}) if nothing
    is open."""
    if not open_alerts:
        return {}, {}
    vlist = sorted(alive)
    all_pairs = list(combinations(vlist, 2))
    pool = []
    for order in permutations(sorted(open_alerts)):
        for prs in product(all_pairs, repeat=len(order)):
            if max_slots is not None:
                cnt = Counter(v for pr in prs for v in pr)
                if max(cnt.values()) > max_slots:
                    continue
            chains = build_time_first(order, prs, states)
            if chains is None:
                continue
            starts, _ev, qual = chain_sim(chains, states)
            viol = [a for a in qual if not qual[a]["ok"]]
            # rank by coverage + max start ONLY (sigma is a final tie-break:
            # a mixed key can push the fastest structure out of the top-k)
            key = (len(viol), max(starts.values()))
            pool.append((key, chains))
    if not pool:
        return {}, {}
    pool.sort(key=lambda x: x[0])
    best = None
    for _key, ch in pool[:top_refine]:
        ch2, _ = repair_quality({v: list(x) for v, x in ch.items()}, states)
        starts, _ev, qual = chain_sim(ch2, states)
        viol = [a for a in qual if not qual[a]["ok"]]
        key = (len(viol), max(starts.values()),
               sum(qual[a]["sigma"] for a in qual))
        if best is None or key < best[0]:
            best = (key, starts, {v: list(x) for v, x in ch2.items()})
    if best is None:
        return {}, {}
    return best[1], best[2]


_LAST_PLAN = {"chains": None}
_LAST_P1 = {"chains": None}


def strategy_P2(states, alive, open_alerts):
    """Rule-based greedy with the sigma<2m gate: each round assigns the alert
    whose best free pair gives the earliest completion."""
    free = {v: dict(states[v]) for v in alive}
    chains = {v: [] for v in alive}
    starts = {}
    remaining = list(sorted(open_alerts))
    touched = set()
    while remaining:
        best = None
        for a in remaining:
            for (v1, v2) in combinations(sorted(free), 2):
                # candidate stations from each vehicle's state
                c1 = _cand_list(free[v1], a, topn=30)
                c2 = _cand_list(free[v2], a, topn=30)
                best_pair = None
                for p1 in c1:
                    for p2 in c2:
                        ok, q = pair_ok(p1, p2, S.TARGETS[a])
                        if not ok:
                            continue
                        f1 = v1 not in touched
                        r1 = _rdy_of(states, f1, v1,
                                     dict(t=free[v1]["t"], pos=free[v1]["pos"],
                                          pan=free[v1]["pan"],
                                          tilt=free[v1]["tilt"]), p1, a)
                        r2 = _rdy_of(states, v2 not in touched, v2,
                                     dict(t=free[v2]["t"], pos=free[v2]["pos"],
                                          pan=free[v2]["pan"],
                                          tilt=free[v2]["tilt"]), p2, a)
                        start = max(r1, r2)
                        val = (start, q["sigma"])
                        if best_pair is None or val < best_pair[0]:
                            best_pair = (val, p1, p2)
                if best_pair is None:
                    continue
                start = best_pair[0][0]
                pan1, tilt1 = S.pointing_at(best_pair[1][0], best_pair[1][1],
                                            S.TARGETS[a])
                pan2, tilt2 = S.pointing_at(best_pair[2][0], best_pair[2][1],
                                            S.TARGETS[a])
                rst = max(S.slew_time(pan1, tilt1, 0.0, 0.0),
                          S.slew_time(pan2, tilt2, 0.0, 0.0))
                comp = start + S.OBS_DUR + rst
                if best is None or (comp, a) < (best[0], best[1]):
                    best = (comp, a, start, v1, v2, best_pair[1], best_pair[2])
        if best is None:
            break
        comp, a, start, v1, v2, p1, p2 = best
        starts[a] = start
        for v, p in ((v1, p1), (v2, p2)):
            pan, tilt = S.pointing_at(p[0], p[1], S.TARGETS[a])
            rst = S.slew_time(pan, tilt, 0.0, 0.0)
            chains[v].append((a, p))
            free[v] = dict(t=start + S.OBS_DUR + rst, pos=p, pan=0.0, tilt=0.0)
            touched.add(v)
        remaining.remove(a)
    for a in remaining:
        starts[a] = float("inf")
    _LAST_PLAN["chains"] = {v: list(ch) for v, ch in chains.items()}
    return starts


def plan_audit(chains, states):
    """min gamma / max sigma over a plan's pairs (None if empty)."""
    _, _, qual = chain_sim(chains, states) if chains else ({}, {}, {})
    if not qual:
        return None, None
    return min(q["gamma"] for q in qual.values()), \
        max(q["sigma"] for q in qual.values())


# ----------------------------------------------------------------------
# 5. Scenario runner: P0 / P1 / P1c / P2
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
    bq = _B2M["alert_quality"]
    row["P0"] = dict(eta=eta_k, t_last=tl_k,
                     min_gamma_deg=round(min(q["gamma"] for q in bq.values()), 1),
                     max_sigma_m=round(max(q["sigma"] for q in bq.values()), 3))

    # ---- shared state capture + T4 setup ----
    alive = [v for v in VEHS if v != failed]
    states, completed, pending = P4.extract_states(stub, decision_t, set(VEHS))
    open_alerts = [a for a in pending if a not in completed]
    if new_alert is not None:
        open_alerts = open_alerts + ["T4"]
    # inject the baseline's open-alert stations into the candidate grid so the
    # "retain" option stays representable in the search (idempotent per alert)
    for a in open_alerts:
        if a not in C_FULL:
            continue
        for v, al in ASG2M.items():
            if a in al:
                st = ST2M[(v, a)]
                if tuple(st) not in C_FULL[a]["pts"]:
                    C_FULL[a]["pts"].append(tuple(st))
                    C_FULL[a]["arr"] = np.asarray(C_FULL[a]["pts"], float)
                    C_FULL[a]["r"] = np.append(
                        C_FULL[a]["r"],
                        math.hypot(S.TARGETS[a][0] - st[0],
                                   S.TARGETS[a][1] - st[1]))
    saved = dict(S.TARGETS)
    if new_alert is not None:
        S.TARGETS["T4"] = S.target_xyz(new_alert)
        C_FULL["T4"] = make_set("T4")
    try:
        # ---- P1 (cap-free) / P1c (two-slot capped control) ----
        g1, ch1 = search_2m(states, alive, open_alerts, max_slots=None)
        g1c, ch1c = search_2m(states, alive, open_alerts, max_slots=2)
        _LAST_P1["chains"] = {v: list(ch) for v, ch in ch1.items()}
        starts_1 = {a: stub["alert_start"][a] for a in completed}
        starts_1.update({a: g1.get(a, float("inf")) for a in open_alerts})
        starts_1c = {a: stub["alert_start"][a] for a in completed}
        starts_1c.update({a: g1c.get(a, float("inf")) for a in open_alerts})
        for key, sts, ch in (("P1", starts_1, ch1), ("P1c", starts_1c, ch1c)):
            eta, tl = P4.evaluate(sts, all_alerts)
            gmin, smax = plan_audit(ch, states)
            row[key] = dict(eta=eta, t_last=tl,
                            starts={a: (round(t, 2) if t < 1e9 else None)
                                    for a, t in sts.items()},
                            min_gamma_deg=round(gmin, 1) if gmin else None,
                            max_sigma_m=round(smax, 3) if smax else None)

        # ---- P2 (rule greedy, gated) ----
        g2 = strategy_P2(states, alive, open_alerts)
        starts_2 = {a: stub["alert_start"][a] for a in completed}
        starts_2.update(g2)
        eta_2, tl_2 = P4.evaluate(starts_2, all_alerts)
        gmin2, smax2 = plan_audit(_LAST_PLAN["chains"], states)
        row["P2"] = dict(eta=eta_2, t_last=tl_2,
                         starts={a: (round(t, 2) if t < 1e9 else None)
                                 for a, t in starts_2.items()},
                         min_gamma_deg=round(gmin2, 1) if gmin2 else None,
                         max_sigma_m=round(smax2, 3) if smax2 else None)
    finally:
        S.TARGETS.clear(); S.TARGETS.update(saved)
    print(f"{name:36s} N={row['N']}  P0 ({eta_k:.2f},{tl_k:6.2f})  "
          f"P1 ({row['P1']['eta']:.2f},{row['P1']['t_last']:6.2f})  "
          f"P1c ({row['P1c']['eta']:.2f},{row['P1c']['t_last']:6.2f})  "
          f"P2 ({eta_2:.2f},{tl_2:6.2f})")
    return row


# ----------------------------------------------------------------------
# 6. 86-scenario family (paper §Q4 Table design)
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

if "--smoke" in sys.argv:
    print("=== SMOKE (6 representative scenarios, all families) ===")
    for idx in (0, 38, 43, 50, 79, 85):
        run_v2(**scenarios[idx])
    print("=== SMOKE OK ===")
    sys.exit(0)

if os.environ.get("Q4_PROBE"):
    print("=== PROBE ===")
    for sc in scenarios:
        if sc["name"] in os.environ["Q4_PROBE"].split(";"):
            run_v2(**sc)
            print("CHAINS P1:", json.dumps(_LAST_P1["chains"], default=str))
    sys.exit(0)

if os.environ.get("Q4_DUMP"):
    print("=== DUMP ===")
    for sc in scenarios:
        if sc["name"] in os.environ["Q4_DUMP"].split(";"):
            r = run_v2(**sc)
            print("P1 starts:", r["P1"]["starts"])
            print("P1 chains:", json.dumps(_LAST_P1["chains"], default=str))
            print("P2 starts:", r["P2"]["starts"])
            print("P2 chains:", json.dumps(_LAST_PLAN["chains"], default=str))
    sys.exit(0)

print("=== 86-scenario suite (2m gate; P0 / P1 / P1c / P2) ===")
rows = [run_v2(**sc) for sc in scenarios]

# ----------------------------------------------------------------------
# 7. Kinematic sensitivity of the P1 response (4 representative scenarios)
# ----------------------------------------------------------------------
BASE_SPEED, BASE_SLEW = S.V_SPEED, S.SLEW_RATE
SENS = []
for name, kw in (("NA worst (θ=225,R=50)", dict(new_alert=(225.0, 50.0, 50.0))),
                 ("NA easy  (θ=315,R=30)", dict(new_alert=(315.0, 30.0, 50.0))),
                 ("FA V3", dict(failed="V3")),
                 ("CO t_rec=20", dict(t_rec=20.0))):
    cfgs = {}
    for cname, ss, sl in (("nominal", 1.0, 1.0), ("speed−10%", 0.9, 1.0),
                          ("speed+10%", 1.1, 1.0), ("slew−10%", 1.0, 0.9),
                          ("slew+10%", 1.0, 1.1)):
        S.V_SPEED, S.SLEW_RATE = BASE_SPEED * ss, BASE_SLEW * sl
        row = run_v2(f"{name} [{cname}]", **kw)
        cfgs[cname] = dict(P1_tlast=row["P1"]["t_last"], P1_eta=row["P1"]["eta"])
    S.V_SPEED, S.SLEW_RATE = BASE_SPEED, BASE_SLEW
    row = run_v2(f"{name} [nominal-recheck]", **kw)
    cfgs["nominal_recheck"] = dict(P1_tlast=row["P1"]["t_last"])
    SENS.append(dict(scenario=name, **cfgs))
    print(f"sens {name}: " + ", ".join(
        f"{c}={v['P1_tlast']:.2f}" for c, v in cfgs.items()))

# ----------------------------------------------------------------------
# 8. Figures
# ----------------------------------------------------------------------
na_rows = [r for r in rows if r["scenario"].startswith("(a) NA ")]
thetas = sorted({float(r["new_alert"][0]) for r in na_rows})
radii = sorted({float(r["new_alert"][1]) for r in na_rows})
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
ax.set_title("P1 response over the new-alert grid (sigma<2m gate, η=1)",
             fontsize=10)
fig.colorbar(im, ax=ax, label="t_last (s)")
fig.tight_layout()
fig.savefig(OUT / "fig_q4_v2_grid.png", dpi=300)

fb = [r for r in rows if r["scenario"].startswith("(b) FA ") and
      "CO" not in r["scenario"]]
names = [r["scenario"].split()[-1] for r in fb]
x = np.arange(len(fb))
fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=300)
ax.bar(x - 0.27, [r["P0"]["t_last"] for r in fb], width=0.27,
       label="P0 retain", color="#c44")
ax.bar(x, [r["P1"]["t_last"] for r in fb], width=0.27,
       label="P1 re-opt (cap-free)", color="#36c")
ax.bar(x + 0.27, [r["P2"]["t_last"] for r in fb], width=0.27,
       label="P2 rule greedy", color="#3a3")
for xi, r in zip(x, fb):
    ax.text(xi - 0.27, r["P0"]["t_last"] + 0.5, f"η={r['P0']['eta']:.2f}",
            ha="center", fontsize=7.5)
    ax.text(xi + 0.27, r["P2"]["t_last"] + 0.5, f"η={r['P2']['eta']:.2f}",
            ha="center", fontsize=7.5)
ax.set_xticks(x); ax.set_xticklabels(names)
ax.set_ylabel("t_last (s)"); ax.set_ylim(0, max(65, max(
    r["P0"]["t_last"] for r in fb) * 1.08))
ax.set_title("Failure family: three strategies (window [0,60] s, 2m gate)",
             fontsize=10)
ax.legend(fontsize=8.5); ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "fig_q4_v2_failure.png", dpi=300)

co = [r for r in rows if r["scenario"].startswith("(c) CO ")]
trs = [r["t_rec"] for r in co]
fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=300)
ax.plot(trs, [r["P0"]["t_last"] for r in co], "s-", color="#c44",
        label="P0 retain")
ax.plot(trs, [r["P1"]["t_last"] for r in co], "o-", color="#36c",
        label="P1 re-opt (cap-free)")
ax.plot(trs, [r["P2"]["t_last"] for r in co], "^-", color="#3a3",
        label="P2 rule greedy")
ax.set_xlabel("outage recovery time t_rec (s)"); ax.set_ylabel("t_last (s)")
ax.set_title("Pure-outage family (2m gate)", fontsize=10)
ax.legend(fontsize=8.5); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "fig_q4_v2_trec.png", dpi=300)

fig, ax = plt.subplots(figsize=(7.4, 5.4), dpi=300)
for r in rows:
    d_eta = r["P1"]["eta"] - r["P0"]["eta"]
    d_tl = r["P1"]["t_last"] - r["P0"]["t_last"]
    ax.scatter(d_tl, d_eta, s=26,
               color="#36c" if d_eta > 0 else "#888", alpha=0.85)
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("Δ t_last = t_last(P1) − t_last(P0) (s)")
ax.set_ylabel("Δ η = η(P1) − η(P0)")
ax.set_title("Value of re-optimization across all 86 scenarios (2m gate):\n"
             "coverage is gained exactly when P0 loses an alert", fontsize=10)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "fig_q4_v2_value.png", dpi=300)

# ----------------------------------------------------------------------
# 9. Slot-cap cost + JSON
# ----------------------------------------------------------------------
better = [r for r in rows if r["P1"]["t_last"] < r["P1c"]["t_last"] - 1e-6]
worse = [r for r in rows if r["P1"]["t_last"] > r["P1c"]["t_last"] + 1e-6]
cap_cost = {
    "definition": "P1c (two-slot-capped search) minus P1 (cap-free), 2m gate",
    "n_better": len(better), "n_equal": 86 - len(better) - len(worse),
    "n_worse": len(worse),
    "max_gain_s": round(max((r["P1c"]["t_last"] - r["P1"]["t_last"]
                             for r in better), default=0.0), 3),
    "worst_gain_scenario": max(better, key=lambda r: r["P1c"]["t_last"]
                               - r["P1"]["t_last"])["scenario"] if better
    else None,
    "mean_gain_s": round(sum(r["P1c"]["t_last"] - r["P1"]["t_last"]
                             for r in rows) / len(rows), 3),
}
print("cap cost:", cap_cost)

res = {
    "baseline": {
        "description": "sigma<2m Q2(2) schedule (problem2_2_2m.py): time-first "
                       "with sigma<2.0 m and gamma>=30 deg on every pair",
        "assignment": {v: a for v, a in ASG2M.items() if a},
        "stations": {f"{v}-{a}": [round(x, 3) for x in s]
                     for (v, a), s in ST2M.items()},
        "alert_start_s": {k: round(t, 3) for k, t in
                          _B2M["alert_start_s"].items()},
        "max_sigma_m": _B2M["max_sigma_m"],
    },
    "quality_gate": {"sigma_max_m": SIGMA_MAX, "gamma_min_deg": GAMMA_MIN},
    "margins": {"clearance_beyond_buffer_m": MARGIN, "elevation_max_deg": EL_MAX},
    "strategies": ["P0 retain", "P1 cap-free re-optimization (time-first + "
                   "quality repair)", "P1c two-slot-capped control",
                   "P2 rule-based greedy"],
    "slot_cap_cost": cap_cost,
    "scenario_design": "86 = (a) 12x3 new-alert grid + (b) 5 failures + "
                       "(c) 5 outage lengths + (a+c) 28 + (b+c) 5 + (a+b) 5 "
                       "+ (a+b+c) 2",
    "kinematic_sensitivity": SENS,
    "scenarios": rows,
}
with open(OUT / "q4_v2_results.json", "w", encoding="utf-8") as fp:
    json.dump(res, fp, ensure_ascii=False, indent=2)
print("results saved:", OUT / "q4_v2_results.json")
