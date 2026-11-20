# -*- coding: utf-8 -*-
"""
Problem 2(1) — Two-station intersection accuracy model + optimal pair of
stations (and vehicles) for alert T2, subject to travel time <= 30 s.

Accuracy model
--------------
Measured bearings  th_i = th_i^true + eps_i,  eps_i ~ N(0, sigma^2) indep, sigma = 1 deg.
Intersection perturbation (linearized, exact to first order):
    dP = ( d1*eps1*u2 - d2*eps2*u1 ) / sin(gamma)
where d_i = horizontal range of station i, u_i = unit LOS vector (station->target),
gamma = angle between the two LOS lines at the target (0..180 deg).
Position-error covariance:
    C = sigma^2/sin^2(gamma) * ( d1^2 u2 u2^T + d2^2 u1 u1^T )
Accuracy metric (RMS horizontal position error):
    sigma_rms = sqrt(trace C) = sigma * sqrt(d1^2 + d2^2) / |sin gamma|   [m]
Rationale: trace C = E[||dP||^2] is the mean squared horizontal miss distance;
it is invariant to axis rotation, standard (GDOP-like), and reveals the two
failure modes (gamma -> 0/180 deg, long ranges).

Optimal pair: minimize sigma_rms over feasible stations s1, s2 with per-vehicle
travel time <= 30 s. Analytic floor: d_i >= (H-CAM_H)/tan(80 deg) and sin(gamma) <= 1
=> sigma_rms >= sigma * sqrt(2) * r_min = 0.1436 m, attained iff both stations sit on
the r_min circle with gamma = 90 deg (1-D feasible family -> tie-breaking rules).

Outputs: results/problem2_1_results.json, fig_p2_1_accuracy.png
"""
import json
import math
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import swb_common as S
from scipy.optimize import minimize

OUT = "results"
os.makedirs(OUT, exist_ok=True)

tgt = S.TARGETS["T2"]
tx, ty, H = tgt
SIG = math.radians(S.SIG_ANGLE)          # 1 deg in radians
planner = S.PathPlanner()

# ----------------------------------------------------------------------
# 0. Accuracy-model helpers
# ----------------------------------------------------------------------
def los_unit(sx, sy):
    d = math.hypot(tx - sx, ty - sy)
    return np.array([(tx - sx) / d, (ty - sy) / d]), d


def gamma_of(s1, s2):
    u1, d1 = los_unit(*s1)
    u2, d2 = los_unit(*s2)
    g = math.degrees(math.acos(np.clip(u1 @ u2, -1, 1)))
    return g, d1, d2, u1, u2


def cov_of(s1, s2):
    g, d1, d2, u1, u2 = gamma_of(s1, s2)
    C = SIG ** 2 / math.sin(math.radians(g)) ** 2 * \
        (d1 ** 2 * np.outer(u2, u2) + d2 ** 2 * np.outer(u1, u1))
    return C, g, d1, d2


def sigma_rms(s1, s2):
    C, g, d1, d2 = cov_of(s1, s2)
    return math.sqrt(max(np.trace(C), 0.0)), g, d1, d2


def metric_of(s1, s2):
    """Pure geometry metric sqrt(d1^2+d2^2)/sin(gamma) (meters per radian)."""
    g, d1, d2, _, _ = gamma_of(s1, s2)
    return math.sqrt(d1 ** 2 + d2 ** 2) / max(math.sin(math.radians(g)), 1e-9), g, d1, d2


# ----------------------------------------------------------------------
# 1. Candidate stations for T2 + travel-time feasibility per vehicle
# ----------------------------------------------------------------------
r_min = S.min_range_for_target(H)
print(f"T2 = ({tx:.3f}, {ty:.3f}, {H}); r_min (el=80deg) = {r_min:.4f} m")

cand = S.candidate_stations(tgt, ang_step=1.0, rad_step=0.5, r_max=70.0)
pts = cand["pts"]
print(f"feasible candidate stations: {len(pts)}")

travel = {}     # (veh, idx) -> path time (s)
veh_ok = {}
for vi, v in enumerate(S.V_NAMES):
    p0 = S.V_INIT[vi]
    times = []
    for i, p in enumerate(pts):
        euc = math.dist(p0, p)
        t = planner.time(p0, p) if euc <= 90.0 else 1e9   # path >= straight line
        travel[(v, i)] = t
    times = np.array([travel[(v, i)] for i in range(len(pts))])
    veh_ok[v] = bool((times <= 30.0).any())
    print(f"{v}: reachable candidates within 30 s: {(times <= 30.0).sum()}, "
          f"min travel {times.min():.1f} s -> {'OK' if veh_ok[v] else 'EXCLUDED'}")

usable = [v for v in S.V_NAMES if veh_ok[v]]
print("vehicles able to reach T2 within 30 s:", usable)

# ----------------------------------------------------------------------
# 2. Continuous optimization per vehicle pair (multi-start Nelder-Mead)
# ----------------------------------------------------------------------
def pair_objective(z, v1, v2):
    s1 = (z[0], z[1]); s2 = (z[2], z[3])
    pen = 0.0
    for s, v in ((s1, v1), (s2, v2)):
        ok, el, d = S.station_feasible(s[0], s[1], tgt)
        if not ok:
            pen += 100.0
        tt = planner.time(S.V_INIT[S.V_NAMES.index(v)], s)
        if tt > 30.0:
            pen += (tt - 30.0)
    m, g, d1, d2 = metric_of(s1, s2)
    return m + pen


pairs_result = []
from itertools import combinations
for v1, v2 in combinations(S.V_NAMES, 2):
    if v1 not in usable or v2 not in usable:
        continue
    i1v = [i for i in range(len(pts)) if travel[(v1, i)] <= 30.0]
    i2v = [i for i in range(len(pts)) if travel[(v2, i)] <= 30.0]
    if not i1v or not i2v:
        continue
    # multi-starts: for each vehicle keep the M candidates nearest to the target
    # (short ranges dominate the metric), then combine M x M
    M = 60
    i1v.sort(key=lambda i: (cand["r"][i], travel[(v1, i)]))
    i2v.sort(key=lambda i: (cand["r"][i], travel[(v2, i)]))
    i1v, i2v = i1v[:M], i2v[:M]
    best_starts = []
    for i in i1v:
        for j in i2v:
            if i == j:
                continue
            m, g, d1, d2 = metric_of(pts[i], pts[j])
            best_starts.append((m, i, j))
    best_starts.sort(key=lambda x: x[0])
    best_val, best_z = None, None
    for (m0, i, j) in best_starts[:12]:
        z0 = np.array([*pts[i], *pts[j]])
        res = minimize(pair_objective, z0, args=(v1, v2), method="Nelder-Mead",
                       options=dict(maxiter=400, xatol=1e-4, fatol=1e-6))
        if best_val is None or res.fun < best_val:
            best_val, best_z = float(res.fun), res.x
    s1o = (float(best_z[0]), float(best_z[1]))
    s2o = (float(best_z[2]), float(best_z[3]))
    pairs_result.append((best_val, v1, v2, s1o, s2o))
    print(f"pair {v1}+{v2}: metric {best_val:.4f} m/rad  "
          f"sigma_rms={SIG*best_val:.4f} m  stations ({s1o[0]:.2f},{s1o[1]:.2f}) ({s2o[0]:.2f},{s2o[1]:.2f})")

pairs_result.sort(key=lambda x: x[0])
best_val, v1, v2, s1o, s2o = pairs_result[0]
print(f"\noptimal vehicle pair: {v1} + {v2}, metric {best_val:.4f} m/rad")

# ----------------------------------------------------------------------
# 3. Analytic optimal family (both stations on r_min circle, gamma = 90 deg)
#    swept over ALL usable vehicle pairs + explicit tie-breaking rules
# ----------------------------------------------------------------------
# Tie-breaking (stated):
#  (T1) accuracy first — all family pairs attain the analytic floor;
#  (T2) minimize total travel time of the two assigned vehicles;
#  (T3) minimize the pair-midpoint distance to the dispatch center O;
#  (T4) lexicographic smallest (x1, y1, x2, y2).
ang_grid = np.arange(0.0, 360.0, 0.05)
a = np.deg2rad(ang_grid)
sA_all = np.stack([tx + r_min * np.cos(a), ty + r_min * np.sin(a)], axis=1)
sB_all = np.stack([tx + r_min * np.cos(a + math.pi / 2),
                   ty + r_min * np.sin(a + math.pi / 2)], axis=1)
okA, _, _ = S.feasibility_mask(sA_all[:, 0], sA_all[:, 1], tgt)
okB, _, _ = S.feasibility_mask(sB_all[:, 0], sB_all[:, 1], tgt)
both = okA & okB
print(f"analytic family: {int(both.sum())}/{len(ang_grid)} rotations feasible (A,B)")

family = []
for v1, v2 in combinations(usable, 2):
    p01, p02 = S.V_INIT[S.V_NAMES.index(v1)], S.V_INIT[S.V_NAMES.index(v2)]
    idxs = np.where(both)[0]
    for i in idxs:
        for (vv1, ss1, vv2, ss2) in ((v1, tuple(sA_all[i]), v2, tuple(sB_all[i])),
                                     (v1, tuple(sB_all[i]), v2, tuple(sA_all[i]))):
            pp1, pp2 = (p01, p02) if vv1 == v1 else (p02, p01)
            t1 = planner.time(pp1, ss1)
            t2 = planner.time(pp2, ss2)
            if t1 <= 30.0 + 1e-9 and t2 <= 30.0 + 1e-9:
                mid = ((ss1[0] + ss2[0]) / 2, (ss1[1] + ss2[1]) / 2)
                family.append((t1 + t2, math.dist(mid, S.CENTER),
                               ss1[0], ss1[1], ss2[0], ss2[1], vv1, vv2,
                               t1, t2, v1, v2))
print(f"family entries with feasible travel: {len(family)}")
if family:
    family.sort()
    f = family[0]
    print("global tie-break winner:")
    print(f"  {f[6]} -> station ({f[2]:.4f}, {f[3]:.4f})  travel {f[8]:.2f} s")
    print(f"  {f[7]} -> station ({f[4]:.4f}, {f[5]:.4f})  travel {f[9]:.2f} s")
    print(f"  total travel {f[0]:.2f} s, midpoint dist to O {f[1]:.2f} m")
    s1F = (f[2], f[3]); s2F = (f[4], f[5])
    v1, v2 = f[6], f[7]
else:
    s1F, s2F, f = s1o, s2o, None
    print("analytic family empty — using NM optimum")

# final answer = analytic family point (matches NM optimum value; use family w/ tie-break)
s1F_, s2F_ = (s1F, s2F) if f else (s1o, s2o)
srms, g, d1, d2 = sigma_rms(s1F_, s2F_)
C, _, _, _ = cov_of(s1F_, s2F_)
B = math.dist(s1F_, s2F_)
# Monte-Carlo validation of the linearized model
rng = np.random.default_rng(7)
u1, dd1 = los_unit(*s1F_)
u2, dd2 = los_unit(*s2F_)
g_, _, _, _, _ = gamma_of(s1F_, s2F_)
M = 200_000
e1 = rng.normal(0, math.radians(S.SIG_ANGLE), M)
e2 = rng.normal(0, math.radians(S.SIG_ANGLE), M)
dP = (dd1 * e1[:, None] * u2[None, :] - dd2 * e2[:, None] * u1[None, :]) / math.sin(math.radians(g_))
rad = np.hypot(dP[:, 0], dP[:, 1])
# end-to-end exact MC: perturb the measured bearings and intersect the two lines
th1 = math.degrees(math.atan2(ty - s1F_[1], tx - s1F_[0]))
th2 = math.degrees(math.atan2(ty - s2F_[1], tx - s2F_[0]))
m1 = rng.normal(0, S.SIG_ANGLE, M) + th1
m2 = rng.normal(0, S.SIG_ANGLE, M) + th2
# solve S1 + t1*u(m1) = S2 + t2*u(m2)
m1r, m2r = np.deg2rad(m1), np.deg2rad(m2)
a11, a12 = np.cos(m1r), -np.cos(m2r)
a21, a22 = np.sin(m1r), -np.sin(m2r)
bx_ = s2F_[0] - s1F_[0]
by_ = s2F_[1] - s1F_[1]
det = a11 * a22 - a12 * a21
t1 = (bx_ * a22 - a12 * by_) / det
Px = s1F_[0] + t1 * np.cos(m1r)
Py = s1F_[1] + t1 * np.sin(m1r)
rad_exact = np.hypot(Px - tx, Py - ty)
print(f"\nfinal: stations {s1F_}, {s2F_}; B = {B:.3f} m; gamma = {g_:.3f} deg; "
      f"d1={d1:.3f} d2={d2:.3f}")
print(f"sigma_rms analytic = {srms:.4f} m; MC-linear mean radial = {rad.mean():.4f} m; "
      f"MC-linear CEP50 = {np.median(rad):.4f} m")
print(f"MC-exact CEP50 = {np.median(rad_exact):.4f} m (validates the linearization)")

# ----------------------------------------------------------------------
# 4. Figure: observable region + chosen stations + error ellipses
# ----------------------------------------------------------------------
xs = np.arange(S.R_V, S.FIELD - S.R_V + 1e-9, 0.5)
X, Y = np.meshgrid(xs, xs)
ok, el, rr = S.feasibility_mask(X.ravel(), Y.ravel(), tgt)
OK = ok.reshape(X.shape)

fig, ax = plt.subplots(figsize=(9, 9), dpi=300)
ax.imshow(OK.astype(float), origin="lower", extent=[0, 100, 0, 100],
          cmap=matplotlib.colors.ListedColormap(["#cfe3ff"]), alpha=0.6, zorder=0)
for (x0, x1, y0, y1, h), name in zip(S.BUILDINGS, S.BUILD_NAMES):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="#555555",
                           edgecolor="k", zorder=3))
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"{name}\n{h} m", color="w",
            ha="center", va="center", fontsize=8, zorder=4)
ax.plot([S.CENTER[0], tx], [S.CENTER[1], ty], "r--", lw=1.0, zorder=5)
ax.plot(tx, ty, "^", color="red", ms=12, zorder=6)
ax.annotate(f"T2 ({tx:.2f}, {ty:.2f}, {H})", (tx, ty), xytext=(tx + 2, ty + 3), fontsize=9)
ax.plot(*S.CENTER, "ko", ms=5, zorder=6)
# error ellipses (1-sigma and 3-sigma) from covariance C
evals, evecs = np.linalg.eigh(C)
for s_ in (s1F_, s2F_):
    ax.plot(*s_, "s", color="blue", ms=9, zorder=6)
ang = math.degrees(math.atan2(evecs[1, 1], evecs[0, 1]))
for scale, ls in ((1.0, "-"), (3.0, "--")):
    ax.add_patch(Ellipse((tx, ty), width=2 * scale * math.sqrt(max(evals[1], 1e-12)),
                         height=2 * scale * math.sqrt(max(evals[0], 1e-12)), angle=ang,
                         fill=False, color="darkorange", lw=1.4, ls=ls, zorder=6,
                         label=f"{'1σ' if scale == 1 else '3σ'} error ellipse"))
ax.plot([s1F_[0], s2F_[0]], [s1F_[1], s2F_[1]], "b-.", lw=1.2, zorder=6,
        label=f"baseline B = {B:.2f} m")
for s_, v in ((s1F_, f[6] if f else v1), (s2F_, f[7] if f else v2)):
    vi = S.V_NAMES.index(v)
    ax.annotate(f"{v}: {planner.time(S.V_INIT[vi], s_):.1f} s", s_,
                xytext=(s_[0] + 3, s_[1] - 4), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.8))
ax.plot([], [], "s", color="blue", label=f"stations ({v1}+{v2})")
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_aspect("equal")
ax.set_xlabel("x (m) east"); ax.set_ylabel("y (m) north")
ax.set_title(f"Problem 2(1) — optimal two-station geometry for T2  "
             f"(γ={g_:.1f}°, σ_RMS={srms:.3f} m)")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_p2_1_accuracy.png", dpi=300)

# ----------------------------------------------------------------------
# 5. Save results
# ----------------------------------------------------------------------
res = {
    "target": "T2", "target_xyz": [tx, ty, H],
    "vehicles_usable_within_30s": usable,
    "pair_metric_table_m_per_rad": [
        {"pair": f"{p[1]}+{p[2]}", "metric": round(p[0], 4)} for p in pairs_result],
    "optimal": {
        "vehicles": [f[6] if f else v1, f[7] if f else v2],
        "stations": [list(s1F_), list(s2F_)],
        "baseline_m": round(B, 4),
        "intersection_angle_deg": round(g_, 3),
        "ranges_m": [round(d1, 4), round(d2, 4)],
        "elevations_deg": [round(S.station_feasible(*s1F_, tgt)[1], 3),
                           round(S.station_feasible(*s2F_, tgt)[1], 3)],
        "travel_times_s": [round(f[8], 3), round(f[9], 3)] if f else None,
        "sigma_rms_m": round(srms, 5),
        "metric_m_per_rad": round(best_val if not f else metric_of(s1F_, s2F_)[0], 5),
        "MC_CEP50_m": round(float(np.median(rad)), 5),
        "MC_mean_radial_m": round(float(rad.mean()), 5),
        "MC_exact_CEP50_m": round(float(np.median(rad_exact)), 5),
    },
    "tie_breaking": "among accuracy-optimal pairs: (1) min total travel time, "
                    "(2) min midpoint distance to dispatch center, (3) lexicographic",
    "degeneracy_notes": "gamma->0/180 deg: sin(gamma)->0, C diverges (near-parallel "
                        "LOS); boundary stations (el=80 deg): use margins and robust "
                        "formulation; see README.",
}
with open(f"{OUT}/problem2_1_results.json", "w", encoding="utf-8") as fp:
    json.dump(res, fp, ensure_ascii=False, indent=2)
print("figure + results saved.")
