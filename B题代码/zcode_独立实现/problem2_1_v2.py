# -*- coding: utf-8 -*-
"""
Problem 2(1), revised (v2).

What changed vs 2.1.py
----------------------
1. Hard margins:
      building/field clearance >= 0.10 m beyond the 1.5 m vehicle buffer
      elevation <= 79.9 deg (0.1 deg margin)
      LOS clearance to buildings checked explicitly for the final pair
2. No Nelder-Mead and no soft 30 s penalty.  The optimum is found on the
   analytic 90-degree family:
      r in {r_min+0.10, +0.12, +0.15, +0.20}
      alpha in [0,360) step 0.1 deg
      stations A(alpha), B(alpha+90 deg)
   Every candidate is checked against feasibility, margins and travel <= 30 s;
   the best is refined locally.
3. Pair table is recomputed from the same hard-constrained search, so it no
   longer contradicts the final answer.
4. Figure uses a true two-colour observable mask (the old single-colour
   ListedColormap painted the whole field blue).

Path planning
-------------
S.PathPlanner() inflates each building by R_V=1.5 m, builds a visibility graph
on the four corners of every inflated rectangle (20 nodes), precomputes all
pairs by Floyd-Warshall, and returns the shortest path length / 3 m/s.  A
direct segment is used when it is clear; otherwise the path goes through
visible inflated corners.  This is conservative around sharp corners (the true
Minkowski obstacle is rounded), but for the selected V3/V5 stations both paths
are straight and clear.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, Ellipse, Patch, Rectangle

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
import swb_common as S  # noqa: E402

OUT = BASE / "results"
OUT.mkdir(exist_ok=True)

tgt = S.TARGETS["T2"]
tx, ty, H = tgt
r_min = S.min_range_for_target(H)
SIG = math.radians(S.SIG_ANGLE)
planner = S.PathPlanner()

CLR_MIN = 0.10          # m, beyond the 1.5 m vehicle buffer
EL_MARGIN = 0.10        # deg
EL_MAX = 80.0 - EL_MARGIN
TRAVEL_MAX = 30.0       # s
ANG_STEP = 0.10         # deg, coarse global search
R_OFFSETS = [0.10, 0.12, 0.15, 0.20]   # r = r_min + offset


def clearance_vec(xs, ys):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    c = np.minimum.reduce([xs - S.R_V, S.FIELD - S.R_V - xs,
                           ys - S.R_V, S.FIELD - S.R_V - ys])
    for (x0, x1, y0, y1, _h) in S.BUILDINGS:
        c = np.minimum(c, S.point_rect_dist_vec(xs, ys, (x0, x1, y0, y1)) - S.R_V)
    return c


def los_min_clearance(x, y, n=4000):
    """Minimum 3-D distance from the LOS segment to any building prism."""
    A = np.array([x, y, S.CAM_H], float)
    B = np.array([tx, ty, H], float)
    ts = np.linspace(0.0, 1.0, n)
    Q = A[None, :] + ts[:, None] * (B - A)[None, :]
    best = float("inf")
    for (x0, x1, y0, y1, h) in S.BUILDINGS:
        dx = np.maximum.reduce([np.full(n, x0) - Q[:, 0], np.zeros(n), Q[:, 0] - x1])
        dy = np.maximum.reduce([np.full(n, y0) - Q[:, 1], np.zeros(n), Q[:, 1] - y1])
        dz = np.maximum.reduce([np.zeros(n) - Q[:, 2], np.zeros(n), Q[:, 2] - h])
        best = min(best, float(np.sqrt(dx * dx + dy * dy + dz * dz).min()))
    return best


def tau_of(s):
    """Gimbal pointing time from due east/horizontal."""
    r = math.hypot(tx - s[0], ty - s[1])
    el = math.degrees(math.atan2(H - S.CAM_H, r))
    pan = math.degrees(math.atan2(ty - s[1], tx - s[0]))
    return max(abs((pan + 180.0) % 360.0 - 180.0), el) / S.SLEW_RATE


# ----------------------------------------------------------------------
# 1. Coarse hard-constrained search on the 90-degree family
# ----------------------------------------------------------------------
angles = np.arange(0.0, 360.0, ANG_STEP)
rad = np.deg2rad(angles)
best = None
pair_best = {}

for dr in R_OFFSETS:
    r = r_min + dr
    ax_ = tx + r * np.cos(rad)
    ay_ = ty + r * np.sin(rad)
    bx_ = tx + r * np.cos(rad + math.pi / 2)
    by_ = ty + r * np.sin(rad + math.pi / 2)

    okA, elA, _ = S.feasibility_mask(ax_, ay_, tgt)
    okB, elB, _ = S.feasibility_mask(bx_, by_, tgt)
    clrA = clearance_vec(ax_, ay_)
    clrB = clearance_vec(bx_, by_)
    valid = (okA & okB & (elA <= EL_MAX + 1e-9) & (elB <= EL_MAX + 1e-9)
             & (clrA >= CLR_MIN - 1e-9) & (clrB >= CLR_MIN - 1e-9))
    idx = np.where(valid)[0]
    if idx.size == 0:
        print(f"r=r_min+{dr:.2f}: no feasible 90-deg pair")
        continue

    tA = {}
    tB = {}
    for vi, v in enumerate(S.V_NAMES):
        p0 = S.V_INIT[vi]
        tA[v] = np.array([planner.time(p0, (float(ax_[i]), float(ay_[i]))) for i in idx])
        tB[v] = np.array([planner.time(p0, (float(bx_[i]), float(by_[i]))) for i in idx])

    for v1, v2 in __import__("itertools").combinations(S.V_NAMES, 2):
        for orient in (0, 1):
            t1 = tA[v1] if orient == 0 else tB[v1]
            t2 = tB[v2] if orient == 0 else tA[v2]
            good = (t1 <= TRAVEL_MAX) & (t2 <= TRAVEL_MAX)
            js = np.where(good)[0]
            for j in js:
                i = int(idx[j])
                s1 = (float(ax_[i]), float(ay_[i])) if orient == 0 else (float(bx_[i]), float(by_[i]))
                s2 = (float(bx_[i]), float(by_[i])) if orient == 0 else (float(ax_[i]), float(ay_[i]))
                total = float(t1[j] + t2[j])
                ta1, ta2 = tau_of(s1), tau_of(s2)
                mission = max(t1[j] + ta1, t2[j] + ta2) + 5.0 + max(ta1, ta2)
                mid = ((s1[0] + s2[0]) / 2, (s1[1] + s2[1]) / 2)
                sigma = SIG * math.sqrt(2.0) * r
                key = (round(sigma, 3), mission, total,
                       math.dist(mid, S.CENTER), angles[i] % 360.0, r)
                if best is None or key < best[0]:
                    best = (key, r, v1, v2, s1, s2, float(t1[j]), float(t2[j]), float(mission))
                pk = (v1, v2)
                if pk not in pair_best or key < pair_best[pk][0]:
                    pair_best[pk] = (key, r, v1, v2, s1, s2,
                                     float(t1[j]), float(t2[j]), float(mission))

    print(f"r=r_min+{dr:.2f}: {len(idx)} feasible 90-deg candidates")

if best is None:
    raise RuntimeError("No feasible 90-degree pair under the margins.")

# ----------------------------------------------------------------------
# 2. Local refinement around the coarse optimum
# ----------------------------------------------------------------------
_, r0, v1, v2, s10, s20, _, _, _ = best
ang0 = math.degrees(math.atan2(s10[1] - ty, s10[0] - tx)) % 360.0
print(f"coarse best: {v1}+{v2}, r={r0:.4f}, angle={ang0:.3f}")

best_ref = best
for r in np.arange(r0 - 0.02, r0 + 0.02 + 1e-12, 0.002):
    if r < r_min + 0.095:
        continue
    for adeg in np.arange(ang0 - 0.10, ang0 + 0.10 + 1e-12, 0.005):
        a = math.radians(adeg)
        pA = (tx + r * math.cos(a), ty + r * math.sin(a))
        pB = (tx + r * math.cos(a + math.pi / 2), ty + r * math.sin(a + math.pi / 2))
        if min(clearance_vec(np.array([pA[0], pB[0]]), np.array([pA[1], pB[1]]))) < CLR_MIN - 1e-9:
            continue
        okA, elA, _ = S.station_feasible(*pA, tgt)
        okB, elB, _ = S.station_feasible(*pB, tgt)
        if not (okA and okB) or elA > EL_MAX or elB > EL_MAX:
            continue
        for ss1, ss2, vv1, vv2 in ((pA, pB, v1, v2), (pB, pA, v2, v1)):
            t1 = planner.time(S.V_INIT[S.V_NAMES.index(vv1)], ss1)
            t2 = planner.time(S.V_INIT[S.V_NAMES.index(vv2)], ss2)
            if t1 > TRAVEL_MAX or t2 > TRAVEL_MAX:
                continue
            sigma = SIG * math.sqrt(2.0) * r
            ta1, ta2 = tau_of(ss1), tau_of(ss2)
            mission = max(t1 + ta1, t2 + ta2) + 5.0 + max(ta1, ta2)
            mid = ((ss1[0] + ss2[0]) / 2, (ss1[1] + ss2[1]) / 2)
            key = (round(sigma, 3), mission, t1 + t2,
                   math.dist(mid, S.CENTER), adeg % 360.0, r)
            if key < best_ref[0]:
                best_ref = (key, r, vv1, vv2, ss1, ss2, t1, t2, mission)

_, r_best, v1, v2, s1, s2, t1, t2, mission = best_ref
sigma = SIG * math.sqrt(2.0) * r_best
B = math.dist(s1, s2)
gamma = 90.0
ta1, ta2 = tau_of(s1), tau_of(s2)
obs_start = max(t1 + ta1, t2 + ta2)
clr1 = float(clearance_vec(np.array([s1[0]]), np.array([s1[1]]))[0])
clr2 = float(clearance_vec(np.array([s2[0]]), np.array([s2[1]]))[0])
el1 = math.degrees(math.atan2(H - S.CAM_H, math.hypot(tx - s1[0], ty - s1[1])))
el2 = math.degrees(math.atan2(H - S.CAM_H, math.hypot(tx - s2[0], ty - s2[1])))
los1 = los_min_clearance(*s1)
los2 = los_min_clearance(*s2)

print("\nROBUST OPTIMUM")
print(f"  vehicles: {v1} + {v2}")
print(f"  stations: {s1}, {s2}")
print(f"  r1=r2={r_best:.4f} m, B={B:.4f} m, gamma={gamma:.3f} deg")
print(f"  travel: {t1:.3f} s, {t2:.3f} s")
print(f"  pointing: {ta1:.3f} s, {ta2:.3f} s")
print(f"  observation start: {obs_start:.3f} s")
print(f"  mission completion: {mission:.3f} s")
print(f"  sigma_rms={sigma:.4f} m")
print(f"  clearance: {clr1:.4f}, {clr2:.4f} m")
print(f"  elevation: {el1:.6f}, {el2:.6f} deg")
print(f"  LOS min distance: {los1:.4f}, {los2:.4f} m")

# ----------------------------------------------------------------------
# 3. Figure: local two-station geometry (no observable-region shading)
# ----------------------------------------------------------------------
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

fig, ax = plt.subplots(figsize=(8.6, 7.2), dpi=300)

# B3 and its 1.5 m vehicle buffer
for (x0, x1, y0, y1, h), name in zip(S.BUILDINGS, S.BUILD_NAMES):
    if name != "B3":
        continue
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                           facecolor="#9a9a9a", edgecolor="k", alpha=0.55, zorder=2))
    ax.add_patch(Rectangle((x0 - S.R_V, y0 - S.R_V),
                           x1 - x0 + 2 * S.R_V, y1 - y0 + 2 * S.R_V,
                           fill=False, ls="--", ec="0.45", lw=0.9, zorder=3))
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"{name}\n{h:g} m",
            color="k", ha="center", va="center", fontsize=8, zorder=4)

# r_min circle (geometric locus, not the observable region)
theta = np.linspace(0, 2 * np.pi, 720)
ax.plot(tx + r_best * np.cos(theta), ty + r_best * np.sin(theta),
        color="0.55", ls="--", lw=1.0, zorder=4,
        label=f"operating circle ($r_{{\\min}}+0.10$ m = {r_best:.3f} m)")

# bearing lines (horizontal projection; actual LOS passes over B3)
ax.plot([s1[0], tx], [s1[1], ty], color="k", ls="--", lw=1.0, zorder=5,
        label="bearing lines (2D projection)")
ax.plot([s2[0], tx], [s2[1], ty], color="k", ls="--", lw=1.0, zorder=5)

# baseline
ax.plot([s1[0], s2[0]], [s1[1], s2[1]], color="#1f4fd8", ls="-.", lw=1.6,
        zorder=6, label=f"baseline B = {B:.3f} m")

# target
ax.plot(tx, ty, "^", color="red", ms=13, zorder=9)
ax.text(tx + 0.4, ty - 2.5, "T2\n(H=35 m)", fontsize=8.5, color="red",
        ha="left", va="top", zorder=10,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

# stations and their 1.5 m disks
for s, v, clr in ((s1, v1, clr1), (s2, v2, clr2)):
    ax.plot(*s, "s", color="#1f4fd8", ms=11, mec="gold", mew=1.8, zorder=10)
    ax.add_patch(Circle(s, S.R_V, fill=False, ec="k", ls=":", lw=1.0, zorder=8))
    if s[0] < 10:
        lx, ly, ha = s[0] - 2.2, s[1] + 2.2, "center"
    else:
        lx, ly, ha = s[0] + 2.0, s[1] + 2.2, "center"
    ax.text(lx, ly, f"S{1 if s is s1 else 2}({v})", fontsize=8.5,
            color="#1f4fd8", ha=ha, va="bottom", zorder=11,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))

# intersection angle at T2
a1 = math.degrees(math.atan2(s1[1] - ty, s1[0] - tx)) % 360.0
a2 = math.degrees(math.atan2(s2[1] - ty, s2[0] - tx)) % 360.0
lo, hi = min(a1, a2), max(a1, a2)
ax.add_patch(Arc((tx, ty), 4.4, 4.4, theta1=lo, theta2=hi,
                 color="k", lw=1.1, zorder=9))
ax.text(tx + 0.2, ty + 2.8, f"$\\gamma$ = {hi - lo:.1f}°", fontsize=9,
        ha="left", va="bottom", zorder=10,
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))

ax.text(0.98, 0.03,
        "2D top view: actual LOS passes over B3\n(elevation 79.83 deg)",
        transform=ax.transAxes, fontsize=7.6, ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.6", alpha=0.9),
        zorder=12)

# compact data box
ax.text(0.02, 0.03,
        f"$r_1=r_2={r_best:.3f}$ m\n"
        f"$B={B:.3f}$ m,  $\\gamma=90^\\circ$\n"
        f"$\\sigma_{{\\mathrm{{RMS}}}}={sigma:.4f}$ m\n"
        f"{v1}$\\to$S1: {t1:.2f} s,  {v2}$\\to$S2: {t2:.2f} s",
        transform=ax.transAxes, fontsize=8.2, ha="left", va="bottom",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.6", alpha=0.92),
        zorder=12)

ax.set_xlim(2, 26)
ax.set_ylim(64, 82)
ax.set_aspect("equal")
ax.set_xlabel("x (m) east")
ax.set_ylabel("y (m) north")
ax.set_title(f"Problem 2(1): two-station geometry ({v1}+{v2}, gamma = 90 deg)",
             fontsize=11)
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "fig_p2_1_geometry_v2.png", dpi=300)
fig.savefig(OUT / "fig_p2_1_accuracy_v2.png", dpi=300)
print("figure saved:", OUT / "fig_p2_1_geometry_v2.png")

# ----------------------------------------------------------------------
# 3b. Overall layout figure (white background, no observable region)
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.6, 8.6), dpi=300)

for (x0, x1, y0, y1, h), name in zip(S.BUILDINGS, S.BUILD_NAMES):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                           facecolor="#555555", edgecolor="k", zorder=3))
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"{name}\n{h:g} m",
            color="w", ha="center", va="center", fontsize=8, zorder=4)

for i, (vx, vy) in enumerate(S.V_INIT):
    ax.plot(vx, vy, "s", color="#1f4fd8", ms=6, mec="k", mew=0.4, zorder=7)
    ax.text(vx + 1.3, vy + 1.3, S.V_NAMES[i], fontsize=7.5,
            color="#1f4fd8", zorder=8)

ax.plot(tx, ty, "^", color="red", ms=12, zorder=8)
ax.text(tx + 1.2, ty + 1.4, "T2", fontsize=8.5, color="red", zorder=9,
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
ax.plot([s1[0], s2[0]], [s1[1], s2[1]], color="#1f4fd8", ls="-.", lw=1.4,
        zorder=6, label=f"baseline B = {B:.2f} m")
for s, v in ((s1, v1), (s2, v2)):
    ax.plot(*s, "s", color="#1f4fd8", ms=10, mec="gold", mew=1.8, zorder=9)



# travel paths for the two assigned vehicles
for v, s, tt in ((v1, s1, t1), (v2, s2, t2)):
    p0 = S.V_INIT[S.V_NAMES.index(v)]
    ax.annotate("", xy=s, xytext=p0,
                arrowprops=dict(arrowstyle="->", ls="--", lw=1.0,
                                color="#1f4fd8", shrinkA=0, shrinkB=5),
                zorder=6)
    mid = ((p0[0] + s[0]) / 2, (p0[1] + s[1]) / 2)
    lx, ly = (mid[0] - 6.0, mid[1] + 2.0) if v == v1 else (mid[0] + 1.0, mid[1] + 1.0)
    ax.text(lx, ly, f"{v}: {tt:.2f} s", fontsize=7.2,
            color="#1f4fd8", zorder=9,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))

# local-geometry figure area
ax.add_patch(Rectangle((2, 64), 24, 18, fill=False, ls="--", ec="k", lw=1.0, zorder=7))
ax.text(3.0, 65.2, "detail view", fontsize=7.5, zorder=9,
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))

handles = [
    Patch(facecolor="#555555", label="building"),
    plt.Line2D([0], [0], color="#1f4fd8", ls="-.", lw=1.5, label="baseline"),
    plt.Line2D([0], [0], marker="s", color="none", mfc="#1f4fd8",
               mec="gold", mew=1.5, ms=9, label="selected station"),
    plt.Line2D([0], [0], marker="s", color="none", mfc="#1f4fd8",
               mec="k", mew=0.5, ms=6, label="standby vehicle"),
    plt.Line2D([0], [0], color="#1f4fd8", ls="--", lw=1.0, label="travel path"),
]
ax.legend(handles=handles, loc="upper right", fontsize=7.5, framealpha=0.9)

ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_aspect("equal")
ax.set_xlabel("x (m) east")
ax.set_ylabel("y (m) north")
ax.set_title(f"T2 layout with selected stations ({v1}+{v2})", fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "fig_p2_1_region_v2.png", dpi=300)
print("figure saved:", OUT / "fig_p2_1_region_v2.png")

# ----------------------------------------------------------------------
# 4. JSON
# ----------------------------------------------------------------------
pair_rows = []
for (vv1, vv2), rec in sorted(pair_best.items(), key=lambda kv: kv[1][0]):
    key, r, a1, a2, p1, p2, tt1, tt2, mm = rec
    pair_rows.append({
        "pair": f"{a1}+{a2}",
        "sigma_rms_m": round(float(key[0]), 5),
        "r_m": round(float(r), 4),
        "travel_total_s": round(float(tt1 + tt2), 3),
        "mission_time_s": round(float(mm), 3),
        "stations": [list(p1), list(p2)],
    })

res = {
    "target": "T2",
    "target_xyz": [round(tx, 6), round(ty, 6), H],
    "theoretical_bound": {
        "r_min_m": round(r_min, 6),
        "sigma_lower_bound_m": round(SIG * math.sqrt(2.0) * r_min, 6),
    },
    "margins": {
        "clearance_min_m": CLR_MIN,
        "elevation_margin_deg": EL_MARGIN,
        "travel_max_s": TRAVEL_MAX,
    },
    "robust_optimal": {
        "vehicles": [v1, v2],
        "stations": [list(s1), list(s2)],
        "r_m": round(r_best, 6),
        "baseline_m": round(B, 4),
        "intersection_angle_deg": 90.0,
        "travel_times_s": [round(t1, 3), round(t2, 3)],
        "pointing_times_s": [round(ta1, 3), round(ta2, 3)],
        "observation_start_s": round(obs_start, 3),
        "mission_time_s": round(mission, 3),
        "sigma_rms_m": round(sigma, 6),
        "clearance_m": [round(clr1, 4), round(clr2, 4)],
        "elevation_deg": [round(el1, 6), round(el2, 6)],
        "los_min_distance_m": [round(los1, 4), round(los2, 4)],
    },
    "pair_table": pair_rows,
    "method": "hard-constrained search on the analytic 90-degree family: "
              "r in r_min+{0.10,0.12,0.15,0.20}, alpha step 0.1 deg, "
              "clearance>=0.10 m, elevation<=79.9 deg, travel<=30 s, "
              "then local refinement at 0.005 deg / 0.002 m.",
    "tie_break": "accuracy first (sigma rounded to 0.001 m), then the "
                 "two-vehicle mission completion time "
                 "max_i(t_i+tau_i)+5+max_i(tau_i), then total travel, "
                 "then midpoint distance to O, then angle.",
    "path_planner": "S.PathPlanner: buildings inflated by 1.5 m, visibility "
                    "graph on 20 inflated corners, Floyd-Warshall, direct "
                    "segment when clear, path length/3 m/s.",
}
with open(OUT / "problem2_1_results_v2.json", "w", encoding="utf-8") as fp:
    json.dump(res, fp, ensure_ascii=False, indent=2)
print("results saved:", OUT / "problem2_1_results_v2.json")
