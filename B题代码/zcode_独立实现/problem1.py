# -*- coding: utf-8 -*-
"""
Problem 1 — Observable region and optimal observation station for target T1.

Model
-----
A station s = (x, y) (vehicle disk center) is FEASIBLE for target T = (tx, ty, H) iff
  (C1) vehicle geometry : R_V <= x,y <= FIELD-R_V  and  dist(s, building footprint) >= R_V
  (C2) line of sight    : segment (s, CAM_H) -> T intersects no building prism
                          (tangency/edge contact counts as occlusion)
  (C3) gimbal elevation : EL_MIN <= el(s) <= EL_MAX,  el = atan2(H-CAM_H, ||T'-s||)
Observation quality metric  Q(s) = 1/r(s),  r = ||T' - s||  (horizontal camera-target range).
"Optimal" = argmax Q over the feasible set  <=>  argmin r  (monotone transformation).
Since el <= EL_MAX forces r >= (H-CAM_H)/tan(EL_MAX), the theoretical optimum is
r* = (H-CAM_H)/tan(80 deg) whenever some point of that circle is feasible.

Outputs: results/problem1_results.json, fig_p1_region.png
"""
import json
import math
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import swb_common as S

OUT = "results"
import os
os.makedirs(OUT, exist_ok=True)

tgt = S.TARGETS["T1"]
tx, ty, H = tgt
alert = S.ALERTS["T1"]

# ----------------------------------------------------------------------
# 1. Observable region on a 0.5 m grid (vectorized feasibility)
# ----------------------------------------------------------------------
xs = np.arange(S.R_V, S.FIELD - S.R_V + 1e-9, 0.5)
ys = np.arange(S.R_V, S.FIELD - S.R_V + 1e-9, 0.5)
X, Y = np.meshgrid(xs, ys)
ok, el, r = S.feasibility_mask(X.ravel(), Y.ravel(), tgt)
OK = ok.reshape(X.shape)
R = r.reshape(X.shape)
Q = np.where(OK, 1.0 / np.maximum(R, 1e-9), np.nan)

print(f"grid points: {OK.size}, feasible: {int(OK.sum())} "
      f"({100*OK.sum()/OK.size:.1f}% of field)")

# ----------------------------------------------------------------------
# 2. Optimal station: minimize r over the feasible set
#    grid minimum -> analytic refinement on the r_min circle
# ----------------------------------------------------------------------
feas_r = np.where(OK, R, np.inf)
i_min = np.unravel_index(np.argmin(feas_r), feas_r.shape)
r_grid = float(feas_r[i_min])
print(f"grid optimum: r = {r_grid:.4f} m at "
      f"({X[i_min]:.2f}, {Y[i_min]:.2f}), el = {el.reshape(X.shape)[i_min]:.3f} deg")

r_theory = S.min_range_for_target(H)
print(f"theoretical lower bound r* = (H-CAM_H)/tan(80deg) = {r_theory:.4f} m")

# feasible arc on the optimal circle: sweep angles finely, pick tie-break set
arc = []
for adeg in np.arange(0.0, 360.0, 0.02):
    a = math.radians(adeg)
    x = tx + r_theory * math.cos(a)
    y = ty + r_theory * math.sin(a)
    f, e, d = S.station_feasible(x, y, tgt)
    if f:
        arc.append((adeg, x, y, e))
print(f"optimal-circle feasible arc: {len(arc)} points "
      f"(angles {arc[0][0]:.2f}..{arc[-1][0]:.2f} deg)") if arc else print("arc empty")

# Tie-breaking rules (stated explicitly):
#   (T1) among all argmax-Q stations, minimize the distance to the dispatch
#        center O (fastest response / shortest re-deployment);
#   (T2) remaining ties: smallest azimuth (CCW from east), i.e. lexicographic
#        (angle, then x, then y).
best = None
for (adeg, x, y, e) in arc:
    dO = math.dist((x, y), S.CENTER)
    key = (round(dO, 6), adeg)
    if best is None or key < best[0]:
        best = (key, adeg, x, y, e, dO)
_, adeg, bx, by, bel, bdO = best
print(f"tie-break optimal station: ({bx:.4f}, {by:.4f}), azimuth {adeg:.2f} deg, "
      f"r = {r_theory:.4f} m, el = {bel:.4f} deg, dist to O = {bdO:.3f} m")
Q_star = 1.0 / r_theory
print(f"Q* = {Q_star:.6f} 1/m")

# ----------------------------------------------------------------------
# 3. Plot
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.6, 8.6), dpi=300)
# feasible region shading
cmap = matplotlib.colors.ListedColormap(["#cfe3ff"])
ax.imshow(OK.astype(float), origin="lower", extent=[0, 100, 0, 100],
          cmap=cmap, alpha=0.55, zorder=0)
# Q heatmap contours on feasible part
Qm = np.ma.masked_invalid(Q)
cf = ax.contourf(X, Y, Qm, levels=14, cmap="viridis", alpha=0.75, zorder=1)
plt.colorbar(cf, ax=ax, label="Observation quality  Q = 1/r  (1/m)", shrink=0.8)
# buildings
for (x0, x1, y0, y1, h), name in zip(S.BUILDINGS, S.BUILD_NAMES):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="#555555",
                           edgecolor="k", alpha=0.95, zorder=3))
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"{name}\n{h} m", color="w",
            ha="center", va="center", fontsize=8, zorder=4)
# target & azimuth line
ax.plot([S.CENTER[0], tx], [S.CENTER[1], ty], "r--", lw=1.2, zorder=5,
        label="alert azimuth line")
ax.plot(tx, ty, marker="^", color="red", ms=12, zorder=6,
        label=f"T1 ({tx:.2f}, {ty:.2f}, H=45)")
ax.plot(*S.CENTER, "ko", ms=6, zorder=6, label="dispatch center O")
# optimal circle and station
th = np.linspace(0, 2 * np.pi, 400)
ax.plot(tx + r_theory * np.cos(th), ty + r_theory * np.sin(th), "w-.", lw=1.0,
        zorder=5, label=f"r* = {r_theory:.3f} m circle (el = 80 deg)")
ax.plot(bx, by, "*", color="gold", ms=18, markeredgecolor="k", zorder=7,
        label=f"optimal station ({bx:.2f}, {by:.2f})")
ax.set_xlim(0, 100); ax.set_ylim(0, 100)
ax.set_aspect("equal")
ax.set_xlabel("x (m) east"); ax.set_ylabel("y (m) north")
ax.set_title("Problem 1 — T1 vehicle-observable region and optimal station")
ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_p1_region.png", dpi=300)
print("figure saved:", f"{OUT}/fig_p1_region.png")

# ----------------------------------------------------------------------
# 4. Save results
# ----------------------------------------------------------------------
res = {
    "target": "T1", "target_xyz": [tx, ty, H],
    "grid_points": int(OK.size), "feasible_points": int(OK.sum()),
    "grid_optimum": {"r_m": r_grid, "x": float(X[i_min]), "y": float(Y[i_min])},
    "theoretical_r_min_m": r_theory,
    "optimal_station": {
        "x": round(bx, 4), "y": round(by, 4),
        "azimuth_from_target_deg": round(adeg, 2),
        "horizontal_range_m": round(r_theory, 4),
        "elevation_deg": round(bel, 4),
        "Q_per_m": round(Q_star, 6),
        "dist_to_center_m": round(bdO, 3),
    },
    "optimization_criterion": "maximize Q(s)=1/r(s) over the feasible set "
                              "(equivalently minimize horizontal range r); "
                              "tie-break: (1) nearest to dispatch center O, "
                              "(2) smallest azimuth CCW from east",
}
with open(f"{OUT}/problem1_results.json", "w", encoding="utf-8") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
print("results saved:", f"{OUT}/problem1_results.json")
