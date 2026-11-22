# -*- coding: utf-8 -*-
"""
Problem 1 (revised, v1.01) — Observable region and optimal station for T1.

What changed vs 1.py
--------------------
1. Figure: the observable region is now a true two-colour mask.  The old
   single-colour ListedColormap painted the whole field light blue, because
   both 0 and 1 mapped to the same colour.
2. Optimal station: the selected point is no longer the raw "nearest to O"
   boundary point (clearance 1.5004 m, i.e. only 0.4 mm above the 1.5 m
   vehicle buffer).  We require a user-set clearance margin (default 0.10 m)
   and then pick the point nearest to O, then the shortest gimbal slew.
   Because the optimal set is an arc with r = r_min, moving along the arc
   costs no Q at all.
3. The optimal SET is reported explicitly: the feasible arc(s) on the
   r = r_min circle, plus the blocked sector caused by B2.
4. Naming: "station polar angle" (around T1') is separated from "gimbal pan"
   (the camera's pointing azimuth).  They differ by 180 deg.
5. Vectorised arc sweep; fallback to the grid optimum if the r_min circle
   has no feasible point (useful when this module is reused for T2/T3/P4).
6. Output paths are anchored to this file; richer JSON; two-panel figure.

Outputs
-------
results/problem1_results_v101.json
results/fig_p1_region_v101.png
"""

import json
import math
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Rectangle

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
import swb_common as S  # noqa: E402

OUT = BASE / "results"
OUT.mkdir(exist_ok=True)

# ----------------------------------------------------------------------
# 0. Parameters of this run
# ----------------------------------------------------------------------
GRID_STEP = 0.20          # m, raster used for the region plot / area estimate
ARC_STEP = 0.02           # deg, sweep on the r = r_min circle
MARGIN = 0.10             # m, required clearance beyond the 1.5 m buffer
PAN_MODE = os.environ.get("PAN_MODE", "shortest")  # "shortest" or "ccw"
DESIGNATED_VEHICLE = "V4"  # P1 has no assignment; V4 is the closest/fastest vehicle

tgt = S.TARGETS["T1"]
tx, ty, H = tgt
r_theory = S.min_range_for_target(H)
Q_star = 1.0 / r_theory
planner = S.PathPlanner()   # geodesic travel times around the 1.5 m-inflated buildings


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------
def clearance_vec(xs, ys):
    """Distance by which the vehicle disk clears the field boundary and all
    buildings, i.e. min(distance_to_obstacle) - R_V.  >= 0 means feasible."""
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    c = np.minimum.reduce([
        xs - S.R_V,
        S.FIELD - S.R_V - xs,
        ys - S.R_V,
        S.FIELD - S.R_V - ys,
    ])
    for (x0, x1, y0, y1, _h) in S.BUILDINGS:
        c = np.minimum(c, S.point_rect_dist_vec(xs, ys, (x0, x1, y0, y1)) - S.R_V)
    return c


def gimbal_tau(xs, ys):
    """Gimbal pointing time from due-east/horizontal, vectorised."""
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    r = np.hypot(tx - xs, ty - ys)
    pan = (np.degrees(np.arctan2(ty - ys, tx - xs)) + 360.0) % 360.0
    tilt = np.degrees(np.arctan2(H - S.CAM_H, r))
    if PAN_MODE == "ccw":
        tau = np.maximum(pan, tilt) / S.SLEW_RATE
    else:
        tau = np.maximum(np.abs((pan + 180.0) % 360.0 - 180.0), tilt) / S.SLEW_RATE
    return pan, tilt, tau


def station_info(name, x, y):
    """Scalar metrics for one station (for JSON tables)."""
    x, y = float(x), float(y)
    r = math.hypot(tx - x, ty - y)
    el = math.degrees(math.atan2(H - S.CAM_H, r))
    pan = (math.degrees(math.atan2(ty - y, tx - x)) + 360.0) % 360.0
    tilt = el
    if PAN_MODE == "ccw":
        tau = max(pan, tilt) / S.SLEW_RATE
    else:
        tau = max(abs((pan + 180.0) % 360.0 - 180.0), tilt) / S.SLEW_RATE
    d_build = min(S.point_rect_dist(x, y, b[:4]) for b in S.BUILDINGS)
    clr = min(
        min(x - S.R_V, S.FIELD - S.R_V - x, y - S.R_V, S.FIELD - S.R_V - y),
        d_build - S.R_V,
    )
    polar = (math.degrees(math.atan2(y - ty, x - tx)) + 360.0) % 360.0
    feasible, _, _ = S.station_feasible(x, y, tgt)

    # total time from each vehicle's start to task end:
    # travel + pointing + 5 s observation + reset (reset = pointing)
    vehicle_times = {}
    for i, v in enumerate(S.V_NAMES):
        tt = planner.time(S.V_INIT[i], (x, y))
        vehicle_times[v] = round(tt + 2.0 * tau + S.OBS_DUR, 3)
    best_v = min(vehicle_times, key=vehicle_times.get)
    best_t = vehicle_times[best_v]
    best_travel = planner.time(S.V_INIT[S.V_NAMES.index(best_v)], (x, y))

    return {
        "name": name,
        "x": round(x, 4),
        "y": round(y, 4),
        "station_polar_angle_deg": round(polar, 3),
        "gimbal_pan_deg": round(pan, 3),
        "elevation_deg": round(el, 4),
        "range_m": round(r, 4),
        "Q_per_m": round(1.0 / r, 6) if r > 0 else None,
        "clearance_m": round(clr, 4),
        "distance_to_nearest_building_m": round(d_build, 4),
        "slew_time_s": round(tau, 3),
        "best_vehicle": best_v,
        "mission_time_s": best_t,
        "travel_time_s": round(best_travel, 3),
        "pointing_time_s": round(tau, 3),
        "observation_time_s": S.OBS_DUR,
        "reset_time_s": round(tau, 3),
        "vehicle_times_s": vehicle_times,
        "distance_to_center_m": round(math.dist((x, y), S.CENTER), 3),
        "feasible": bool(feasible),
    }


def feasible_intervals(angles, ok, step):
    """Contiguous feasible intervals on a 0..360 circle.

    Intervals that touch 0/360 are merged and returned unwrapped, e.g.
    [(236.9, 520.68)] for a single arc crossing due east.
    """
    a = angles[ok]
    if a.size == 0:
        return []
    gaps = np.where(np.diff(a) > 2.5 * step)[0]
    parts = list(np.split(a, gaps + 1))
    if (len(parts) >= 2 and parts[0][0] <= 2.5 * step
            and parts[-1][-1] >= 360.0 - 2.5 * step):
        merged = np.concatenate([parts[-1], parts[0] + 360.0])
        parts = [merged] + parts[1:-1]
    return [(float(p[0]), float(p[-1])) for p in parts]


def normalise_interval(a0, a1):
    """Shift an unwrapped interval into a readable range."""
    while a0 >= 180.0:
        a0 -= 360.0
        a1 -= 360.0
    while a0 < -180.0:
        a0 += 360.0
        a1 += 360.0
    return float(a0), float(a1)


# ----------------------------------------------------------------------
# 2. Observable region on a raster
# ----------------------------------------------------------------------
xs = np.arange(S.R_V, S.FIELD - S.R_V + 1e-9, GRID_STEP)
ys = np.arange(S.R_V, S.FIELD - S.R_V + 1e-9, GRID_STEP)
X, Y = np.meshgrid(xs, ys)

ok, el, rng = S.feasibility_mask(X.ravel(), Y.ravel(), tgt)
OK = ok.reshape(X.shape)
R = rng.reshape(X.shape)
Q = np.where(OK, 1.0 / np.maximum(R, 1e-12), np.nan)

n_grid = int(OK.size)
n_feas = int(OK.sum())
frac = n_feas / n_grid
print(f"T1 = ({tx:.4f}, {ty:.4f}, H={H})")
print(f"grid {GRID_STEP} m: {n_grid} points, feasible {n_feas} ({100*frac:.2f}%)")

# grid optimum (sanity check for the analytic value)
feas_r = np.where(OK, R, np.inf)
i_min = np.unravel_index(np.argmin(feas_r), feas_r.shape)
grid_r = float(feas_r[i_min])
grid_x, grid_y = float(X[i_min]), float(Y[i_min])
grid_el = float(el.reshape(X.shape)[i_min])
print(f"grid optimum: r = {grid_r:.4f} m at ({grid_x:.2f}, {grid_y:.2f}), "
      f"elevation = {grid_el:.3f} deg")
print(f"theoretical lower bound: r_min = {r_theory:.6f} m, Q* = {Q_star:.6f} 1/m")

# ----------------------------------------------------------------------
# 3. Feasible arc on the r = r_min circle (vectorised)
# ----------------------------------------------------------------------
ang = np.arange(0.0, 360.0, ARC_STEP)
ax_ = tx + r_theory * np.cos(np.deg2rad(ang))
ay_ = ty + r_theory * np.sin(np.deg2rad(ang))
ok_arc, el_arc, r_arc = S.feasibility_mask(ax_, ay_, tgt)
clr_arc = clearance_vec(ax_, ay_)

print(f"r_min circle: {int(ok_arc.sum())}/{ang.size} sampled angles feasible")

intervals = feasible_intervals(ang, ok_arc, ARC_STEP)
normalised = [normalise_interval(a, b) for a, b in intervals]
print("optimal-set intervals (unwrapped): "
      + ", ".join(f"[{a:.2f}, {b:.2f}]" for a, b in intervals))
print("optimal-set intervals (normalised): "
      + ", ".join(f"[{a:.2f}, {b:.2f}]" for a, b in normalised))

# blocked sector(s) = complement of the feasible arc
blocked = []
for k, (a0, a1) in enumerate(intervals):
    nxt = intervals[(k + 1) % len(intervals)][0]
    if k == len(intervals) - 1:
        nxt += 360.0
    if nxt > a1 + 1e-9:
        blocked.append((a1 % 360.0, nxt % 360.0))
print("blocked sector(s): "
      + ", ".join(f"({a:.2f}, {b:.2f})" for a, b in blocked))

# ----------------------------------------------------------------------
# 4. Tie-break: shortest total mission time from a standby vehicle
#    T_i(p) = travel_i(p) + tau(p) + 5 s observation + tau(p) reset
#    T(p)   = min_i T_i(p)   (the fastest vehicle that can do the task)
# ----------------------------------------------------------------------
def mission_terms(x, y):
    """(best total time, vehicle, travel time, pointing time) at one station."""
    _pan, _tilt, _tau = gimbal_tau(np.array([x]), np.array([y]))
    tau = float(_tau[0])
    best = None
    for i, v in enumerate(S.V_NAMES):
        trav = planner.time(S.V_INIT[i], (x, y))
        total = trav + 2.0 * tau + S.OBS_DUR
        if best is None or total < best[0]:
            best = (total, v, trav)
    return best[0], best[1], best[2], tau


def arc_time_table(angles_deg, margin):
    """Evaluate T(p) at all angles on the r = r_min circle."""
    a = np.asarray(angles_deg, float)
    x = tx + r_theory * np.cos(np.deg2rad(a))
    y = ty + r_theory * np.sin(np.deg2rad(a))
    ok, el, _rr = S.feasibility_mask(x, y, tgt)
    clr = clearance_vec(x, y)
    pan, tilt, tau = gimbal_tau(x, y)
    dO = np.hypot(x - S.CENTER[0], y - S.CENTER[1])
    T = np.full(a.shape, np.inf)
    veh = np.array([""] * a.size, dtype=object)
    trav = np.full(a.shape, np.inf)
    valid = ok & (clr >= margin - 1e-12)

    if DESIGNATED_VEHICLE is None:
        vehicles = list(enumerate(S.V_NAMES))
    else:
        vehicles = [(S.V_NAMES.index(DESIGNATED_VEHICLE), DESIGNATED_VEHICLE)]

    for i, v in vehicles:
        p0 = S.V_INIT[i]
        # Euclidean lower bound on the mission time; the path is never shorter.
        lb = (np.hypot(x - p0[0], y - p0[1]) / S.V_SPEED
              + 2.0 * tau + S.OBS_DUR)
        order = np.argsort(lb)
        best = np.inf
        for j in order:
            if not valid[j]:
                continue
            if lb[j] >= best - 1e-12:
                break
            tt = planner.time(p0, (float(x[j]), float(y[j])))
            total = tt + 2.0 * float(tau[j]) + S.OBS_DUR
            if total < T[j]:
                T[j], veh[j], trav[j] = total, v, tt
            if total < best:
                best = total

    return dict(angles=a, x=x, y=y, ok=ok, el=el, clr=clr,
                pan=pan, tilt=tilt, tau=tau, dO=dO, T=T, veh=veh, trav=trav)


def pick_station(table, j, margin):
    return {
        "index": int(j),
        "x": float(table["x"][j]),
        "y": float(table["y"][j]),
        "r": float(r_theory),
        "el": float(table["el"][j]),
        "clr": float(table["clr"][j]),
        "pan": float(table["pan"][j]),
        "tilt": float(table["tilt"][j]),
        "tau": float(table["tau"][j]),
        "dO": float(table["dO"][j]),
        "T": float(table["T"][j]),
        "veh": str(table["veh"][j]),
        "trav": float(table["trav"][j]),
        "margin": margin,
    }


used_margin = MARGIN
coarse = arc_time_table(np.arange(0.0, 360.0, 0.05), used_margin)
if not np.isfinite(coarse["T"]).any():
    used_margin = 0.0
    coarse = arc_time_table(np.arange(0.0, 360.0, 0.05), used_margin)
    print("WARNING: no arc point with the requested margin; relaxed to 0.")

cidx = np.where(np.isfinite(coarse["T"]))[0]
corder = np.lexsort((coarse["angles"][cidx],
                     -coarse["clr"][cidx],
                     coarse["T"][cidx]))
j0 = int(cidx[corder[0]])

# local refinement around the coarse optimum
a0 = float(coarse["angles"][j0])
fine = arc_time_table(np.arange(a0 - 0.05, a0 + 0.05 + 1e-12, 0.002),
                      used_margin)
fidx = np.where(np.isfinite(fine["T"]))[0]
if fidx.size:
    forder = np.lexsort((fine["angles"][fidx],
                         -fine["clr"][fidx],
                         fine["T"][fidx]))
    sel = pick_station(fine, int(fidx[forder[0]]), used_margin)
else:
    sel = pick_station(coarse, j0, used_margin)

print(f"fastest mission-time station: T = {sel['T']:.3f} s with {sel['veh']} "
      f"(travel {sel['trav']:.3f} + 2*tau {2*sel['tau']:.3f} + {S.OBS_DUR:g} s)")
print(f"  at polar angle "
      f"{(np.degrees(np.arctan2(sel['y'] - ty, sel['x'] - tx)) % 360):.3f} deg, "
      f"clearance {sel['clr']:.3f} m")

fallback = not np.isfinite(sel["T"])
if fallback:
    # No feasible point on the r_min circle (should not happen for T1).
    # Fall back to the best feasible grid point.
    bx, by = grid_x, grid_y
    br, bel = grid_r, grid_el
    bclr = float(clearance_vec(np.array([bx]), np.array([by]))[0])
    _pan, _tilt, _tau = gimbal_tau(np.array([bx]), np.array([by]))
    bpan, btilt, btau = float(_pan[0]), float(_tilt[0]), float(_tau[0])
    bdO = float(math.dist((bx, by), S.CENTER))
    bT, bveh, btrav, _ = mission_terms(bx, by)
    bmargin = 0.0
    print("WARNING: r_min circle has no feasible point; using grid optimum.")
else:
    bx, by = sel["x"], sel["y"]
    br, bel = sel["r"], sel["el"]
    bclr, bpan, btilt = sel["clr"], sel["pan"], sel["tilt"]
    btau, bdO = sel["tau"], sel["dO"]
    bT, bveh, btrav = sel["T"], sel["veh"], sel["trav"]
    bmargin = sel["margin"]

# scalar verification
feas_sel, el_sel, r_sel = S.station_feasible(bx, by, tgt)
los_clear = not S.los_blocked(bx, by, tx, ty, H)
print(f"selected station: ({bx:.4f}, {by:.4f})")
print(f"  range {br:.6f} m | elevation {bel:.6f} deg | Q {1/br:.6f} 1/m")
print(f"  clearance beyond 1.5 m buffer: {bclr:.6f} m")
print(f"  gimbal pan {bpan:.3f} deg | slew {btau:.3f} s | distance to O {bdO:.3f} m")
print(f"  fastest mission time {bT:.3f} s via {bveh} (travel {btrav:.3f} s)")
print(f"  feasible={feas_sel}, LOS clear={los_clear}")

assert abs(br - r_theory) < 1e-6
assert abs(bel - S.EL_MAX) < 1e-6
assert feas_sel and los_clear
assert bclr >= bmargin - 1e-9

# ----------------------------------------------------------------------
# 5. Alternative stations (for the paper's tie-break table)
# ----------------------------------------------------------------------
def at_angle(a_deg):
    a = math.radians(a_deg)
    return tx + r_theory * math.cos(a), ty + r_theory * math.sin(a)


alternatives = []
if intervals:
    sw_a = intervals[0][0]
    for name, a_deg in (
        ("SW arc endpoint", sw_a),
        ("N", 90.0),
        ("E", 0.0),
        ("S", 270.0),
    ):
        x, y = at_angle(a_deg)
        alternatives.append(station_info(name, x, y))
alternatives.append(station_info("selected (fastest time)", bx, by))

# the unconstrained nearest-O direction is inside B2; keep it as a caveat
a_O = math.degrees(math.atan2(S.CENTER[1] - ty, S.CENTER[0] - tx)) % 360.0
xO, yO = at_angle(a_O)
unconstrained = station_info("nearest-O direction (infeasible)", xO, yO)

# ----------------------------------------------------------------------
# 6. Figure
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 7.2), dpi=300)

# ---------- (a) full field ----------
ax1.contourf(X, Y, OK.astype(float), levels=[0.5, 1.5],
             colors=["#cfe3ff"], alpha=0.95, zorder=0)

for (x0, x1, y0, y1, h), name in zip(S.BUILDINGS, S.BUILD_NAMES):
    ax1.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                            facecolor="#555555", edgecolor="k", zorder=3))
    ax1.text((x0 + x1) / 2, (y0 + y1) / 2, f"{name}\n{h:g} m",
             color="w", ha="center", va="center", fontsize=8, zorder=4)

ax1.plot([S.CENTER[0], tx], [S.CENTER[1], ty], "r--", lw=1.1, zorder=2,
         label="alert azimuth (not a LOS)")
ax1.plot(tx, ty, marker="^", color="red", ms=12, zorder=6,
         label=f"T1 ({tx:.2f}, {ty:.2f}, H={H:g})")
ax1.plot(*S.CENTER, "ko", ms=6, zorder=6, label="dispatch centre O")

theta = np.linspace(0, 2 * np.pi, 720)
ax1.plot(tx + r_theory * np.cos(theta), ty + r_theory * np.sin(theta),
         color="0.35", ls="-.", lw=1.0, zorder=4, label=f"r = r_min = {r_theory:.3f} m")

for a0, a1 in intervals:
    aa = np.deg2rad(np.linspace(a0, a1, 900))
    ax1.plot(tx + r_theory * np.cos(aa), ty + r_theory * np.sin(aa),
             color="gold", lw=3.2, solid_capstyle="round", zorder=6,
             label="optimal set S*" if a0 == intervals[0][0] else None)

# standby vehicles at their initial positions (blue squares, as in the problem figure)
for i, (vx, vy) in enumerate(S.V_INIT):
    ax1.plot(vx, vy, "s", color="#1f4fd8", ms=6.0, mec="k", mew=0.4, zorder=7)
    ax1.text(vx + 1.2, vy + 1.2, S.V_NAMES[i], fontsize=7,
             color="#1f4fd8", zorder=8)

# selected station: same blue-square vehicle icon, gold edge,
# plus the physical 1.5 m vehicle disk (dotted circle)
ax1.plot(bx, by, "s", color="#1f4fd8", ms=13, mec="gold", mew=2.0, zorder=9)
ax1.add_patch(Circle((bx, by), S.R_V, fill=False, ec="k", ls=":",
                     lw=1.0, zorder=8))
ax1.text(bx + 1.0, by - 2.0, "V4 station", fontsize=7,
         color="#1f4fd8", zorder=9)

ax1.add_patch(Rectangle((72, 55), 28, 30, fill=False, ls="--", ec="k",
                        lw=1.0, zorder=7))
ax1.text(72.6, 56.2, "detail shown in (b)", fontsize=7.2, zorder=8,
         bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.8))

ax1.set_xlim(0, 100)
ax1.set_ylim(0, 100)
ax1.set_aspect("equal")
ax1.set_xlabel("x (m) east")
ax1.set_ylabel("y (m) north")
ax1.set_title(f"(a) T1 observable region ({GRID_STEP:g} m grid; "
              f"{100*frac:.1f}% of sampled centers)")
handles = [
    Patch(facecolor="#cfe3ff", label="observable region"),
    Patch(facecolor="white", edgecolor="0.7", label="infeasible"),
    Patch(facecolor="#555555", label="building"),
    Line2D([0], [0], color="gold", lw=3, label="optimal set S*"),
    Line2D([0], [0], marker="s", color="none", mfc="#1f4fd8", mec="k",
           ms=7, label="standby vehicle V1–V5"),
    Line2D([0], [0], marker="s", color="none", mfc="#1f4fd8", mec="gold",
           mew=1.5, ms=11, label="selected station (V4)"),
]
labels = [h.get_label() for h in handles]
# legend is placed at figure level below the panels (see save section)

# ---------- (b) zoom near T1 ----------
Qm = np.ma.masked_invalid(Q)
cf = ax2.contourf(X, Y, Qm, levels=14, cmap="viridis", alpha=0.88, zorder=0)
cbar = fig.colorbar(cf, ax=ax2, shrink=0.86, pad=0.02)
cbar.set_label("observation quality  Q = 1/r  (1/m)")

for (x0, x1, y0, y1, h), name in zip(S.BUILDINGS, S.BUILD_NAMES):
    ax2.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                            facecolor="#555555", edgecolor="k", zorder=3))
    if x1 >= 72:  # only B2 is visible in this zoom
        ax2.text(max((x0 + x1) / 2, 73.0), (y0 + y1) / 2, f"{name}\n{h:g} m",
                 color="w", ha="center", va="center", fontsize=8, zorder=4)

ax2.plot(tx + r_theory * np.cos(theta), ty + r_theory * np.sin(theta),
         color="0.35", ls="-.", lw=1.0, zorder=4)
for a0, a1 in intervals:
    aa = np.deg2rad(np.linspace(a0, a1, 900))
    ax2.plot(tx + r_theory * np.cos(aa), ty + r_theory * np.sin(aa),
             color="gold", lw=3.2, solid_capstyle="round", zorder=6)
ax2.plot(tx, ty, marker="^", color="red", ms=11, zorder=9)
ax2.text(tx - 1.3, ty + 1.0, "T1", fontsize=8.5, color="red",
         ha="right", zorder=10)

# dashed radius from T1 to the r_min circle (eastward, clear of buildings/labels)
ax2.plot([tx, tx + r_theory], [ty, ty], color="k", ls="--", lw=1.1, zorder=7)
ax2.text(tx + 2.6, ty + 0.2, f"$r_{{\\min}}$ = {r_theory:.3f} m",
         fontsize=8.0, color="k", ha="left", va="bottom", zorder=10,
         bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))
v4x, v4y = S.V_INIT[S.V_NAMES.index(DESIGNATED_VEHICLE)]
ax2.annotate("", xy=(bx, by), xytext=(v4x, v4y),
             arrowprops=dict(arrowstyle="->", ls="--", lw=1.0, color="#1f4fd8"),
             zorder=7)
ax2.plot(v4x, v4y, "s", color="#1f4fd8", ms=6.5, mec="k", mew=0.4, zorder=8)
ax2.text(v4x + 0.7, v4y + 0.7, "V4 start", fontsize=7,
         color="#1f4fd8", zorder=9)
ax2.plot(bx, by, "s", color="#1f4fd8", ms=12, mec="gold", mew=2.0, zorder=9)
ax2.add_patch(Circle((bx, by), S.R_V, fill=False, ec="k", ls=":",
                     lw=1.2, zorder=8))

ax2.annotate(f"({bx:.2f}, {by:.2f})\nclearance {bclr:.2f} m\n"
             f"T = {bT:.2f} s ({bveh})",
             xy=(bx, by), xytext=(86.2, 57.6),
             fontsize=8.5, zorder=9, ha="left", va="bottom",
             bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.4", alpha=0.9),
             arrowprops=dict(arrowstyle="->", lw=0.9, color="k"))


ax2.set_xlim(72, 100)
ax2.set_ylim(55, 85)
ax2.set_aspect("equal")
ax2.set_xlabel("x (m) east")
ax2.set_ylabel("y (m) north")
ax2.set_title("(b) zoom near T1: Q field and optimal arc", fontsize=10)

fig.suptitle("Problem 1 — observable region and optimal station for T1",
             fontsize=13, y=0.975)
fig.subplots_adjust(left=0.055, right=0.905, top=0.865, bottom=0.155, wspace=0.22)
fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False,
           fontsize=7.5, bbox_to_anchor=(0.5, 0.014))
fig.savefig(OUT / "fig_p1_region_v101.png", dpi=300)
print("figure saved:", OUT / "fig_p1_region_v101.png")

# ----------------------------------------------------------------------
# 7. JSON
# ----------------------------------------------------------------------
selected_info = station_info("selected (fastest time)", bx, by)
selected_info["clearance_margin_used_m"] = round(bmargin, 4)
selected_info["vehicle_disk_radius_m"] = S.R_V

res = {
    "target": "T1",
    "target_xyz": [round(tx, 6), round(ty, 6), H],
    "pan_mode": PAN_MODE,
    "designated_vehicle": DESIGNATED_VEHICLE,
    "grid": {
        "step_m": GRID_STEP,
        "points": n_grid,
        "feasible_points": n_feas,
        "feasible_fraction": round(frac, 6),
    },
    "theoretical_r_min_m": round(r_theory, 6),
    "Q_star_per_m": round(Q_star, 6),
    "grid_optimum": {
        "r_m": round(grid_r, 6),
        "x": round(grid_x, 4),
        "y": round(grid_y, 4),
        "elevation_deg": round(grid_el, 4),
    },
    "optimal_set": {
        "description": "Q = 1/r depends only on range, so every feasible point "
                       "on the r = r_min circle is optimal; the optimum is an arc, "
                       "not a single point.",
        "intervals_deg_unwrapped": [[round(a, 3), round(b, 3)] for a, b in intervals],
        "intervals_deg_normalised": [[round(a, 3), round(b, 3)] for a, b in normalised],
        "blocked_sectors_deg": [[round(a, 3), round(b, 3)] for a, b in blocked],
        "arc_length_deg": round(sum(b - a for a, b in intervals), 3),
    },
    "selected_station": selected_info,
    "unconstrained_nearest_O_direction": unconstrained,
    "alternatives": alternatives,
    "clearance_definition": "min(field-boundary margin, distance to nearest "
                            "building - 1.5 m); 0 means the vehicle disk touches "
                            "the 1.5 m buffer boundary.",
    "time_model": "P1 has no vehicle assignment; V4 is the closest/fastest "
                  "standby vehicle and is treated as the designated observer. "
                  "Total time from V4's start to task end = geodesic travel + "
                  "gimbal pointing + 5 s observation + gimbal reset; reset time "
                  "equals pointing time. tau uses the shortest pan direction "
                  "unless pan_mode='ccw'.",
    "optimization_criterion": "maximise Q(p)=1/r(p) over the feasible set; "
                              "equivalent to minimise horizontal range r.",
    "tie_breaking": "within the Q-optimal arc: require clearance >= 0.10 m beyond "
                    "the 1.5 m vehicle buffer, then minimise the total mission "
                    "time T(alpha)=d_V4(alpha)/3 + 2*tau(alpha) + 5 s over the "
                    "arc angle alpha; ties by larger clearance, then smaller "
                    "distance to the dispatch centre O, then polar angle. The "
                    "search uses a Euclidean lower bound plus exact-path "
                    "branch-and-bound, so not every arc point is path-planned.",
    "assumptions": [
        "Vehicle disk must stay inside the field and outside building footprints.",
        "Camera is at the disk centre, 2 m above ground.",
        "Each vehicle starts at its initial position with the gimbal due east and "
        "horizontal, so pointing time = tau(p) and reset time = tau(p).",
        "LOS tangency/contact with a building counts as occlusion.",
        "The vehicle disk may touch (not enter) the 1.5 m buffer boundary; a "
        "0.10 m margin is imposed for the reported station.",
    ],
}

with open(OUT / "problem1_results_v101.json", "w", encoding="utf-8") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
print("results saved:", OUT / "problem1_results_v101.json")
