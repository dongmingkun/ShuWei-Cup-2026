# -*- coding: utf-8 -*-
"""Q1: 可观测域 + r80 临界圆最优弧 + 任务周期决胜选站 + 独立验证.

运行: python -X utf8 -B q1_observable_domain.py
口径与 geometry.py 一致: 净距 >= 1.5 m 贴边允许; 俯仰 [-10°,80°] 闭区间;
视线切触即遮挡. 论文采用口径: 最优弧上先要求 0.10 m 硬净距裕量, 再最小化
指定车辆任务周期 = 绕障行驶 + 转位 + 5 s 观测 + 复位; 平局取更大裕量、
更近调度中心 O、更小方位角. 指定车辆为距目标投影最近的待命车.
站点坐标输出 6 位小数, 避免 4 位舍入把 80° 边界点打印成不可行.
"""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Rectangle

from geometry import (BUILDINGS, DISPATCH, FIELD_SIZE, SPEED, TARGETS,
                      VEHICLE_INIT, VEHICLE_RADIUS, clearance_ok,
                      dist_point_to_rect, elev_min_radius, elevation_angle,
                      is_station_feasible, los_ok, travel_time)
from simulator import OBS_DURATION, pointing, rotate_duration

OUT = Path(__file__).with_name("q1_output")
OUT.mkdir(exist_ok=True)

CLEARANCE_MARGIN = 0.10
ARC_STEP = 0.01
COARSE_STEP = 0.05
FINE_STEP = 0.002
INTERIOR_STEP = 0.5


def _los_blocked_vec(X, Y, Tx, Ty, H, b):
    dx, dy, dz = Tx - X, Ty - Y, H - 2.0
    z0 = np.full_like(X, 2.0)
    tmin = np.zeros_like(X)
    tmax = np.ones_like(X)
    parallel_inside = np.ones(X.shape, dtype=bool)
    for d, lo, hi, p in ((dx, b.x0, b.x1, X), (dy, b.y0, b.y1, Y),
                         (dz, 0.0, b.h, z0)):
        small = np.abs(d) < 1e-15
        parallel_inside &= ~small | ((p >= lo) & (p <= hi))
        safe = np.where(small, 1.0, d)
        t1, t2 = (lo - p) / safe, (hi - p) / safe
        tmin = np.maximum(tmin, np.where(small, -np.inf, np.minimum(t1, t2)))
        tmax = np.minimum(tmax, np.where(small, np.inf, np.maximum(t1, t2)))
    return parallel_inside & (tmin <= tmax)


def clearance_margin_vec(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    c = np.minimum.reduce([x - VEHICLE_RADIUS, FIELD_SIZE - VEHICLE_RADIUS - x,
                           y - VEHICLE_RADIUS, FIELD_SIZE - VEHICLE_RADIUS - y])
    for b in BUILDINGS:
        ddx = np.maximum(np.maximum(b.x0 - x, x - b.x1), 0.0)
        ddy = np.maximum(np.maximum(b.y0 - y, y - b.y1), 0.0)
        c = np.minimum(c, np.hypot(ddx, ddy) - VEHICLE_RADIUS)
    return c


def clearance_margin(p):
    x, y = p
    return float(min(x - VEHICLE_RADIUS, FIELD_SIZE - VEHICLE_RADIUS - x,
                     y - VEHICLE_RADIUS, FIELD_SIZE - VEHICLE_RADIUS - y,
                     min(dist_point_to_rect(x, y, b.x0, b.x1, b.y0, b.y1)
                         - VEHICLE_RADIUS for b in BUILDINGS)))


def feasible_vec(x, y, target):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    r = np.hypot(x - target.x, y - target.y)
    ok = ((x >= VEHICLE_RADIUS - 1e-9) & (x <= FIELD_SIZE - VEHICLE_RADIUS + 1e-9)
          & (y >= VEHICLE_RADIUS - 1e-9) & (y <= FIELD_SIZE - VEHICLE_RADIUS + 1e-9))
    ok &= r >= elev_min_radius(target) - 1e-9
    ok &= clearance_margin_vec(x, y) >= -1e-9
    for b in BUILDINGS:
        ok &= ~_los_blocked_vec(x, y, target.x, target.y, target.H, b)
    return ok, r


def region_mask(target, res):
    n = int(round(FIELD_SIZE / res)) + 1
    xs = np.linspace(0.0, FIELD_SIZE, n)
    X, Y = np.meshgrid(xs, xs)
    ok, _ = feasible_vec(X, Y, target)
    return xs, ok


def region_stats(target, res):
    from scipy import ndimage
    xs, m = region_mask(target, res)
    labels, ncomp = ndimage.label(m)
    sizes = ndimage.sum_labels(m, labels, index=np.arange(1, ncomp + 1))
    return {
        "res_m": res,
        "area_m2": float(m.sum()) * res * res,
        "n_components": int(ncomp),
        "largest_component_area_m2": float(sizes.max()) * res * res if ncomp else 0.0,
        "mask": m, "xs": xs,
    }


def validate_against_scalar(target, res, n_samples=400, seed=2026):
    rng = np.random.default_rng(seed)
    xs, m = region_mask(target, res)
    ok = True
    for _ in range(n_samples):
        p = (float(rng.uniform(0, FIELD_SIZE)), float(rng.uniform(0, FIELD_SIZE)))
        i = min(max(int(round(p[0] / res)), 0), len(xs) - 1)
        j = min(max(int(round(p[1] / res)), 0), len(xs) - 1)
        p_snap = (float(xs[i]), float(xs[j]))
        if bool(m[j, i]) != is_station_feasible(p_snap, target):
            ok = False
            print(f"  MISMATCH at {p_snap}")
    return ok


def designated_vehicle(target):
    return min(VEHICLE_INIT,
               key=lambda v: math.dist(VEHICLE_INIT[v], (target.x, target.y)))


def mission_at(p, target, vehicle):
    pan, tilt = pointing(p, target)
    slew = rotate_duration(0.0, 0.0, pan, tilt)
    reset = rotate_duration(pan, tilt, 0.0, 0.0)
    tt = travel_time(VEHICLE_INIT[vehicle], p)
    if tt is None:
        return None
    return {"pan_deg": pan, "tilt_deg": tilt, "slew_s": slew, "reset_s": reset,
            "travel_time_s": tt, "travel_distance_m": tt * SPEED,
            "mission_cycle_s": tt + slew + OBS_DURATION + reset}


def _circle_points(target, r80, angles):
    a = np.radians(np.asarray(angles, float))
    return target.x + r80 * np.cos(a), target.y + r80 * np.sin(a)


def feasible_arcs(target, r80, step_deg=ARC_STEP):
    n = int(round(360.0 / step_deg))
    angles = np.arange(n) * step_deg
    x, y = _circle_points(target, r80, angles)
    ok, _ = feasible_vec(x, y, target)
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return [], False
    if idx.size == n:
        return [(0.0, 360.0)], True
    runs, start, prev = [], idx[0], idx[0]
    for k in idx[1:]:
        if k != prev + 1:
            runs.append((start, prev))
            start = k
        prev = k
    runs.append((start, prev))
    if len(runs) > 1 and runs[0][0] == 0 and runs[-1][1] == n - 1:
        s_last, e_last = runs.pop()
        s_first, e_first = runs.pop(0)
        runs.append((s_last - n, e_first))
        runs.sort()
    return [(float(s * step_deg), float(e * step_deg)) for s, e in runs], False


def arc_endpoint_causes(target, r80, arcs, full_circle):
    if full_circle:
        return []
    causes = []
    for a0, a1 in arcs:
        for a in (a0 - ARC_STEP, a1 + ARC_STEP):
            p = _circle_points(target, r80, [a % 360.0])
            p = (float(p[0][0]), float(p[1][0]))
            reason = []
            if not clearance_ok(p):
                reason.append("clearance")
            e = elevation_angle(p, target)
            if e > 80.0 + 1e-9 or e < -10.0 - 1e-9:
                reason.append("elevation")
            if not is_station_feasible(p, target) and not reason:
                reason.append("los")
            causes.append({"angle_deg": round(float(a % 360.0), 2),
                           "blocked_by": reason})
    return causes


def arc_time_table(target, r80, angles, margin, vehicle):
    a = np.asarray(angles, float)
    x, y = _circle_points(target, r80, a)
    ok, r = feasible_vec(x, y, target)
    clr = clearance_margin_vec(x, y)
    pan = (np.degrees(np.arctan2(target.y - y, target.x - x)) + 360.0) % 360.0
    tilt = np.degrees(np.arctan2(target.H - 2.0, r))
    dpan = np.abs((pan + 180.0) % 360.0 - 180.0)
    tau = np.maximum(dpan, tilt) / 30.0
    dO = np.hypot(x - DISPATCH[0], y - DISPATCH[1])
    valid = ok & (clr >= margin - 1e-12)
    p0 = VEHICLE_INIT[vehicle]
    lb = np.hypot(x - p0[0], y - p0[1]) / SPEED + 2.0 * tau + OBS_DURATION
    T = np.full(a.shape, np.inf)
    trav = np.full(a.shape, np.inf)
    best = np.inf
    for j in np.argsort(lb):
        if not valid[j]:
            continue
        if lb[j] >= best - 1e-12:
            break
        tt = travel_time(p0, (float(x[j]), float(y[j])))
        if tt is None:
            continue
        total = tt + 2.0 * float(tau[j]) + OBS_DURATION
        T[j], trav[j] = total, tt
        if total < best:
            best = total
    return {"angles": a, "x": x, "y": y, "ok": ok, "clr": clr, "pan": pan,
            "tilt": tilt, "tau": tau, "dO": dO, "T": T, "trav": trav}


def _best_index(table):
    idx = np.flatnonzero(np.isfinite(table["T"]))
    if idx.size == 0:
        raise RuntimeError("r80 圆上不存在满足硬裕量的可行站")
    order = np.lexsort((table["angles"][idx], table["dO"][idx],
                        -table["clr"][idx], table["T"][idx]))
    return int(idx[order[0]])


def pick_station(target, r80, vehicle):
    coarse = arc_time_table(target, r80, np.arange(0.0, 360.0, COARSE_STEP),
                            CLEARANCE_MARGIN, vehicle)
    j0 = _best_index(coarse)
    a0 = float(coarse["angles"][j0])
    fine = arc_time_table(target, r80,
                          np.arange(a0 - COARSE_STEP, a0 + COARSE_STEP + 1e-12,
                                    FINE_STEP),
                          CLEARANCE_MARGIN, vehicle)
    j = _best_index(fine)
    p = (float(fine["x"][j]), float(fine["y"][j]))
    m = mission_at(p, target, vehicle)
    return {
        "x": p[0], "y": p[1],
        "azimuth_on_circle_deg": float(fine["angles"][j]) % 360.0,
        "pan_deg": m["pan_deg"], "tilt_deg": m["tilt_deg"],
        "elevation_angle_deg": elevation_angle(p, target),
        "range_m": math.dist(p, (target.x, target.y)),
        "clearance_margin_m": clearance_margin(p),
        "slew_s": m["slew_s"], "reset_s": m["reset_s"],
        "travel_time_s": m["travel_time_s"],
        "travel_distance_m": m["travel_distance_m"],
        "mission_cycle_s": m["mission_cycle_s"],
        "dist_to_O_m": math.dist(p, DISPATCH),
        "feasible": is_station_feasible(p, target),
    }


def nearest_O_alternative(target, r80, arcs, full_circle, vehicle):
    angles = np.arange(0.0, 360.0, ARC_STEP)
    x, y = _circle_points(target, r80, angles)
    ok, _ = feasible_vec(x, y, target)
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return None
    dO = np.hypot(x[idx] - DISPATCH[0], y[idx] - DISPATCH[1])
    j = int(idx[np.argmin(dO)])
    p = (float(x[j]), float(y[j]))
    m = mission_at(p, target, vehicle)
    return {"x": p[0], "y": p[1],
            "azimuth_on_circle_deg": float(angles[j]),
            "clearance_margin_m": clearance_margin(p),
            "dist_to_O_m": float(dO[np.argmin(dO)]),
            "mission_cycle_s": m["mission_cycle_s"] if m else None,
            "note": "earlier nearest-to-O tie-break without the 0.10 m margin"}


def plot_grid(target, res):
    xs = np.arange(VEHICLE_RADIUS, FIELD_SIZE - VEHICLE_RADIUS + 1e-9, res)
    X, Y = np.meshgrid(xs, xs)
    ok, r = feasible_vec(X, Y, target)
    Q = np.where(ok, 1.0 / np.maximum(r, 1e-12), np.nan)
    return X, Y, ok, Q


def _draw_field(ax, target, r80, arcs, station, vehicle, label_station=True):
    for k, b in enumerate(BUILDINGS, start=1):
        ax.add_patch(Rectangle((b.x0, b.y0), b.x1 - b.x0, b.y1 - b.y0,
                               facecolor="#555555", edgecolor="k", zorder=3))
        ax.text((b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2, f"B{k}\n{b.h:g} m",
                color="w", ha="center", va="center", fontsize=8, zorder=4)
    th = np.linspace(0, 2 * np.pi, 721)
    ax.plot(target.x + r80 * np.cos(th), target.y + r80 * np.sin(th),
            color="0.35", ls="-.", lw=1.0, zorder=4)
    for a0, a1 in arcs:
        aa = np.radians(np.linspace(a0, a1, 900))
        ax.plot(target.x + r80 * np.cos(aa), target.y + r80 * np.sin(aa),
                color="gold", lw=2.6, solid_capstyle="round", zorder=6)
    ax.plot(target.x, target.y, "^", color="red", ms=10, zorder=8)
    if station:
        ax.plot(station["x"], station["y"], "s", color="#1f4fd8", ms=12,
                mec="gold", mew=2.0, zorder=9)
        ax.add_patch(Circle((station["x"], station["y"]), VEHICLE_RADIUS,
                            fill=False, ec="k", ls=":", lw=1.0, zorder=8))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m) east")
    ax.set_ylabel("y (m) north")


def plot_t1(target, r80, arcs, station, vehicle):
    X, Y, ok, Q = plot_grid(target, 0.2)
    frac = float(ok.mean())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 7.2), dpi=300)
    ax1.contourf(X, Y, ok.astype(float), levels=[0.5, 1.5],
                 colors=["#cfe3ff"], alpha=0.95, zorder=0)
    _draw_field(ax1, target, r80, arcs, station, vehicle)
    ax1.plot([DISPATCH[0], target.x], [DISPATCH[1], target.y], "r--", lw=1.1,
             zorder=2)
    ax1.plot(*DISPATCH, "ko", ms=6, zorder=6)
    ax1.text(target.x + 1.2, target.y + 1.2,
             f"T1 ({target.x:.2f}, {target.y:.2f}, H={target.H:g})",
             color="red", fontsize=8, zorder=9)
    for v, p in VEHICLE_INIT.items():
        ax1.plot(*p, "s", color="#1f4fd8", ms=6, mec="k", mew=0.4, zorder=7)
        ax1.text(p[0] + 1.2, p[1] + 1.2, f"V{v}", fontsize=7,
                 color="#1f4fd8", zorder=8)
    ax1.text(station["x"] + 1.0, station["y"] - 2.0, f"V{vehicle} station",
             fontsize=7, color="#1f4fd8", zorder=9)
    ax1.add_patch(Rectangle((72, 55), 28, 30, fill=False, ls="--", ec="k",
                            lw=1.0, zorder=7))
    ax1.text(72.6, 56.2, "detail shown in (b)", fontsize=7.2, zorder=8,
             bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.8))
    ax1.set_title(f"(a) T1 observable region (0.2 m grid; "
                  f"{100 * frac:.1f}% of sampled centers)")

    cf = ax2.contourf(X, Y, np.ma.masked_invalid(Q), levels=14, cmap="viridis",
                      alpha=0.88, zorder=0)
    cbar = fig.colorbar(cf, ax=ax2, shrink=0.86, pad=0.02)
    cbar.set_label("observation quality  Q = 1/r  (1/m)")
    _draw_field(ax2, target, r80, arcs, station, vehicle)
    ax2.text(target.x - 1.3, target.y + 1.0, "T1", fontsize=8.5, color="red",
             ha="right", zorder=10)
    ax2.plot([target.x, target.x + r80], [target.y, target.y], "k--", lw=1.1,
             zorder=7)
    ax2.text(target.x + 2.6, target.y + 0.2, f"$r_{{min}}$ = {r80:.3f} m",
             fontsize=8, ha="left", va="bottom", zorder=10,
             bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))
    v0 = VEHICLE_INIT[vehicle]
    ax2.annotate("", xy=(station["x"], station["y"]), xytext=v0,
                 arrowprops=dict(arrowstyle="->", ls="--", lw=1.0, color="#1f4fd8"),
                 zorder=7)
    ax2.plot(*v0, "s", color="#1f4fd8", ms=6.5, mec="k", mew=0.4, zorder=8)
    ax2.text(v0[0] + 0.7, v0[1] + 0.7, f"V{vehicle} start", fontsize=7,
             color="#1f4fd8", zorder=9)
    ax2.annotate(f"({station['x']:.2f}, {station['y']:.2f})\n"
                 f"clearance {station['clearance_margin_m']:.2f} m\n"
                 f"T = {station['mission_cycle_s']:.2f} s (V{vehicle})",
                 xy=(station["x"], station["y"]), xytext=(86.2, 57.6),
                 fontsize=8.5, zorder=9, ha="left", va="bottom",
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.4", alpha=0.9),
                 arrowprops=dict(arrowstyle="->", lw=0.9, color="k"))
    ax2.set_xlim(72, 100)
    ax2.set_ylim(55, 85)
    ax2.set_title("(b) zoom near T1: Q field and optimal arc", fontsize=10)
    fig.suptitle("Problem 1 — observable region and optimal station for T1",
                 fontsize=13, y=0.975)
    handles = [
        Patch(facecolor="#cfe3ff", label="observable region"),
        Patch(facecolor="white", edgecolor="0.7", label="infeasible"),
        Patch(facecolor="#555555", label="building"),
        Line2D([0], [0], color="gold", lw=3, label="optimal set S*"),
        Line2D([0], [0], marker="s", color="none", mfc="#1f4fd8", mec="k",
               ms=7, label="standby vehicle V1–V5"),
        Line2D([0], [0], marker="s", color="none", mfc="#1f4fd8", mec="gold",
               mew=1.5, ms=11, label=f"selected station (V{vehicle})"),
    ]
    fig.subplots_adjust(left=0.055, right=0.905, top=0.865, bottom=0.155, wspace=0.22)
    fig.legend(handles, [h.get_label() for h in handles], loc="lower center",
               ncol=6, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, 0.014))
    fig.savefig(OUT / "q1_T1_region.png", dpi=300)
    plt.close(fig)
    return frac


def plot_t23(targets_info):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.6), dpi=300)
    fracs = {}
    for ax, (name, info) in zip(axes, targets_info.items()):
        target = TARGETS[name]
        X, Y, ok, _ = plot_grid(target, 0.25)
        frac = float(ok.mean())
        fracs[name] = frac
        ax.contourf(X, Y, ok.astype(float), levels=[0.5, 1.5],
                    colors=["#cfe3ff"], alpha=0.95, zorder=0)
        _draw_field(ax, target, info["r80_m"], info["arcs"],
                    info["station"], info["vehicle"])
        ax.text(target.x + 2.0, target.y - 4.0, f"{name} (H={target.H:g})",
                color="red", fontsize=8, zorder=9)
        st = info["station"]
        xytext = (16, 84) if name == "T2" else (64, 6)
        ax.annotate(f"({st['x']:.2f}, {st['y']:.2f})",
                    xy=(st["x"], st["y"]), xytext=xytext, fontsize=8,
                    zorder=9, arrowprops=dict(arrowstyle="->", lw=0.8, color="k"))
        ax.set_title(f"{name} observable region ({100 * frac:.1f}% of sampled centers)")
    fig.tight_layout()
    fig.savefig(OUT / "q1_T2_T3_regions.png", dpi=300)
    plt.close(fig)
    return fracs


def main():
    summary = {
        "model": {
            "objective": "min horizontal range r over clearance ∩ LOS ∩ elevation",
            "lower_bound": "r >= r80 = (H-2)/tan(80 deg)",
            "tie_breaking": "on the r80 optimal arc: clearance margin >= 0.10 m "
                            "beyond the 1.5 m buffer, then minimise designated "
                            "vehicle mission cycle = travel + slew + 5 s + reset; "
                            "ties: larger margin, nearer O, smaller azimuth",
            "designated_vehicle_rule": "nearest standby vehicle to target projection",
            "clearance_margin_m": CLEARANCE_MARGIN,
            "printed_station_precision_decimals": 6,
            "arc_search": {"feasible_arc_step_deg": ARC_STEP,
                           "coarse_step_deg": COARSE_STEP,
                           "fine_step_deg": FINE_STEP},
        },
        "targets": {}, "validation": {},
    }
    regions, infos = {}, {}
    for name in ("T1", "T2", "T3"):
        t = TARGETS[name]
        stats_list = [region_stats(t, res) for res in (0.2, 0.1, 0.05)]
        stats = stats_list[1]
        regions[name] = stats
        r80 = elev_min_radius(t)
        arcs, full = feasible_arcs(t, r80)
        vehicle = designated_vehicle(t)
        station = pick_station(t, r80, vehicle)
        causes = arc_endpoint_causes(t, r80, arcs, full)
        alt = nearest_O_alternative(t, r80, arcs, full, vehicle)
        areas = {f"{s['res_m']}": round(s["area_m2"], 1) for s in stats_list}
        a01 = areas["0.1"]
        half_rel = abs(areas["0.05"] - a01) / a01
        double_rel = abs(areas["0.2"] - a01) / a01
        span = 360.0 if full else sum(a1 - a0 for a0, a1 in arcs)
        summary["targets"][name] = {
            "projection": [round(t.x, 6), round(t.y, 6)],
            "H_m": t.H,
            "r80_m": round(r80, 6),
            "Q_star_per_m": round(1.0 / r80, 6),
            "r80_circle_feasible": bool(arcs),
            "full_circle": full,
            "feasible_arcs_deg": [[round(a0, 2), round(a1, 2)] for a0, a1 in arcs],
            "arc_span_deg": round(span, 2),
            "arc_endpoint_causes": causes,
            "area_m2_by_res": areas,
            "n_components_0.1m": stats["n_components"],
            "largest_component_area_m2": round(stats["largest_component_area_m2"], 1),
            "designated_vehicle": f"V{vehicle}",
            "optimal_station": {
                "x": round(station["x"], 6), "y": round(station["y"], 6),
                "azimuth_on_circle_deg": round(station["azimuth_on_circle_deg"], 3),
                "pan_deg": round(station["pan_deg"], 3),
                "tilt_deg": round(station["tilt_deg"], 3),
                "elevation_angle_deg": round(station["elevation_angle_deg"], 6),
                "range_m": round(station["range_m"], 6),
                "clearance_margin_m": round(station["clearance_margin_m"], 6),
                "slew_s": round(station["slew_s"], 3),
                "reset_s": round(station["reset_s"], 3),
                "travel_time_s": round(station["travel_time_s"], 6),
                "travel_distance_m": round(station["travel_distance_m"], 6),
                "observation_s": OBS_DURATION,
                "mission_cycle_s": round(station["mission_cycle_s"], 6),
                "dist_to_O_m": round(station["dist_to_O_m"], 6),
                "feasible": station["feasible"],
                "clearance_margin_used_m": CLEARANCE_MARGIN,
            },
            "alternative_nearest_O": alt and {
                "x": round(alt["x"], 6), "y": round(alt["y"], 6),
                "azimuth_on_circle_deg": round(alt["azimuth_on_circle_deg"], 3),
                "clearance_margin_m": round(alt["clearance_margin_m"], 6),
                "dist_to_O_m": round(alt["dist_to_O_m"], 6),
                "mission_cycle_s": alt["mission_cycle_s"] and round(alt["mission_cycle_s"], 6),
                "note": alt["note"],
            },
        }
        infos[name] = {"r80_m": r80, "arcs": arcs, "station": station,
                       "vehicle": vehicle}
        print(f"[{name}] r80={r80:.4f} m, arcs={summary['targets'][name]['feasible_arcs_deg']}, "
              f"V{vehicle} station=({station['x']:.6f},{station['y']:.6f}), "
              f"cycle={station['mission_cycle_s']:.3f} s, "
              f"margin={station['clearance_margin_m']:.3f} m")

    half_max, double_max = 0.0, 0.0
    for name in ("T1", "T2", "T3"):
        t = TARGETS[name]
        st = summary["targets"][name]["optimal_station"]
        p = (st["x"], st["y"])
        scalar_ok = validate_against_scalar(t, 0.1)
        pointwise = {
            "feasible": is_station_feasible(p, t),
            "los_clear": los_ok(p, t),
            "elevation_deg": round(elevation_angle(p, t), 6),
            "clearance_margin_m": round(clearance_margin(p), 6),
            "margin_satisfied": clearance_margin(p) >= CLEARANCE_MARGIN - 1e-9,
        }
        a = st["azimuth_on_circle_deg"]
        arcs = summary["targets"][name]["feasible_arcs_deg"]
        if summary["targets"][name]["full_circle"]:
            a_in = a + INTERIOR_STEP
        else:
            a0, a1 = arcs[0]
            cand = [a, a + 360.0, a - 360.0]
            au = next(c for c in cand if a0 - 1e-9 <= c <= a1 + 1e-9)
            a_in = au + INTERIOR_STEP if au - a0 < a1 - au else au - INTERIOR_STEP
        x_in, y_in = _circle_points(t, elev_min_radius(t), [a_in % 360.0])
        p_in = (float(x_in[0]), float(y_in[0]))
        interior = {
            "step_deg": INTERIOR_STEP,
            "station": [round(p_in[0], 6), round(p_in[1], 6)],
            "feasible": is_station_feasible(p_in, t),
            "clearance_margin_m": round(clearance_margin(p_in), 6),
            "range_m": round(math.dist(p_in, (t.x, t.y)), 6),
        }
        areas = summary["targets"][name]["area_m2_by_res"]
        half_rel = abs(areas["0.05"] - areas["0.1"]) / areas["0.1"]
        double_rel = abs(areas["0.2"] - areas["0.1"]) / areas["0.1"]
        half_max = max(half_max, half_rel)
        double_max = max(double_max, double_rel)
        summary["validation"][name] = {
            "vectorised_vs_scalar_400pts": scalar_ok,
            "pointwise": pointwise,
            "interior_arc_step": interior,
            "area_rel_change_halving_step": round(half_rel, 6),
            "area_rel_change_doubling_step": round(double_rel, 6),
        }
        print(f"[{name}] scalar互证={scalar_ok}, pointwise={pointwise['feasible']}, "
              f"interior={interior['feasible']}")

    summary["validation"]["resolution_max_rel_change"] = {
        "halving_step_0.1_to_0.05": round(half_max, 6),
        "doubling_step_0.1_to_0.2": round(double_max, 6),
    }
    frac1 = plot_t1(TARGETS["T1"], infos["T1"]["r80_m"], infos["T1"]["arcs"],
                    infos["T1"]["station"], infos["T1"]["vehicle"])
    fracs = plot_t23({k: infos[k] for k in ("T2", "T3")})
    summary["figure_rasters"] = {
        "T1": {"step_m": 0.2, "feasible_fraction": round(frac1, 6)},
        "T2": {"step_m": 0.25, "feasible_fraction": round(fracs["T2"], 6)},
        "T3": {"step_m": 0.25, "feasible_fraction": round(fracs["T3"], 6)},
    }
    (OUT / "q1_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"outputs written to {OUT}")


if __name__ == "__main__":
    main()