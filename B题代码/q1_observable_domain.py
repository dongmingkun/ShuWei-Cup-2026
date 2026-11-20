# -*- coding: utf-8 -*-
"""Q1: 可观测域网格计算(T1 主线, T2/T3 对比) + 最优站求解 + 独立验证.

运行: python -X utf8 -B q1_observable_domain.py
输出: q1_output/ 下 PNG 图件与 q1_summary.json, 控制台打印关键结论.

口径与 geometry.py 一致: 净距 >= 1.5 m 贴边允许; 俯仰 [-10°,80°] 闭区间;
视线切触即遮挡. 网格法为 numpy 向量化实现, 与 geometry.py 标量谓词
随机抽样互证; 最优弧按 0.01° 角步长(在 r≈7.58 m 处约 1.3 mm 弧长)扫描.
"""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from geometry import (BUILDINGS, DISPATCH, FIELD_SIZE, TARGETS,
                      VEHICLE_RADIUS, clearance_ok, elev_min_radius,
                      elevation_angle, is_station_feasible)
from simulator import rotate_duration

OUT = Path(__file__).with_name("q1_output")
OUT.mkdir(exist_ok=True)


# ---------------------------------------------------------------- 向量化网格
def _los_blocked_vec(X, Y, Tx, Ty, H, b):
    """网格点上 (X,Y,2)->(Tx,Ty,H) 是否与建筑棱柱相交(闭集, 切触算遮挡)."""
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


def region_mask(target, res):
    """可观测域示性函数网格: 三约束(净距/俯仰/无遮挡)同时成立."""
    n = int(round(FIELD_SIZE / res)) + 1
    xs = np.linspace(0.0, FIELD_SIZE, n)
    X, Y = np.meshgrid(xs, xs)
    m = ((X >= VEHICLE_RADIUS - 1e-9) & (X <= FIELD_SIZE - VEHICLE_RADIUS + 1e-9)
         & (Y >= VEHICLE_RADIUS - 1e-9) & (Y <= FIELD_SIZE - VEHICLE_RADIUS + 1e-9))
    r = np.hypot(X - target.x, Y - target.y)
    m &= r >= elev_min_radius(target) - 1e-9          # 俯仰 <= 80°(闭区间)
    for b in BUILDINGS:
        ddx = np.maximum(np.maximum(b.x0 - X, X - b.x1), 0.0)
        ddy = np.maximum(np.maximum(b.y0 - Y, Y - b.y1), 0.0)
        m &= np.hypot(ddx, ddy) >= VEHICLE_RADIUS - 1e-9   # 净距贴边允许
        m &= ~_los_blocked_vec(X, Y, target.x, target.y, target.H, b)
    return xs, m


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
    """随机抽样: 向量化网格判定 vs geometry.py 标量谓词必须一致."""
    rng = np.random.default_rng(seed)
    xs, m = region_mask(target, res)
    ok = True
    for _ in range(n_samples):
        p = (float(rng.uniform(0, FIELD_SIZE)), float(rng.uniform(0, FIELD_SIZE)))
        i = int(round(p[0] / res)); j = int(round(p[1] / res))
        i = min(max(i, 0), len(xs) - 1); j = min(max(j, 0), len(xs) - 1)
        vec = bool(m[j, i])
        p_snap = (float(xs[i]), float(xs[j]))
        sca = is_station_feasible(p_snap, target)
        if vec != sca:
            ok = False
            print(f"  MISMATCH at {p_snap}: vec={vec} scalar={sca}")
    return ok


# ---------------------------------------------------------------- 最优站
def _circle_point(target, r80, a_deg):
    a = math.radians(a_deg)
    return (target.x + r80 * math.cos(a), target.y + r80 * math.sin(a))


def feasible_arcs(target, r80, step_deg=0.01):
    """r80 圆上的可行角区间(闭集, 与三约束一致); 返回 [(a0,a1), ...]."""
    n = int(round(360.0 / step_deg))
    feas = np.array([is_station_feasible(_circle_point(target, r80, i * step_deg),
                                         target) for i in range(n)])
    idx = np.flatnonzero(feas)
    if idx.size == 0:
        return [], feas
    runs, start, prev = [], idx[0], idx[0]
    for k in idx[1:]:
        if k != prev + 1:
            runs.append((start, prev))
            start = k
        prev = k
    runs.append((start, prev))
    if len(runs) > 1 and runs[0][0] == 0 and runs[-1][1] == n - 1:  # 跨 0° 环绕
        s_last, e_last = runs.pop()
        s_first, e_first = runs.pop(0)
        runs.append((s_last - n, e_first))
        runs.sort()
    return [(float(s * step_deg), float(e * step_deg)) for s, e in runs], feas


def pick_station(target, r80, arcs, step_deg=0.01):
    """破平规则: ①转位+复位耗时最小 ②距调度中心 O 最近."""
    best = None
    for a0, a1 in arcs:
        k = int(round((a1 - a0) / step_deg))
        for i in range(k + 1):
            a = a0 + i * step_deg
            p = _circle_point(target, r80, a % 360.0)
            pan = math.degrees(math.atan2(target.y - p[1], target.x - p[0])) % 360.0
            slew = rotate_duration(0.0, 0.0, pan, 80.0)
            reset = rotate_duration(pan, 80.0, 0.0, 0.0)
            key = (slew + reset, math.dist(p, DISPATCH))
            if best is None or key < best["key"]:
                best = {"key": key, "angle_deg": a % 360.0, "station": p,
                        "pan_deg": pan, "tilt_deg": 80.0,
                        "slew_s": slew, "reset_s": reset,
                        "dist_to_O_m": math.dist(p, DISPATCH)}
    best.pop("key")
    return best


def arc_endpoint_causes(target, r80, arcs, step_deg=0.01):
    """每个可行弧端点外侧 0.01° 处的不可行原因(净距/俯仰/遮挡)."""
    causes = []
    for a0, a1 in arcs:
        for a in (a0 - step_deg, a1 + step_deg):
            p = _circle_point(target, r80, a % 360.0)
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


# ---------------------------------------------------------------- 绘图
def plot_region(ax, target, stats, station, r80, arcs, title):
    xs = stats["xs"]
    ax.imshow(stats["mask"].T, origin="lower", extent=(0, 100, 0, 100),
              cmap=matplotlib.colors.ListedColormap(["#ffffff", "#c9e7cf"]),
              vmin=0, vmax=1, interpolation="nearest")
    for k, b in enumerate(BUILDINGS, start=1):
        ax.add_patch(Rectangle((b.x0, b.y0), b.x1 - b.x0, b.y1 - b.y0,
                               facecolor="#9e9e9e", edgecolor="k", zorder=3))
        ax.text((b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2, f"B{k}", ha="center",
                va="center", fontsize=8, zorder=4)
    ax.plot(target.x, target.y, "^", color="crimson", ms=9, zorder=5)
    ax.annotate(f"{target.name} ({target.x:.2f},{target.y:.2f},{target.H:.0f})",
                (target.x, target.y), textcoords="offset points", xytext=(6, 6),
                fontsize=8, zorder=5)
    th = np.linspace(0, 2 * np.pi, 721)
    ax.plot(target.x + r80 * np.cos(th), target.y + r80 * np.sin(th), "--",
            color="gray", lw=0.8, zorder=4)
    for a0, a1 in arcs:
        aa = np.radians(np.linspace(a0, a1, max(2, int((a1 - a0) / 0.002))))
        ax.plot(target.x + r80 * np.cos(aa), target.y + r80 * np.sin(aa),
                color="crimson", lw=2.5, zorder=6)
    if station:
        ax.plot(*station["station"], "*", color="royalblue", ms=14, zorder=7)
        ax.annotate(f"optimal r={r80:.3f} m", station["station"],
                    textcoords="offset points", xytext=(8, -12), fontsize=8,
                    color="royalblue", zorder=7)
    ax.plot(*DISPATCH, "x", color="k", ms=6, zorder=5)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_aspect("equal"); ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
    ax.set_title(title, fontsize=10)


def main():
    summary = {"targets": {}, "validation": {}}
    regions = {}
    for name in ("T1", "T2", "T3"):
        t = TARGETS[name]
        stats_list = [region_stats(t, res) for res in (0.2, 0.1, 0.05)]
        stats = stats_list[1]                       # 0.1 m 主线
        regions[name] = stats
        r80 = elev_min_radius(t)
        arcs, feas = feasible_arcs(t, r80)
        station = pick_station(t, r80, arcs) if arcs else None
        causes = arc_endpoint_causes(t, r80, arcs)
        summary["targets"][name] = {
            "projection": [round(t.x, 4), round(t.y, 4)],
            "H_m": t.H, "r80_m": round(r80, 4),
            "r80_circle_feasible": bool(feas.any()),
            "feasible_arcs_deg": [[round(a0, 2), round(a1, 2)] for a0, a1 in arcs],
            "arc_endpoint_causes": causes,
            "optimal_station": station and {
                "x": round(station["station"][0], 4),
                "y": round(station["station"][1], 4),
                "azimuth_on_circle_deg": round(station["angle_deg"], 2),
                "pan_deg": round(station["pan_deg"], 3),
                "tilt_deg": station["tilt_deg"],
                "elevation_angle_deg": round(elevation_angle(station["station"], t), 3),
                "slew_s": round(station["slew_s"], 3),
                "reset_s": round(station["reset_s"], 3),
                "dist_to_O_m": round(station["dist_to_O_m"], 3),
            },
            "area_m2_by_res": {f"{s['res_m']}": round(s["area_m2"], 1) for s in stats_list},
            "n_components_0.1m": stats["n_components"],
            "largest_component_area_m2": round(stats["largest_component_area_m2"], 1),
        }
        print(f"[{name}] r80={r80:.4f} m, 可行弧={summary['targets'][name]['feasible_arcs_deg']}")
        if station:
            print(f"      最优站 {tuple(round(v, 3) for v in station['station'])}, "
                  f"pan={station['pan_deg']:.2f}°, slew+reset={station['slew_s'] + station['reset_s']:.2f} s, "
                  f"|O|={station['dist_to_O_m']:.2f} m")

    for name in ("T1", "T2", "T3"):
        summary["validation"][name] = validate_against_scalar(TARGETS[name], 0.1)
        print(f"[{name}] 向量化 vs 标量谓词 400 点互证: "
              f"{'一致' if summary['validation'][name] else '不一致!'}")

    fig, ax = plt.subplots(figsize=(7, 6.8))
    t1 = TARGETS["T1"]
    plot_region(ax, t1, regions["T1"],
                summary["targets"]["T1"]["optimal_station"] and
                {"station": (summary["targets"]["T1"]["optimal_station"]["x"],
                             summary["targets"]["T1"]["optimal_station"]["y"])},
                elev_min_radius(t1),
                summary["targets"]["T1"]["feasible_arcs_deg"],
                f"Observable region of T1 (grid 0.1 m, area "
                f"{regions['T1']['area_m2']:.0f} m$^2$)")
    fig.tight_layout()
    fig.savefig(OUT / "q1_T1_region.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.6))
    for ax, name in zip(axes, ("T2", "T3")):
        t = TARGETS[name]
        st = summary["targets"][name]["optimal_station"]
        plot_region(ax, t, regions[name],
                    st and {"station": (st["x"], st["y"])},
                    elev_min_radius(t),
                    summary["targets"][name]["feasible_arcs_deg"],
                    f"Observable region of {name} (grid 0.1 m, area "
                    f"{regions[name]['area_m2']:.0f} m$^2$)")
    fig.tight_layout()
    fig.savefig(OUT / "q1_T2_T3_regions.png", dpi=300)
    plt.close(fig)

    (OUT / "q1_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n输出已写入 {OUT}")


if __name__ == "__main__":
    main()
