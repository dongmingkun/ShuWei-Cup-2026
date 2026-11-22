# -*- coding: utf-8 -*-
"""Q2(2) 精度审计：保留时间主方案并复核 30° 保底变体；运行 python -X utf8 -B q2_2_accuracy_audit.py。"""
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent / "q2_output" / "accuracy_audit"
OUT.mkdir(parents=True, exist_ok=True)
SIGMA_B = math.radians(1.0)
Q21_REFERENCE_RMSE_M = 0.14362
ALERTS = ("T1", "T2", "T3")
TARGET_SPECS = {
    "T1": (30.0, 40.0, 45.0),
    "T2": (150.0, 45.0, 35.0),
    "T3": (285.0, 38.0, 55.0),
}
TARGETS = {
    name: (50.0 + R * math.cos(math.radians(theta)),
           50.0 + R * math.sin(math.radians(theta)), H)
    for name, (theta, R, H) in TARGET_SPECS.items()
}
BUILDINGS = (
    (22.0, 34.0, 16.0, 26.0, 20.0),
    (63.0, 79.0, 59.0, 71.0, 28.0),
    (10.0, 24.0, 65.0, 76.0, 12.0),
    (72.0, 86.0, 22.0, 34.0, 24.0),
    (41.0, 51.0, 38.0, 46.0, 32.0),
)
PLANS = {
    "time_first": {
        "label": "Time-first",
        "azimuth_guard_deg": 1.0,
        "max_start_s": 11.046,
        "makespan_s": 21.705,
        "alert_start_s": {"T1": 2.943, "T2": 10.730, "T3": 11.046},
        "reported_gamma_deg": {"T1": 17.65, "T2": 50.23, "T3": 59.21},
        "stations": {
            "T1": (("V3", (18.0, 90.0)), ("V5", (51.0, 93.0))),
            "T2": (("V2", (70.0, 3.5)), ("V4", (87.5, 73.5))),
            "T3": (("V3", (18.0, 90.0)), ("V1", (15.0, 15.0))),
        },
        "source": "q2_2_summary.json, azimuth_guard_deg=1.0",
    },
    "accuracy_guarded": {
        "label": "Accuracy-guarded",
        "azimuth_guard_deg": 30.0,
        "max_start_s": 11.649,
        "makespan_s": 18.791,
        "alert_start_s": {"T1": 3.545, "T2": 7.816, "T3": 11.649},
        "reported_gamma_deg": {"T1": 108.92, "T2": 31.80, "T3": 59.21},
        "stations": {
            "T1": (("V3", (18.0, 90.0)), ("V2", (87.0, 9.0))),
            "T2": (("V4", (87.5, 73.5)), ("V5", (45.5, 94.5))),
            "T3": (("V3", (18.0, 90.0)), ("V1", (15.0, 15.0))),
        },
        "source": "q2_2_scheduler.py frozen result, azimuth_guard_deg=30.0",
    },
}


def dist_point_to_rect(p, building):
    x, y = p
    x0, x1, y0, y1, _ = building
    return math.hypot(max(x0 - x, 0.0, x - x1), max(y0 - y, 0.0, y - y1))


def clearance_ok(p):
    x, y = p
    if not (1.5 - 1e-9 <= x <= 98.5 + 1e-9 and 1.5 - 1e-9 <= y <= 98.5 + 1e-9):
        return False
    return all(dist_point_to_rect(p, b) >= 1.5 - 1e-9 for b in BUILDINGS)


def segment_hits_box(p0, p1, building):
    x0, x1, y0, y1, h = building
    bounds = ((x0, x1), (y0, y1), (0.0, h))
    direction = tuple(p1[i] - p0[i] for i in range(3))
    tmin, tmax = 0.0, 1.0
    for axis, (lo, hi) in enumerate(bounds):
        d = direction[axis]
        if abs(d) < 1e-15:
            if p0[axis] < lo - 1e-9 or p0[axis] > hi + 1e-9:
                return False
            continue
        t1, t2 = (lo - p0[axis]) / d, (hi - p0[axis]) / d
        tmin = max(tmin, min(t1, t2))
        tmax = min(tmax, max(t1, t2))
        if tmin > tmax + 1e-9:
            return False
    return True


def station_feasible(p, alert):
    tx, ty, H = TARGETS[alert]
    if not clearance_ok(p):
        return False
    r = math.hypot(tx - p[0], ty - p[1])
    elevation = math.degrees(math.atan2(H - 2.0, r))
    if not (-10.0 - 1e-9 <= elevation <= 80.0 + 1e-9):
        return False
    p0, p1 = (p[0], p[1], 2.0), (tx, ty, H)
    return not any(segment_hits_box(p0, p1, b) for b in BUILDINGS)


def bearing_data(station, target):
    dx, dy = target[0] - station[0], target[1] - station[1]
    distance = math.hypot(dx, dy)
    phi = math.atan2(dy, dx)
    return distance, phi, (-math.sin(phi), math.cos(phi)), 1.0 / distance**2


def alert_metrics(alert, station_records):
    target = TARGETS[alert]
    (v1, s1), (v2, s2) = station_records
    d1, phi1, u1, w1 = bearing_data(s1, target)
    d2, phi2, u2, w2 = bearing_data(s2, target)
    angle_difference = math.degrees(phi2 - phi1) % 360.0
    gamma = angle_difference if angle_difference <= 180.0 else 360.0 - angle_difference
    rmse = SIGMA_B * math.sqrt(d1 * d1 + d2 * d2) / math.sin(math.radians(gamma))
    a11 = w1 * u1[0] ** 2 + w2 * u2[0] ** 2
    a22 = w1 * u1[1] ** 2 + w2 * u2[1] ** 2
    a12 = w1 * u1[0] * u1[1] + w2 * u2[0] * u2[1]
    fisher_rmse = SIGMA_B * math.sqrt((a11 + a22) / (a11 * a22 - a12 * a12))
    assert math.isclose(rmse, fisher_rmse, rel_tol=1e-12, abs_tol=1e-12)
    assert station_feasible(s1, alert) and station_feasible(s2, alert)
    return {
        "alert": alert,
        "vehicles": [v1, v2],
        "stations_xy_m": [list(s1), list(s2)],
        "target_xy_m": [target[0], target[1]],
        "station_distances_m": [d1, d2],
        "baseline_m": math.dist(s1, s2),
        "crossing_angle_deg": gamma,
        "rmse_m": rmse,
        "fisher_rmse_crosscheck_m": fisher_rmse,
        "stations_feasible": True,
    }


def audit_plan(plan_id, plan):
    alerts = {}
    for alert in ALERTS:
        metrics = alert_metrics(alert, plan["stations"][alert])
        assert abs(metrics["crossing_angle_deg"] - plan["reported_gamma_deg"][alert]) < 0.01
        assert metrics["crossing_angle_deg"] >= plan["azimuth_guard_deg"] - 0.01
        alerts[alert] = metrics
    assert math.isclose(max(plan["alert_start_s"].values()), plan["max_start_s"], abs_tol=1e-9)
    return {
        "label": plan["label"],
        "source": plan["source"],
        "azimuth_guard_deg": plan["azimuth_guard_deg"],
        "max_start_s": plan["max_start_s"],
        "makespan_s": plan["makespan_s"],
        "alert_start_s": plan["alert_start_s"],
        "alerts": alerts,
        "worst_rmse_m": max(metrics["rmse_m"] for metrics in alerts.values()),
        "min_crossing_angle_deg": min(metrics["crossing_angle_deg"] for metrics in alerts.values()),
    }


def write_csv(results):
    path = OUT / "q2_2_accuracy_audit.csv"
    with path.open("w", encoding="utf-8-sig", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["plan", "guard_deg", "max_start_s", "alert", "vehicles",
                         "d1_m", "d2_m", "gamma_deg", "rmse_m"])
        for result in results.values():
            for alert in ALERTS:
                m = result["alerts"][alert]
                writer.writerow([result["label"], result["azimuth_guard_deg"],
                                 f"{result['max_start_s']:.3f}", alert, "+".join(m["vehicles"]),
                                 f"{m['station_distances_m'][0]:.4f}",
                                 f"{m['station_distances_m'][1]:.4f}",
                                 f"{m['crossing_angle_deg']:.3f}", f"{m['rmse_m']:.4f}"])


def write_figure(results):
    original = results["time_first"]
    guarded = results["accuracy_guarded"]
    colors = ("#4C78A8", "#F58518")
    x = np.arange(len(ALERTS))
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3))
    for axis, key, ylabel in ((axes[0], "rmse_m", "Bearing-only RMSE (m)"),
                              (axes[1], "crossing_angle_deg", "Crossing angle (deg)")):
        for offset, result, color, label in (
                (-0.18, original, colors[0], "Time-first, max start 11.046 s"),
                (0.18, guarded, colors[1], "Accuracy-guarded, max start 11.649 s")):
            values = [result["alerts"][a][key] for a in ALERTS]
            bars = axis.bar(x + offset, values, width=0.36, color=color, label=label)
            axis.bar_label(bars, fmt="%.2f", padding=2, fontsize=8)
        axis.set_xticks(x, ALERTS)
        axis.set_ylabel(ylabel)
        axis.grid(axis="y", alpha=0.25)
        axis.set_axisbelow(True)
    axes[0].set_ylim(0, 5.4)
    axes[1].axhline(30.0, color="black", ls="--", lw=0.9, label="30 deg guard")
    axes[1].set_ylim(0, 125)
    axes[0].legend(fontsize=8, loc="upper right")
    axes[1].legend(fontsize=8, loc="upper right")
    fig.suptitle("Q2(2) accuracy audit: retain the time-first plan and report a guarded variant")
    fig.tight_layout()
    fig.savefig(OUT / "q2_2_accuracy_tradeoff.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_latex(results, comparison):
    a = results["time_first"]["alerts"]
    b = results["accuracy_guarded"]["alerts"]
    distances = [d for result in results.values() for m in result["alerts"].values()
                 for d in m["station_distances_m"]]
    text = rf"""% Generated by q2_2_accuracy_audit.py; the 11.046 s schedule remains the main Q2(2)/Q4 baseline.
\textbf{{Accuracy audit and guarded variant.}} Using the same $1^\circ$ independent
bearing-error model as Q2(1), the time-first schedule has RMSE values of
{a['T1']['rmse_m']:.3f}, {a['T2']['rmse_m']:.3f}, and {a['T3']['rmse_m']:.3f} m for T1--T3.
The largest error occurs at T1 because its crossing angle is only
{a['T1']['crossing_angle_deg']:.2f}$^\circ$. We therefore retain the
$11.046$ s schedule as the time-objective solution, but additionally report a
$\gamma\ge30^\circ$ guarded variant. Its max start increases by only
{comparison['max_start_delta_s']:.3f} s to $11.649$ s, while T1's crossing angle rises to
{b['T1']['crossing_angle_deg']:.2f}$^\circ$ and its RMSE falls to
{b['T1']['rmse_m']:.3f} m; the worst-alert RMSE falls from
{results['time_first']['worst_rmse_m']:.3f} to {results['accuracy_guarded']['worst_rmse_m']:.3f} m.
T2 changes from {a['T2']['rmse_m']:.3f} to {b['T2']['rmse_m']:.3f} m, so this is a
speed--accuracy trade-off, not a uniform improvement. The remaining metre-level errors
mainly reflect the {min(distances):.1f}--{max(distances):.1f} m station--target distances,
unlike the close-range $90^\circ$ geometry specifically optimized in Q2(1)
(${Q21_REFERENCE_RMSE_M:.5f}$ m).

\begin{{table}}[!htb]\centering\footnotesize\setlength{{\tabcolsep}}{{3.5pt}}
\caption{{Q2(2) accuracy audit. Each alert cell reports crossing angle $\gamma$ (deg) /
RMSE (m); the time-first schedule remains the main time-objective solution.}}
\label{{tab:q22accuracy}}
\begin{{tabular}}{{lccccc}}
\toprule
Schedule & max start (s) & T1 $\gamma$/RMSE & T2 $\gamma$/RMSE & T3 $\gamma$/RMSE & worst RMSE (m) \\
\midrule
Time-first ($\gamma\ge1^\circ$) & 11.046 &
{a['T1']['crossing_angle_deg']:.2f}/{a['T1']['rmse_m']:.3f} &
{a['T2']['crossing_angle_deg']:.2f}/{a['T2']['rmse_m']:.3f} &
{a['T3']['crossing_angle_deg']:.2f}/{a['T3']['rmse_m']:.3f} & {results['time_first']['worst_rmse_m']:.3f} \\
Accuracy-guarded ($\gamma\ge30^\circ$) & 11.649 &
{b['T1']['crossing_angle_deg']:.2f}/{b['T1']['rmse_m']:.3f} &
{b['T2']['crossing_angle_deg']:.2f}/{b['T2']['rmse_m']:.3f} &
{b['T3']['crossing_angle_deg']:.2f}/{b['T3']['rmse_m']:.3f} & {results['accuracy_guarded']['worst_rmse_m']:.3f} \\
\bottomrule
\end{{tabular}}
\end{{table}}
"""
    (OUT / "q2_2_accuracy_for_paper.tex").write_text(text, encoding="utf-8")


def main():
    results = {plan_id: audit_plan(plan_id, plan) for plan_id, plan in PLANS.items()}
    original, guarded = results["time_first"], results["accuracy_guarded"]
    comparison = {
        "max_start_delta_s": guarded["max_start_s"] - original["max_start_s"],
        "max_start_delta_pct": (guarded["max_start_s"] / original["max_start_s"] - 1.0) * 100.0,
        "worst_rmse_before_m": original["worst_rmse_m"],
        "worst_rmse_after_m": guarded["worst_rmse_m"],
        "worst_rmse_reduction_pct": (1.0 - guarded["worst_rmse_m"] / original["worst_rmse_m"]) * 100.0,
        "t1_rmse_reduction_pct": (1.0 - guarded["alerts"]["T1"]["rmse_m"]
                                  / original["alerts"]["T1"]["rmse_m"]) * 100.0,
        "t2_rmse_delta_m": guarded["alerts"]["T2"]["rmse_m"] - original["alerts"]["T2"]["rmse_m"],
    }
    summary = {
        "scope": "Accuracy and station-feasibility audit of two frozen Q2(2) schedules; no re-optimization.",
        "main_decision": "Retain the 11.046 s time-first schedule as the Q2(2) result and Q4 baseline; report the 30 deg variant as a speed-accuracy trade-off.",
        "sigma_b_deg": 1.0,
        "rmse_formula": "sigma_b * sqrt(d1^2 + d2^2) / sin(gamma)",
        "q2_1_reference": {
            "rmse_m": Q21_REFERENCE_RMSE_M,
            "note": "Different objective and geometry: close-range T2 stations with crossing angle 90 deg; not a Q2(2) feasibility threshold.",
        },
        "results": results,
        "comparison": comparison,
    }
    (OUT / "q2_2_accuracy_audit.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    write_csv(results)
    write_figure(results)
    write_latex(results, comparison)

    print("Q2(2) 精度审计（原方案保留，30° 方案仅作保底变体）")
    for result in results.values():
        print(f"\n{result['label']}: max start={result['max_start_s']:.3f} s, "
              f"worst RMSE={result['worst_rmse_m']:.3f} m")
        for alert in ALERTS:
            m = result["alerts"][alert]
            print(f"  {alert}: {'+'.join(m['vehicles'])}, "
                  f"d=({m['station_distances_m'][0]:.2f},{m['station_distances_m'][1]:.2f}) m, "
                  f"gamma={m['crossing_angle_deg']:.2f} deg, RMSE={m['rmse_m']:.3f} m")
    print(f"\n时间代价: +{comparison['max_start_delta_s']:.3f} s "
          f"({comparison['max_start_delta_pct']:.2f}%)")
    print(f"T1 RMSE 降低: {comparison['t1_rmse_reduction_pct']:.2f}%")
    print(f"最坏 RMSE 降低: {comparison['worst_rmse_reduction_pct']:.2f}%")
    print(f"T2 RMSE 变化: {comparison['t2_rmse_delta_m']:+.3f} m（并非三项全部改善）")
    print(f"输出目录: {OUT}")


if __name__ == "__main__":
    main()
