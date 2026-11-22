# -*- coding: utf-8 -*-
"""Q2(1) accuracy audit — pair comparison vs the Fisher bound.

Reads results/problem2_1_results_v2.json (the margin-compliant optimum from
problem2_1_v2.py) and emits the correctly-NAMED audit artefacts:
    results/fig_q2_1_accuracy_audit.png   (bar chart: RMSE + travel/mission, bound line)
    results/q2_1_accuracy_audit.json      (audit data)

NOTE on naming: an earlier draft of this artefact was mislabelled "q2_2_*";
it is Q2(1) content (two-station accuracy), NOT Q2(2).  Q2(2)'s scheduler is
untouched by this audit.
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

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
import swb_common as S  # noqa: E402

OUT = BASE / "results"
SRC = OUT / "problem2_1_results_v2.json"

d = json.load(open(SRC, encoding="utf-8"))
bound = d["theoretical_bound"]["sigma_lower_bound_m"]
rows = d["pair_table"]
sel = d["robust_optimal"]

pairs = [r["pair"] for r in rows]
sigma = [r["sigma_rms_m"] for r in rows]
travel = [r["travel_total_s"] for r in rows]
mission = [r["mission_time_s"] for r in rows]
best_pair = sel["vehicles"][0] + "+" + sel["vehicles"][1]

# ---------------------------------------------------------------- figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.0), dpi=300)
x = np.arange(len(pairs))

b1 = ax1.bar(x, sigma, width=0.55, color="#1f4fd8",
             label="margin-compliant RMSE")
ax1.axhline(bound, color="#c62828", ls="--", lw=1.4,
            label=f"Fisher bound = {bound:.5f} m")
for xi, yi in zip(x, sigma):
    ax1.text(xi, yi + 0.0015, f"{yi:.4f}", ha="center", fontsize=8)
ax1.set_xticks(x)
ax1.set_xticklabels(pairs, fontsize=8.5)
ax1.set_ylabel("RMSE (m)")
ax1.set_title("(a) Q2(1) accuracy audit: every margin-compliant pair sits "
              f"+{(sigma[0]/bound-1)*100:.1f}% over the bound", fontsize=10)
ax1.legend(fontsize=8.5)
ax1.grid(axis="y", alpha=0.3)

w = 0.38
b2 = ax2.bar(x - w/2, travel, width=w, color="#ef8a15", label="total travel (s)")
b3 = ax2.bar(x + w/2, mission, width=w, color="#1f4fd8",
             label="mission completion (s)")
best_i = pairs.index(best_pair)
for xi, yi in zip(x, mission):
    ax2.text(xi + w/2, yi + 0.6, f"{yi:.1f}", ha="center", fontsize=8)
ax2.set_xticks(x)
ax2.set_xticklabels(pairs, fontsize=8.5)
ax2.set_ylabel("time (s)")
ax2.set_title(f"(b) time audit: mission-time tie-break selects {best_pair}",
              fontsize=10)
ax2.legend(fontsize=8.5)
ax2.grid(axis="y", alpha=0.3)

fig.suptitle("Q2(1) accuracy audit — margin-compliant pairs vs the Fisher bound "
             "(Q2(1), not Q2(2))", fontsize=11, y=0.99)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(OUT / "fig_q2_1_accuracy_audit.png", dpi=300)
print("figure saved:", OUT / "fig_q2_1_accuracy_audit.png")

# ---------------------------------------------------------------- json
audit = {
    "artefact": "q2_1_accuracy_audit",
    "belongs_to": "Q2(1) two-station accuracy (NOT Q2(2); the scheduler is untouched)",
    "source": "problem2_1_results_v2.json",
    "bound": {
        "fisher_lower_bound_m": bound,
        "attained_when": "margins relaxed to zero, gamma=90 deg, d1=d2=r80",
    },
    "margins": d["margins"],
    "margin_cost_pct": round((sigma[0] / bound - 1) * 100, 3),
    "pairs": [
        {"pair": p, "sigma_rms_m": s, "over_bound_pct": round((s/bound-1)*100, 2),
         "total_travel_s": t, "mission_s": m}
        for p, s, t, m in zip(pairs, sigma, travel, mission)
    ],
    "selected": {
        "pair": best_pair,
        "stations": sel["stations"],
        "sigma_rms_m": sel["sigma_rms_m"],
        "travel_times_s": sel["travel_times_s"],
        "mission_s": sel["mission_time_s"],
        "margins": {"clearance_m": sel["clearance_m"],
                    "elevation_deg": sel["elevation_deg"]},
    },
    "verdict": ("all margin-compliant pairs (excluding V2) sit at the same "
                "RMSE = bound + margin cost; the mission-time tie-break selects "
                f"{best_pair}"),
}
with open(OUT / "q2_1_accuracy_audit.json", "w", encoding="utf-8") as fp:
    json.dump(audit, fp, ensure_ascii=False, indent=2)
print("results saved:", OUT / "q2_1_accuracy_audit.json")
