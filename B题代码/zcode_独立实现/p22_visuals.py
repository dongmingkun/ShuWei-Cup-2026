# -*- coding: utf-8 -*-
"""Minimal visuals for the sigma<2m Q2(2) script (draw_gantt / draw_trajectory).

Implements the two entry points used by problem2_2_2m.py against the zcode
swb_common engine (field layout, buildings, timeline events)."""
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import swb_common as S

COL = {"depart": "#1f4fd8", "arrive": "#1f4fd8", "ready": "#7a9bd9",
       "reset_done": "#c9d6ec", "available": "#d8d8d8"}


def draw_gantt(sim, title, path):
    vehs = [v for v in sim["timeline"] if sim["timeline"][v]]
    fig, ax = plt.subplots(figsize=(9.5, 2.6 + 0.45 * len(vehs)), dpi=300)
    ypos = {v: i for i, v in enumerate(vehs)}
    for v in vehs:
        evs = sim["timeline"][v]
        for k in range(len(evs) - 1):
            e, n = evs[k], evs[k + 1]
            label = e.get("alert", "")
            color = {"depart": "#1f4fd8", "arrive": "#1f4fd8",
                     "ready": "#e0a33a", "reset_done": "#c9d6ec"}.get(
                e["kind"], "#c9c9c9")
            if e["kind"] in ("depart", "ready") and n["t"] > e["t"]:
                ax.barh(ypos[v], n["t"] - e["t"], left=e["t"], height=0.55,
                        color=color, edgecolor="k", lw=0.4)
                if e["kind"] == "ready":
                    ax.text(e["t"] + (n["t"] - e["t"]) / 2, ypos[v], label,
                            ha="center", va="center", fontsize=7)
            if e["kind"] == "depart":
                ax.text(e["t"] + 0.05, ypos[v] + 0.36, label, fontsize=7,
                        va="bottom")
    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels(vehs)
    ax.set_xlabel("time (s)")
    ax.set_xlim(0, max(22.0, max((e["t"] for v in vehs for e in
                 sim["timeline"][v]), default=0) * 1.08))
    ax.set_title(title, fontsize=10)
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def draw_trajectory(st, asg, info, max_start, path, title):
    fig, ax = plt.subplots(figsize=(8.4, 8.4), dpi=300)
    for (x0, x1, y0, y1, h), name in zip(S.BUILDINGS, S.BUILD_NAMES):
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                               facecolor="#5a5a5a", edgecolor="k", zorder=3))
        ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"{name}\n{h:g} m", color="w",
                ha="center", va="center", fontsize=8, zorder=4)
    for i, (vx, vy) in enumerate(S.V_INIT):
        ax.plot(vx, vy, "s", color="#1f4fd8", ms=6, mec="k", mew=0.4, zorder=6)
        ax.text(vx + 1.4, vy + 1.4, S.V_NAMES[i], fontsize=7.5,
                color="#1f4fd8", zorder=7)
    colors = {"T1": "#c0392b", "T2": "#1e8449", "T3": "#8e44ad"}
    for a, (tx, ty, H) in S.TARGETS.items():
        ax.plot(tx, ty, "^", color=colors.get(a, "red"), ms=11, zorder=7)
        ax.text(tx + 1.2, ty + 1.6, a, fontsize=8, color=colors.get(a, "red"),
                zorder=8, bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                    ec="none", alpha=0.8))
    for (v, a), p in st.items():
        ax.plot(*p, "o", color=colors.get(a, "gray"), ms=7, mec="k", mew=0.5,
                zorder=8)
        tx, ty, _ = S.TARGETS[a]
        ax.plot([p[0], tx], [p[1], ty], color=colors.get(a, "gray"), ls="--",
                lw=0.8, alpha=0.7, zorder=5)
    for v, al in asg.items():
        if not al:
            continue
        p0 = S.V_INIT[S.V_NAMES.index(v)]
        for a in al:
            p = st[(v, a)]
            ax.annotate("", xy=p, xytext=p0,
                        arrowprops=dict(arrowstyle="->", ls="--", lw=0.9,
                                        color=colors.get(a, "gray")), zorder=6)
            p0 = p
    if info:
        lines = [f"max start = {max_start:.2f} s"]
        for a in ("T1", "T2", "T3"):
            q = info.get(a)
            if q:
                lines.append(f"{a}: gamma={q['gamma']:.0f} deg, "
                             f"sigma={q['sigma']:.3f} m")
        ax.text(0.02, 0.02, "\n".join(lines), transform=ax.transAxes,
                fontsize=8, va="bottom",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.6",
                          alpha=0.92), zorder=12)
    ax.set_xlim(-2, 102)
    ax.set_ylim(-2, 102)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m) east")
    ax.set_ylabel("y (m) north")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)
