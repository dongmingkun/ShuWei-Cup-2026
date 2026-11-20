# -*- coding: utf-8 -*-
"""论文补充图: 全文技术路线图 + 场区几何与三约束原理示意图.
AI 辅助实现(QoderCN/GLM-5.3-Flash, 2026-11-21): 图面布局与配色由脚本生成,
内容口径与《B题-整卷解题框架》一致, 数值均来自已验证的求解输出.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle, Polygon, Wedge

import sys
sys.path.insert(0, ".")
from geometry import BUILDINGS, VEHICLE_INIT, DISPATCH, TARGETS

OUT = "../论文/figures"

BLUE = "#2166AC"
ORANGE = "#D95F02"
GREEN = "#1B7837"
GRAY = "#4D4D4D"
LBLUE = "#D6E4F0"
LORANGE = "#FDE0C8"
LGREEN = "#DCEBD9"
LGRAY = "#EFEFEF"


def box(ax, x, y, w, h, text, fc, ec, fs=8.5, weight="normal"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.06", fc=fc, ec=ec, lw=1.1))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            color="#1a1a1a", weight=weight, linespacing=1.3)


def arrow(ax, x0, y0, x1, y1, color=GRAY, style="-", lw=1.2, rad=0.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=11, color=color, lw=lw,
                                 linestyle=style,
                                 connectionstyle=f"arc3,rad={rad}"))


def fig_framework():
    fig, ax = plt.subplots(figsize=(10.5, 7.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    box(ax, 50, 94, 96, 8.2,
        "Common inputs & preprocessing\n"
        "target conversion $(\\theta,R,H)\\to(x,y,H)$;  field, buildings $B_1$-$B_5$, vehicles $V_1$-$V_5$;\n"
        "timing rules: slew 30$^\\circ$/s, joint observation 5 s, evaluation window $[0,60]$ s",
        LGRAY, GRAY, fs=8.0, weight="bold")

    y_eng = 77.5
    box(ax, 18, y_eng, 32, 9.0,
        "Geometry engine\nthree feasibility constraints;\ninflated-obstacle shortest paths", LBLUE, BLUE, fs=8.0)
    box(ax, 50, y_eng, 27, 9.0,
        "Error model library\nbearing $\\sigma_b=1^\\circ$;\npointing $\\sigma_{loc}(s)$, $0.8^\\circ$, $0.6^\\circ$",
        LBLUE, BLUE, fs=8.0)
    box(ax, 81, y_eng, 30, 9.0,
        "Event-timeline simulator\ndrive $\\to$ slew $\\to$ joint observe\n$\\to$ reset (east, level)",
        LBLUE, BLUE, fs=8.0)
    for xx in (18, 50, 81):
        arrow(ax, xx, 89.7, xx, y_eng + 4.9)

    y_q = 57.5
    box(ax, 12.5, y_q, 22, 9.4,
        "Q1\nobservable region;\n$\\min r$ optimal station", LGREEN, GREEN, fs=8.0)
    box(ax, 37.5, y_q, 24, 9.4,
        "Q2(1)\nFisher accuracy model;\noptimal vehicle-station pair", LGREEN, GREEN, fs=8.0)
    box(ax, 63, y_q, 23, 9.4,
        "Q2(2)\nstructure enumeration;\nmin-max schedule", LGREEN, GREEN, fs=8.0)
    box(ax, 88, y_q, 21, 9.4,
        "Q3\nerror propagation;\n$s_{\\max}$, $(a,b)$ region", LGREEN, GREEN, fs=8.0)
    arrow(ax, 16, y_q, 25.5, y_q, lw=1.0)
    arrow(ax, 42.5, y_q, 51.5, y_q, lw=1.0)
    arrow(ax, 76.5, y_q + 1.8, 77.5, y_q + 1.8, style="--", color=ORANGE, lw=1.0)
    ax.text(91, 45.0, "deployment insight", fontsize=6.6, color=ORANGE, ha="center")

    y_q4 = 34
    box(ax, 50, y_q4, 46, 9.8,
        "Q4  event-triggered rescheduling\nscenarios: (a) 36, (a+c) 28, (b)/(c)/(b+c) 15, combinations 7\n"
        "$P_0$ retain  vs  $P_1$ re-optimize  vs  $P_2$ rule-based",
        LGREEN, GREEN, fs=8.0)
    arrow(ax, 61, y_q - 4.9, 56, y_q4 + 5.2, lw=1.0)
    ax.text(61.5, y_q4 + 8.0, "baseline plan", fontsize=6.8, color=GRAY, ha="center")
    arrow(ax, 88, y_q - 4.9, 70, y_q4 + 5.2, style="--", color=ORANGE, lw=1.0, rad=-0.15)

    y_out = 12
    box(ax, 50, y_out, 96, 10.0,
        "Final deliverables\nobservable-region maps; optimal stations & vehicle-station pair; per-vehicle timelines;\n"
        "accuracy / pointing-reliability bounds; resilience $(\\eta,t_{last})$ and value boundary of rescheduling",
        LGRAY, GRAY, fs=8.0, weight="bold")
    arrow(ax, 13, y_q - 4.9, 20, y_out + 5.3)
    arrow(ax, 38, y_q - 4.9, 38, y_out + 5.3)
    arrow(ax, 62, y_q - 4.9, 60, y_out + 5.3)
    arrow(ax, 50, y_q4 - 5.2, 50, y_out + 5.3)

    ax.text(2.0, 46.5, "Monte Carlo validation and\nindependent re-checks\napply to every question",
            fontsize=7.2, color=BLUE, ha="left", va="center", style="italic")
    arrow(ax, 7.5, 51.5, 11.0, y_q - 1.5, style="--", color=BLUE, lw=1.0, rad=0.2)
    arrow(ax, 7.5, 42.5, 15.5, y_out + 4.6, style="--", color=BLUE, lw=1.0, rad=-0.25)

    plt.tight_layout(pad=0.3)
    fig.savefig(f"{OUT}/fig_framework.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_schematic():
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.6),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})

    # ---- panel (a): field map ----
    ax = axes[0]
    ax.set_xlim(-2, 102)
    ax.set_ylim(-2, 102)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)", fontsize=9)
    ax.set_ylabel("y (m)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(True, lw=0.3, color="#cccccc", alpha=0.7)
    for i, b in enumerate(BUILDINGS, 1):
        ax.add_patch(Rectangle((b.x0, b.y0), b.x1 - b.x0, b.y1 - b.y0,
                               fc="#d9d9d9", ec=GRAY, lw=1.0))
        ax.add_patch(Rectangle((b.x0 - 1.5, b.y0 - 1.5),
                               b.x1 - b.x0 + 3, b.y1 - b.y0 + 3,
                               fc="none", ec=ORANGE, lw=0.9, ls="--"))
        ax.text((b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2 + 1.6, f"B{i}",
                ha="center", fontsize=8.5, weight="bold", color=GRAY)
        ax.text((b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2 - 2.4, f"h={b.h} m",
                ha="center", fontsize=6.8, color=GRAY)
    for vid, (vx, vy) in VEHICLE_INIT.items():
        ax.add_patch(Circle((vx, vy), 2.3, fc=LBLUE, ec=BLUE, lw=1.0, zorder=5))
        dx, dy = 0, -4.8
        if vid == 1:
            dx, dy = 5.0, 0.4
        if vid == 3:
            dx, dy = 4.8, 0.4
        if vid == 4:
            dy = -5.0
        ax.text(vx + dx, vy + dy, f"V{vid}", ha="center", fontsize=8,
                color=BLUE, weight="bold")
    tinfo = {"T1": (30, 45, 45), "T2": (150, 45, 35), "T3": (285, 38, 55)}
    lab_pos = {"T1": (72, 79.5), "T2": (11.0, 82.5), "T3": (63, 10.5)}
    for name, (th, R, H) in tinfo.items():
        t = TARGETS[name]
        ax.plot(t.x, t.y, marker="D", ms=6.5, mfc="#C0392B", mec="black",
                mew=0.7, zorder=6)
        ax.text(*lab_pos[name],
                f"{name} ($\\theta$={th}$^\\circ$, R={R} m, H={H} m)",
                ha="center", fontsize=7.4, color="#7B241C", weight="bold")
    ax.plot(*DISPATCH, marker="+", ms=10, color="black", mew=1.6, zorder=6)
    ax.text(DISPATCH[0] + 2.0, DISPATCH[1] + 1.5, "O (dispatch center)", fontsize=8)
    ax.plot([0, 11], [0, 0], color="black", lw=0.8)
    ax.plot([0, 0], [0, 11], color="black", lw=0.8)
    ax.text(12.2, -0.2, "east", fontsize=7.5, va="center")
    ax.text(-1.4, 12.2, "north", fontsize=7.5)
    ax.add_patch(Wedge((0, 0), 7.5, 0, 30, width=2.2, fc="#BDC3C7", ec="none"))
    ax.text(10.5, 5.6, "$\\theta$ (CCW)", fontsize=7.5, color=GRAY)
    ax.text(1.5, 97.5, "dashed outline: 1.5 m clearance inflation",
            fontsize=7, color=ORANGE, ha="left")
    ax.set_title("(a) Field, buildings, vehicles and alert targets", fontsize=9.5)

    # ---- panel (b): three constraints around T1 ----
    ax = axes[1]
    ax.set_xlim(-27, 36)
    ax.set_ylim(-27, 31)
    ax.set_aspect("equal")
    ax.axis("off")
    tx, ty = 0.0, 0.0
    r80 = 7.5821
    bx0, bx1, by0, by1 = 14, 27, -10, 6

    ax.add_patch(Polygon([(bx0, by0), (30, -21.4), (30, 12.9), (bx0, by1)],
                         closed=True, fc="#BDC3C7", ec="none", alpha=0.7, zorder=1))
    ax.text(23.5, -1.5, "occlusion shadow\n(tangency\n= blocked)", fontsize=7.3,
            ha="center", va="center", color="#333333")

    ax.add_patch(Circle((tx, ty), r80, fc="#F5B7B1", ec="none", alpha=0.9, zorder=2))
    ax.add_patch(Circle((tx, ty), r80, fc="none", ec="#C0392B", lw=1.1, ls=":",
                        zorder=5))

    ax.add_patch(Polygon([(bx0 - 2, by0 - 2), (bx1 + 2, by0 - 2), (bx1 + 2, by1 + 2),
                          (bx0 - 2, by1 + 2)], closed=True, fc=LORANGE, ec=ORANGE,
                         lw=1.0, ls="--", zorder=3))
    ax.add_patch(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fc="#d9d9d9",
                           ec=GRAY, lw=1.0, zorder=4))
    ax.text((bx0 + bx1) / 2, (by0 + by1) / 2, "building\n(B, h)", ha="center",
            va="center", fontsize=7.5, color=GRAY)

    ax.plot([tx, bx0], [ty, by0], color="#4A4A4A", lw=0.9, zorder=7)
    ax.plot([tx, bx0], [ty, by1], color="#4A4A4A", lw=0.9, zorder=7)

    ax.add_patch(Wedge((tx, ty), r80, 118, 162, width=2.0, fc=GREEN, ec="none",
                       zorder=6))
    ax.text(-17.5, 16.5, "feasible arc on $r=r_{80}$\n(min-$r$ attained)", fontsize=7.5,
            color=GREEN, ha="center", weight="bold")
    arrow(ax, -13.0, 13.0, -5.2, 6.6, color=GREEN, lw=1.0)

    ax.add_patch(Circle((tx, ty), 0.7, fc="#C0392B", ec="black", lw=0.6, zorder=8))
    ax.text(tx - 1.2, ty - 4.0, "T1 ground\nprojection", fontsize=7.5, ha="center",
            color="#7B241C", va="top")

    ax.text(-25.0, -13.0, "elevation infeasible:\n$r<r_{80}=\\dfrac{H-2}{\\tan 80^\\circ}$",
            fontsize=8, ha="left", color="#7B241C")
    arrow(ax, -6.5, -12.0, -4.8, -6.6, color="#C0392B", lw=1.0)

    ax.text(31.5, 24.0, "dashed: clearance\n$\\geq 1.5$ m", fontsize=7.5, color=ORANGE,
            ha="center")
    arrow(ax, 30.0, 20.5, 28.3, 9.0, color=ORANGE, lw=1.0)

    ax.text(-25.5, 27.0, "camera at disc center, $h_c=2$ m;\n"
            "elevation $e=\\arctan((H-h_c)/r)$", fontsize=7.5, ha="left")
    ax.set_title("(b) Three feasibility constraints (not to scale)", fontsize=9.5)

    plt.tight_layout(pad=0.5)
    fig.savefig(f"{OUT}/fig_geometry_schematic.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_framework()
    fig_schematic()
    print("done")
