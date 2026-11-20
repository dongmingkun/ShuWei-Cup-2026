# -*- coding: utf-8 -*-
"""B题问题3: 三误差源 -> Δψ 传播模型 + 最大安全行程 + (a,b) 约束域.

口径(框架 §4.4, 题面给定):
- 位置误差 (δx,δy) 逐轴独立 N(0, σ_loc²(s)), σ_loc(s)=0.3+0.02s m, 各向同性;
- 云台零偏 0.8°、控制 0.6° 独立, 直接以度加到方位角上;
- 目标已知, 水平距离固定 r=40 m; Δψ 卷绕到 [-180°,180°);
- 可靠性约束 P(|Δψ|≤5°) ≥ 0.95.

模型链(M1 解析主线):
- δ⊥ = δ·n ~ N(0, σ_loc²) (各向同性, 任意视线法向同方差);
- Δψ_loc ≈ (180/π)·δ⊥/r  =>  σ₁(s) = (180/π)·σ_loc(s)/40;
- Δψ ≈ N(0, σ_tot²), σ_tot² = σ₁² + 0.8² + 0.6²;
- P(|Δψ|≤5°) = 2Φ(5/σ_tot)−1 ≥ 0.95  ⇔  σ_tot ≤ 5/z₀.₉₇₅.

M2 验证轨: 蒙特卡洛不做线性化——精确 atan2 + wrap180 + 三源采样.
M3 附注: Chebyshev / 对称 Cantelli 分布无关下界(仅供稳健性讨论).

声明的建模假设(论文假设章须写明): 小角近似、卷绕在 ±5° 量级不起作用、
云台两角误差按度直接相加、s 为自出发起的累计行程(不引入重定位清零).
"""
import json
import math

import numpy as np
from scipy.stats import norm

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = __import__("pathlib").Path(__file__).parent / "q3_output"
OUT.mkdir(exist_ok=True)

R_HORIZ = 40.0        # m, 题面给定水平距离
SIGMA_BIAS = 0.8      # °, 云台零偏
SIGMA_CTRL = 0.6      # °, 云台控制
SIGMA_LOC0, SIGMA_LOC_K = 0.3, 0.02   # σ_loc(s) = 0.3 + 0.02s m
PSI_TOL = 5.0         # °
RELIABILITY = 0.95
Z975 = float(norm.ppf(0.975))
DEG = 180.0 / math.pi
S_GRID = np.arange(0.0, 100.0 + 1e-9, 1.0)
N_SWEEP_MC = 100_000
N_HEADLINE_MC = 400_000
RNG_SEED = 2026


def sigma_loc(s):
    return SIGMA_LOC0 + SIGMA_LOC_K * s


def sigma1(s):
    """位置误差项的角度标准差(度): σ₁ = (180/π)·σ_loc(s)/r."""
    return DEG * sigma_loc(s) / R_HORIZ


def sigma_tot(s):
    return np.sqrt(sigma1(s) ** 2 + SIGMA_BIAS ** 2 + SIGMA_CTRL ** 2)


def p_analytic(s):
    return 2.0 * norm.cdf(PSI_TOL / sigma_tot(s)) - 1.0


def s_max_analytic():
    """σ_tot ≤ 5/z₀.₉₇₅ → σ₁_max → σ_loc_max → s_max."""
    stot_max = PSI_TOL / Z975
    s1_max = math.sqrt(stot_max ** 2 - SIGMA_BIAS ** 2 - SIGMA_CTRL ** 2)
    sloc_max = s1_max * R_HORIZ / DEG
    return stot_max, s1_max, sloc_max, (sloc_max - SIGMA_LOC0) / SIGMA_LOC_K


def wilson(p_hat, n, z=Z975):
    ctr = p_hat + z * z / (2 * n)
    half = z * math.sqrt(p_hat * (1 - p_hat) / n + z * z / (4 * n * n))
    den = 1 + z * z / n
    return (ctr - half) / den, (ctr + half) / den


def mc_batch(s, n, rng):
    """一次蒙特卡洛: 精确 atan2 + wrap180 + 三源采样.

    布局: 目标 T=(40,0), 真实位置 P=(0,0), 估计位置 P̂=P+δ.
    真实方位 φ=0; 指令方位 ψ̂=atan2(−δy, 40−δx); Δψ=wrap180(ψ̂−φ)+ε₀+εc.
    同时返回线性化几何项 (180/π)·δy/40 供近似误差核算.
    """
    sloc = sigma_loc(s)
    dx = rng.normal(0.0, sloc, n)
    dy = rng.normal(0.0, sloc, n)
    e0 = rng.normal(0.0, SIGMA_BIAS, n)
    ec = rng.normal(0.0, SIGMA_CTRL, n)
    dpsi_exact = np.degrees(np.arctan2(-dy, R_HORIZ - dx))
    dpsi_lin = -dy / R_HORIZ * DEG
    dpsi_total = (dpsi_exact + e0 + ec + 180.0) % 360.0 - 180.0
    big = np.abs(dpsi_exact) >= 1.0
    ex_rad = np.radians(dpsi_exact[big])
    lin_rad = np.radians(dpsi_lin[big])
    rel = np.abs(ex_rad - lin_rad) / np.abs(ex_rad)
    return {
        "hit_exact": np.abs(dpsi_total) <= PSI_TOL,
        "hit_lin": np.abs(dpsi_lin + e0 + ec) <= PSI_TOL,
        "geom_diff_deg": np.abs(dpsi_total - (dpsi_lin + e0 + ec)),
        "lin_err_rel": rel,
    }


def mc_sweep(rng):
    ps, los = [], []
    for s in S_GRID:
        b = mc_batch(float(s), N_SWEEP_MC, rng)
        ps.append(float(b["hit_exact"].mean()))
        los.append(float(b["hit_lin"].mean()))
    return np.array(ps), np.array(los)


def mc_smax_interp(ps):
    """P_MC–s 曲线上 0.95 水平的线性插值交点(下降沿)."""
    below = np.where(ps < RELIABILITY)[0]
    if len(below) == 0:
        return None
    i = below[0]
    if i == 0:
        return float(S_GRID[0])
    s0, s1 = S_GRID[i - 1], S_GRID[i]
    p0, p1 = ps[i - 1], ps[i]
    return float(s0 + (RELIABILITY - p0) / (p1 - p0) * (s1 - s0))


def conservative_bounds():
    """M3: 分布无关下界反推 (仅二阶矩).

    Chebyshev: P(|X|≤5) ≥ 1 − σ_tot²/25 ≥ 0.95 → σ_tot ≤ 5√0.05;
    对称 Cantelli: P(|X|≥5) ≤ 2σ²/(σ²+25) ≤ 0.05 → σ_tot² ≤ 25·0.05/1.95.
    """
    stot_cheb = PSI_TOL * math.sqrt(1.0 - RELIABILITY)
    s1_cheb = math.sqrt(stot_cheb ** 2 - 1.0)
    sloc_cheb = s1_cheb * R_HORIZ / DEG
    s_cheb = (sloc_cheb - SIGMA_LOC0) / SIGMA_LOC_K
    var_cant = PSI_TOL ** 2 * (1 - RELIABILITY) / (2 - (1 - RELIABILITY))
    return {
        "chebyshev": {"sigma_tot_max_deg": stot_cheb,
                      "s1_max_deg": s1_cheb,
                      "sigma_loc_max_m": sloc_cheb, "s_max_m": s_cheb},
        "cantelli": {"sigma_tot_max_deg": math.sqrt(var_cant),
                     "feasible_with_fixed_gimbal": var_cant >= 1.0,
                     "note": "σ_tot² ≤ 0.641°² < 云台零偏+控制方差 1.0°², "
                             "即使定位误差为零也不可行"},
    }


def contribution(s):
    """三源方差占比; s 可为标量或数组(JSON 侧自行转 float)."""
    st2 = sigma_tot(s) ** 2
    return {"loc": sigma1(s) ** 2 / st2,
            "gimbal_bias": SIGMA_BIAS ** 2 / st2,
            "gimbal_ctrl": SIGMA_CTRL ** 2 / st2}


def crossover(g1, g2):
    """σ₁(s)² = g² 的行程交点(定位项方差追上固定项)."""
    sloc = g1 * R_HORIZ / DEG
    return (sloc - SIGMA_LOC0) / SIGMA_LOC_K


def main():
    rng = np.random.default_rng(RNG_SEED)
    stot_max, s1_max, sloc_max, s_max = s_max_analytic()
    print(f"z0.975 = {Z975:.6f}")
    print(f"sigma_tot_max = {stot_max:.6f} deg")
    print(f"sigma1_max = {s1_max:.6f} deg, sigma_loc_max = {sloc_max:.6f} m")
    print(f"s_max = {s_max:.4f} m")

    # 线性化相对误差上界: |atan x − x|/|x| ≈ x²/3, x = σ₁_max 弧度
    x_rad = s1_max / DEG
    lin_rel = x_rad ** 2 / 3.0

    ps_mc, ps_lin = mc_sweep(rng)
    mc_sm = mc_smax_interp(ps_mc)
    print(f"MC sweep: max|P_mc − P_ana| = "
          f"{np.max(np.abs(ps_mc - p_analytic(S_GRID))):.5f}")
    print(f"MC s_max(插值) = {mc_sm:.3f} m")

    # 主验证点 s = s_max: 大样本 + Wilson CI
    b = mc_batch(s_max, N_HEADLINE_MC, rng)
    p_head = float(b["hit_exact"].mean())
    lo, hi = wilson(p_head, N_HEADLINE_MC)
    p_head_lin = float(b["hit_lin"].mean())
    print(f"MC@s_max: P_exact = {p_head:.5f} [{lo:.5f},{hi:.5f}], "
          f"P_linearized = {p_head_lin:.5f}")

    # 近似误差实证: 线性化几何项与精确 atan2 的最大绝对差
    print(f"lin vs exact: max abs geom diff = {b['geom_diff_deg'].max():.3e} deg, "
          f"max rel err = {b['lin_err_rel'].max():.3e}")

    m3 = conservative_bounds()
    print(f"M3 Chebyshev s_max = {m3['chebyshev']['s_max_m']:.3f} m; "
          f"Cantelli 可行? {m3['cantelli']['feasible_with_fixed_gimbal']}")

    cross_gb = crossover(SIGMA_BIAS, None)
    cross_gc = crossover(SIGMA_CTRL, None)
    print(f"定位项方差追上云台零偏/控制: s = {cross_gb:.3f} / {cross_gc:.3f} m")

    # (a,b) 约束域: a + 80b ≤ σ_loc_max, a,b ≥ 0
    ab_area = sloc_max ** 2 / 160.0
    baseline_in = SIGMA_LOC0 + 80 * SIGMA_LOC_K <= sloc_max + 1e-12

    # wrap 敏感性: |几何方位差| > 180° 的概率上界(正态尾, 双精度下溢记 0)
    wrap_p = float(math.erfc((180.0 / sigma1(80.0)) / math.sqrt(2.0)))

    summary = {
        "z975": round(Z975, 6),
        "chain": {"sigma_tot_max_deg": round(stot_max, 6),
                  "sigma1_max_deg": round(s1_max, 6),
                  "sigma_loc_max_m": round(sloc_max, 6),
                  "s_max_m": round(s_max, 4)},
        "model": "σ_tot² = [(180/π)·(0.3+0.02s)/40]² + 0.8² + 0.6²; "
                 "P(|Δψ|≤5°) = 2Φ(5/σ_tot)−1 ≥ 0.95",
        "linearization": {
            "theory_rel_err_upper_1sigma": lin_rel,
            "mc_max_rel_err_deg_above_1deg": float(b["lin_err_rel"].max()),
            "mc_max_abs_geom_diff_deg": float(b["geom_diff_deg"].max()),
            "mc_rel_err_filter": "|Δψ| < 1° 的样本不参与相对误差统计 "
                                 "(比值在零点附近无意义)",
            "p_at_smax_exact": p_head, "p_at_smax_wilson95": [lo, hi],
            "p_at_smax_linearized": p_head_lin,
        },
        "mc_validation": {
            "n_per_point": N_SWEEP_MC, "seed": RNG_SEED,
            "max_abs_p_diff_over_grid": float(np.max(np.abs(ps_mc - p_analytic(S_GRID)))),
            "s_max_mc_interp_m": mc_sm,
            "sweep": [{"s": float(s), "p_analytic": float(pa), "p_mc": float(pm)}
                      for s, pa, pm in zip(S_GRID, p_analytic(S_GRID), ps_mc)],
        },
        "conservative_bounds_m3": m3,
        "contribution_crossovers_m": {"loc_vs_gimbal_bias": round(cross_gb, 3),
                                      "loc_vs_gimbal_ctrl": round(cross_gc, 3)},
        "contribution_at": {str(s): {k: float(v)
                                     for k, v in contribution(s).items()}
                            for s in (0, 40, 80)},
        "ab_region": {"constraint": "a + 80b ≤ σ_loc_max, a,b ≥ 0",
                      "sigma_loc_max_m": round(sloc_max, 6),
                      "vertices": [[0, 0], [round(sloc_max, 6), 0], [0, round(sloc_max / 80, 6)]],
                      "area_m2": round(ab_area, 6),
                      "baseline_0.3+0.02s_inside": baseline_in,
                      "baseline_note": "题面基线 (a,b)=(0.3,0.02) 在 80 m 处超差, "
                                       "与 s_max≈66.9<80 一致"},
        "wrap_sensitivity": {"p_geom_gt_180deg_at_s80": wrap_p,
                             "note": "双精度下溢为 0, 卷绕在可行域内不起作用"},
        "cross_question": {
            "q2_1_travel_budget_m": 90.0,
            "note": "指向可靠性行程上限 ≈66.93 m < Q2(1) 90 m 行驶预算, "
                    "可靠性预算先约束; 部署应控制单车累计行程或中途换车",
        },
    }
    with open(OUT / "q3_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    _plot_p_s(s_max, ps_mc, m3)
    _plot_contribution(s_max, cross_gb, cross_gc)
    _plot_ab(sloc_max)
    print(f"输出已写入 {OUT}")


def _plot_p_s(s_max, ps_mc, m3):
    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    ss = np.linspace(0, 100, 1001)
    ax.plot(ss, p_analytic(ss), lw=2, color="#1f618d",
            label=r"analytic  $2\Phi(5/\sigma_{tot})-1$")
    ax.plot(S_GRID, ps_mc, "o", ms=3.5, color="#e67e22", alpha=0.75,
            label=f"MC exact ({N_SWEEP_MC:,}/pt)")
    st = sigma_tot(ss)
    ax.plot(ss, 1 - st ** 2 / PSI_TOL ** 2, ls="--", lw=1.2, color="#7f8c8d",
            label="Chebyshev lower bound")
    ax.plot(ss, 1 - 2 * st ** 2 / (st ** 2 + PSI_TOL ** 2), ls=":", lw=1.2,
            color="#95a5a6", label="symmetric Cantelli lower bound")
    ax.axhline(RELIABILITY, color="crimson", lw=1.0, ls="--")
    ax.axvline(s_max, color="crimson", lw=1.0, ls="--")
    ax.text(s_max + 1.2, 0.62, f"$s_{{max}}$ = {s_max:.2f} m", color="crimson")
    ax.text(2, RELIABILITY - 0.045, "0.95", color="crimson")
    ax.set_xlabel("cumulative travel distance s (m)")
    ax.set_ylabel(r"$P(|\Delta\psi|\leq 5°)$")
    ax.set_ylim(0.5, 1.02)
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title("Q3 reliability vs cumulative travel: analytic, MC and "
                 "distribution-free bounds", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "q3_p_s_curve.png", dpi=300)
    plt.close(fig)


def _plot_contribution(s_max, cross_gb, cross_gc):
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    ss = np.linspace(0, 100, 1001)
    c = contribution(ss)
    ax.stackplot(ss, c["loc"], c["gimbal_bias"], c["gimbal_ctrl"],
                 labels=["positioning (σ₁)", "gimbal bias (0.8°)",
                         "gimbal control (0.6°)"],
                 colors=["#5dade2", "#f5b041", "#58d68d"], alpha=0.9)
    for x, txt in ((cross_gc, f"s = {cross_gc:.1f} m\nloc = ctrl"),
                   (cross_gb, f"s = {cross_gb:.1f} m\nloc = bias")):
        ax.axvline(x, color="k", lw=0.8, ls=":")
        ax.text(x + 1, 0.52, txt, fontsize=8)
    ax.axvline(s_max, color="crimson", lw=1.2, ls="--")
    ax.text(s_max + 1.2, 0.08, f"$s_{{max}}$ = {s_max:.1f} m", color="crimson")
    ax.set_xlabel("cumulative travel distance s (m)")
    ax.set_ylabel("variance share of σ_tot²")
    ax.set_ylim(0, 1.0)
    ax.set_xlim(0, 100)
    ax.legend(loc="center right", fontsize=8)
    ax.set_title("Q3 variance contribution decomposition", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "q3_contribution.png", dpi=300)
    plt.close(fig)


def _plot_ab(sloc_max):
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.9))
    ax = axes[0]
    bmax = sloc_max / 80.0
    tri = plt.Polygon([[0, 0], [sloc_max, 0], [0, bmax]],
                      closed=True, color="#abeba1", alpha=0.8)
    ax.add_patch(tri)
    xs = np.array([0.0, sloc_max])
    ax.plot(xs, (sloc_max - xs) / 80.0, color="#1e8449", lw=1.6,
            label=r"$a+80b=\sigma_{loc}^{max}$")
    ax.plot(0.3, 0.02, "x", color="crimson", ms=9, mew=2)
    ax.annotate(r"baseline $(0.3, 0.02)$" + "\ninfeasible at s=80 m",
                (0.3, 0.02), (0.42, 0.0165), fontsize=8, color="crimson",
                arrowprops=dict(arrowstyle="->", color="crimson", lw=0.8))
    ax.plot(0.3, 0.008, "o", color="#1f618d", ms=5)
    ax.annotate("e.g. $(0.3, 0.008)$ feasible", (0.3, 0.008), (0.5, 0.006),
                fontsize=8, color="#1f618d")
    ax.set_xlim(-0.02, sloc_max * 1.18)
    ax.set_ylim(-0.0006, bmax * 1.25)
    ax.set_xlabel("intercept a (m)")
    ax.set_ylabel("slope b (m/m)")
    ax.set_title(f"(a,b) feasible region (area = {sloc_max ** 2 / 160:.4f} m²)",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper right")

    ax2 = axes[1]
    labels = ["pointing reliability\nbudget (Q3)", "travel budget\n(Q2(1))"]
    vals = [(sloc_max - SIGMA_LOC0) / SIGMA_LOC_K, 90.0]
    bars = ax2.barh(labels, vals, color=["#e67e22", "#5dade2"], height=0.5)
    for bar, v in zip(bars, vals):
        ax2.text(v + 1.5, bar.get_y() + bar.get_height() / 2, f"{v:.1f} m",
                 va="center", fontsize=9)
    ax2.set_xlim(0, 105)
    ax2.set_xlabel("maximum cumulative travel (m)")
    ax2.set_title("Q3 reliability budget vs Q2(1) travel budget", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "q3_ab_region.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
