# -*- coding: utf-8 -*-
"""Q2(1): 双站纯方位定位精度模型 + T2 最优车-站对求解 + MC/退化验证.

运行: python -X utf8 -B q2_1_localization.py
输出: q2_output/ 下 PNG 图件与 q2_1_summary.json, 控制台打印五项关键输出.

模型口径(与框架 §4.2 一致):
- 观测方程 phi_i = atan2(Ty-s_iy, Tx-s_ix) + eps_i, eps_i ~ N(0, sigma_b^2) 独立,
  sigma_b = 1° = pi/180 rad;
- H_i = d(phi)/d(Tx,Ty) = u_i^T/d_i, u_i = (-sin phi_i, cos phi_i) 为横向单位矢量;
- Fisher I = sigma^-2 * sum(u_i u_i^T / d_i^2), C = I^-1;
- 主指标 RMSE = sqrt(tr C) (A-最优), 辅指标 det C (D-最优) 与 cond(C) (各向异性);
- 解析下界 RMSE^2 = sigma^2 (d1^2+d2^2)/sin^2(gamma) >= 2 sigma^2 r80^2,
  等号 iff d1=d2=r80 且 gamma=90° —— 全局最优可用两站几何直接证明;
- 约束: 站点 ∈ T2 可观测域(geometry.py 三约束), 每车初始位->站 绕障行程 <= 90 m
  (30 s x 3 m/s), 行程用 shortest_path 精确核算(有限图上界性质见交接文档 §5).
"""
import heapq
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from scipy.optimize import minimize

from geometry import (BUILDINGS, FIELD_SIZE, TARGETS, VEHICLE_INIT,
                      VEHICLE_RADIUS, shortest_path, is_station_feasible,
                      elevation_angle)
from q1_observable_domain import region_mask

OUT = Path(__file__).with_name("q2_output")
OUT.mkdir(exist_ok=True)

T2T = TARGETS["T2"]
T2 = (T2T.x, T2T.y)                  # 水平投影坐标元组
SIGMA_B = math.radians(1.0)          # 1° 方位标准差
BUDGET_M = 30.0 * 3.0                # 行程预算 90 m
RES = 0.25                           # 路网/筛选网格 (m)


# ---------------------------------------------------------------- 解析精度模型
def bearing_terms(s, t=T2):
    """站点 s 对目标投影 t 的 (phi, u=(ux,uy), w=1/d^2)."""
    dx, dy = t[0] - s[0], t[1] - s[1]
    d = math.hypot(dx, dy)
    phi = math.atan2(dy, dx)
    return phi, (-math.sin(phi), math.cos(phi)), 1.0 / d**2


def rmse_analytic(s1, s2, sigma=SIGMA_B):
    """sqrt(tr C) 闭式: sigma^2 (w1+w2)/(detM), detM = w1 w2 sin^2(gamma)."""
    _, u1, w1 = bearing_terms(s1)
    _, u2, w2 = bearing_terms(s2)
    a11 = w1 * u1[0] ** 2 + w2 * u2[0] ** 2
    a22 = w1 * u1[1] ** 2 + w2 * u2[1] ** 2
    a12 = w1 * u1[0] * u1[1] + w2 * u2[0] * u2[1]
    det = a11 * a22 - a12 * a12
    if det <= 1e-18:
        return float("inf")
    return sigma * math.sqrt((a11 + a22) / det)


def metrics_full(s1, s2, sigma=SIGMA_B):
    """返回 RMSE, detC, condC, 误差椭圆半轴(sqrt 特征值)."""
    _, u1, w1 = bearing_terms(s1)
    _, u2, w2 = bearing_terms(s2)
    M = np.array([[w1 * u1[0] ** 2 + w2 * u2[0] ** 2,
                   w1 * u1[0] * u1[1] + w2 * u2[0] * u2[1]],
                  [w1 * u1[0] * u1[1] + w2 * u2[0] * u2[1],
                   w1 * u1[1] ** 2 + w2 * u2[1] ** 2]])
    C = sigma**2 * np.linalg.inv(M)
    ev = np.linalg.eigvalsh(C)
    tr_c = float(np.trace(C))
    det_c = float(np.linalg.det(C))
    cond_c = float(ev.max() / ev.min())
    axes = (math.sqrt(ev.max()), math.sqrt(ev.min()))
    return {"rmse": math.sqrt(tr_c), "detC": det_c, "condC": cond_c,
            "axes_1sigma": axes}


def crossing_angle(s1, s2):
    """两条视线交角 gamma ∈ [0°,180°]: |wrap180(phi2-phi1)|."""
    p1, _, _ = bearing_terms(s1)
    p2, _, _ = bearing_terms(s2)
    d = math.degrees(p2 - p1) % 360.0
    return d if d <= 180.0 else 360.0 - d


# ---------------------------------------------------------------- 可达域(路网 Dijkstra)
def _clearance_mask(xs):
    n = len(xs)
    X, Y = np.meshgrid(xs, xs)
    m = ((X >= VEHICLE_RADIUS - 1e-9) & (X <= FIELD_SIZE - VEHICLE_RADIUS + 1e-9)
         & (Y >= VEHICLE_RADIUS - 1e-9) & (Y <= FIELD_SIZE - VEHICLE_RADIUS + 1e-9))
    for b in BUILDINGS:
        ddx = np.maximum(np.maximum(b.x0 - X, X - b.x1), 0.0)
        ddy = np.maximum(np.maximum(b.y0 - Y, Y - b.y1), 0.0)
        m &= np.hypot(ddx, ddy) >= VEHICLE_RADIUS - 1e-9
    return m


def geodesic_field(road, start_xy, xs):
    """8 连通(禁斜切角)网格 Dijkstra: 从 start_xy 到各格的绕障最短路长(m)."""
    n = len(xs)
    si = int(round(np.searchsorted(xs, start_xy[1])))
    sj = int(round(np.searchsorted(xs, start_xy[0])))
    si = min(max(si, 0), n - 1); sj = min(max(sj, 0), n - 1)
    if not road[si, sj]:
        raise RuntimeError(f"起点 {start_xy} 不在路网掩码内")
    dist = np.full(road.shape, np.inf)
    dist[si, sj] = 0.0
    sq2 = math.sqrt(2.0) * RES
    heap = [(0.0, si, sj)]
    while heap:
        d, i, j = heapq.heappop(heap)
        if d > dist[i, j] + 1e-12:
            continue
        for di, dj, w in ((-1, 0, RES), (1, 0, RES), (0, -1, RES), (0, 1, RES),
                          (-1, -1, sq2), (-1, 1, sq2), (1, -1, sq2), (1, 1, sq2)):
            ni, nj = i + di, j + dj
            if not (0 <= ni < n and 0 <= nj < n) or not road[ni, nj]:
                continue
            if di and dj and not (road[i + di, j] and road[i, j + dj]):
                continue                       # 对角线不许切角
            nd = d + w
            if nd < dist[ni, nj] - 1e-12:
                dist[ni, nj] = nd
                heapq.heappush(heap, (nd, ni, nj))
    return dist


def lookup_field(dist, xs, p):
    i = min(max(int(round(np.searchsorted(xs, p[1]))), 0), len(xs) - 1)
    j = min(max(int(round(np.searchsorted(xs, p[0]))), 0), len(xs) - 1)
    return float(dist[i, j])


# ---------------------------------------------------------------- 批量 RMSE 筛选
def _pair_rmse_chunk(ua_x, ua_y, wa, ub_x, ub_y, wb, sigma):
    """ua/wa: (n,1) 广播列, ub/wb: (1,m) 广播行 -> (n,m) RMSE 矩阵."""
    a11 = wa * ua_x**2 + wb * ub_x**2
    a22 = wa * ua_y**2 + wb * ub_y**2
    a12 = wa * ua_x * ua_y + wb * ub_x * ub_y
    det = a11 * a22 - a12 * a12
    det = np.where(det > 1e-18, det, np.nan)
    return sigma * np.sqrt((a11 + a22) / det)


def screen_pair(cand_a, cand_b, sigma=SIGMA_B):
    """cand_*: dict(ux, uy, w, xy) -> (best_rmse, (ia, ib))."""
    n, m = len(cand_a["w"]), len(cand_b["w"])
    ua_x = cand_a["ux"][:, None]; ua_y = cand_a["uy"][:, None]
    wa = cand_a["w"][:, None]
    ub_x = cand_b["ux"][None, :]; ub_y = cand_b["uy"][None, :]
    wb = cand_b["w"][None, :]
    best = (float("inf"), (0, 0))
    step = 256
    for i0 in range(0, n, step):
        r = _pair_rmse_chunk(ua_x[i0:i0 + step], ua_y[i0:i0 + step],
                             wa[i0:i0 + step], ub_x, ub_y, wb, sigma)
        k = np.nanargmin(r)
        ki, kj = np.unravel_index(k, r.shape)
        v = float(r[ki, kj])
        if v < best[0]:
            best = (v, (i0 + ki, kj))
    return best


# ---------------------------------------------------------------- MC 验证
def mc_rmse(s1, s2, n=200000, sigma=SIGMA_B, seed=2026):
    """精确直线交点的蒙特卡洛位置 RMSE (M2 验证 M1 线性化)."""
    rng = np.random.default_rng(seed)
    p1, _, _ = bearing_terms(s1); p2, _, _ = bearing_terms(s2)
    eps = rng.normal(0.0, sigma, (n, 2))
    phi1 = p1 + eps[:, 0]; phi2 = p2 + eps[:, 1]
    v1x, v1y = np.cos(phi1), np.sin(phi1)
    v2x, v2y = np.cos(phi2), np.sin(phi2)
    cross = v1x * v2y - v1y * v2x
    ok = np.abs(cross) > 1e-12
    t1 = ((s2[0] - s1[0]) * v2y - (s2[1] - s1[1]) * v2x) / np.where(ok, cross, 1.0)
    px = s1[0] + t1 * v1x
    py = s1[1] + t1 * v1y
    err2 = (px - T2[0]) ** 2 + (py - T2[1]) ** 2
    return math.sqrt(float(err2[ok].mean()))


# ---------------------------------------------------------------- 精化
def refine(s1, s2, geod_a, geod_b, xs, veh_a, veh_b, budget=BUDGET_M):
    """Nelder-Mead 精化: 罚函数处理可观测域与行程约束."""
    def geod_penalty(veh, s):
        d = lookup_field(geod_a if veh == veh_a else geod_b, xs, s)
        return 0.0 if d <= budget else (d - budget)

    def obj(x):
        sa, sb = (x[0], x[1]), (x[2], x[3])
        viol = 0.0
        for s, veh in ((sa, veh_a), (sb, veh_b)):
            if not is_station_feasible(s, T2T):
                viol += 5.0
            viol += geod_penalty(veh, s)
        if viol:
            return 10.0 + 100.0 * viol
        return rmse_analytic(sa, sb)

    res = minimize(obj, np.array([*s1, *s2]), method="Nelder-Mead",
                   options={"maxiter": 1200, "xatol": 1e-4, "fatol": 1e-10})
    sa, sb = (res.x[0], res.x[1]), (res.x[2], res.x[3])
    feas = (is_station_feasible(sa, T2T) and is_station_feasible(sb, T2T))
    return sa, sb, feas


# ---------------------------------------------------------------- 绘图
def plot_reachability(ax, xs, mask_region, dist, veh_id, chosen, budget):
    ax.imshow(np.where(mask_region, 0.5, np.nan), origin="lower",
              extent=(0, FIELD_SIZE, 0, FIELD_SIZE),
              cmap=matplotlib.colors.ListedColormap(["#eaf2f8"]),
              vmin=0, vmax=1, interpolation="nearest")
    im = np.where(np.isfinite(dist) & (dist <= budget) & mask_region, 1.0, np.nan)
    ax.imshow(im, origin="lower", extent=(0, FIELD_SIZE, 0, FIELD_SIZE),
              cmap=matplotlib.colors.ListedColormap(["#7fb3d5"]),
              vmin=0, vmax=1, interpolation="nearest")
    for k, b in enumerate(BUILDINGS, start=1):
        ax.add_patch(Rectangle((b.x0, b.y0), b.x1 - b.x0, b.y1 - b.y0,
                               facecolor="#9e9e9e", edgecolor="k", zorder=3))
        ax.text((b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2, f"B{k}", ha="center",
                va="center", fontsize=7, zorder=4)
    ax.plot(T2[0], T2[1], "^", color="crimson", ms=8, zorder=5)
    ax.plot(*VEHICLE_INIT[veh_id], "s", color="darkorange", ms=6, zorder=5)
    ax.annotate(f"V{veh_id}", VEHICLE_INIT[veh_id], textcoords="offset points",
                xytext=(4, 4), fontsize=8, zorder=5)
    for vid, s in chosen.items():
        ax.plot(*s, "*", color="royalblue", ms=13, zorder=6)
    ax.set_xlim(0, FIELD_SIZE); ax.set_ylim(0, FIELD_SIZE)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"V{veh_id} reachable & feasible for T2", fontsize=9)


def plot_gamma_sweep(ax, r80, sigma, chosen_gamma):
    g = np.linspace(1.0, 179.0, 721)
    rmse_g = sigma * math.sqrt(2 * r80**2) / np.sin(np.radians(g))
    cond_g = (1 + np.abs(np.cos(np.radians(g)))) / (1 - np.abs(np.cos(np.radians(g))))
    ax.plot(g, rmse_g, color="royalblue", label="RMSE (d$_1$=d$_2$=r$_{80}$)")
    ax.axhline(sigma * math.sqrt(2) * r80, color="gray", lw=0.8, ls=":")
    ax.axvline(chosen_gamma, color="crimson", lw=1.0,
               label=f"chosen $\\gamma$={chosen_gamma:.1f}°")
    ax.axvspan(1, 30, color="#f5b7b1", alpha=0.35)
    ax.axvspan(150, 179, color="#f5b7b1", alpha=0.35)
    ax.set_xlabel("crossing angle $\\gamma$ (deg)")
    ax.set_ylabel("RMSE (m)")
    ax.set_yscale("log")
    ax2 = ax.twinx()
    ax2.plot(g, cond_g, color="seagreen", ls="--", label="cond(C)")
    ax2.set_ylabel("cond(C)")
    ax2.set_yscale("log")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper center")
    ax.set_title("Degeneracy sweep: RMSE & conditioning vs $\\gamma$", fontsize=10)


# ---------------------------------------------------------------- main
def main():
    from geometry import elev_min_radius
    r80 = elev_min_radius(T2T)
    print(f"T2 投影 {T2[0]:.4f},{T2[1]:.4f}, H={T2T.H} m, r80={r80:.4f} m")
    print(f"解析下界 RMSE* = sigma*sqrt(2)*r80 = {SIGMA_B * math.sqrt(2) * r80:.5f} m")

    xs = np.linspace(0.0, FIELD_SIZE, int(round(FIELD_SIZE / RES)) + 1)
    road = _clearance_mask(xs)
    _, region = region_mask(T2T, RES)         # T2 可观测域(含俯仰/LOS/净距)

    # --- 每车可达 ∩ 可观测域候选
    cands = {}
    reach_info = {}
    for vid, p in VEHICLE_INIT.items():
        dist = geodesic_field(road, p, xs)
        ok = region & np.isfinite(dist) & (dist <= BUDGET_M)
        idx = np.argwhere(ok)
        reach_info[vid] = int(len(idx))
        print(f"V{vid}: 可达∩可观测域格点 {len(idx)} (网格 {RES} m)")
        if len(idx) == 0:
            cands[vid] = None
            continue
        ii, jj = idx[:, 0], idx[:, 1]
        rng = np.random.default_rng(2026)
        if len(idx) > 900:                    # 筛选子采样(等距优先)
            stride = max(1, int(math.sqrt(len(idx) / 900)))
            sel = (ii % stride == 0) & (jj % stride == 0)
            if sel.sum() > 1200:
                keep = rng.choice(np.flatnonzero(sel), 1200, replace=False)
                selm = np.zeros_like(sel); selm[keep] = True; sel = selm
            ii, jj = ii[sel], jj[sel]
        sx, sy = xs[jj], xs[ii]
        phi = np.arctan2(T2[1] - sy, T2[0] - sx)
        cands[vid] = {"xy": np.column_stack([sx, sy]),
                      "ux": -np.sin(phi), "uy": np.cos(phi),
                      "w": 1.0 / ((T2[0] - sx) ** 2 + (T2[1] - sy) ** 2)}

    # --- 10 车对粗筛
    vids = [v for v in VEHICLE_INIT if cands[v] is not None]
    print(f"\n可达车辆: {vids}  (不可达: {[v for v in VEHICLE_INIT if v not in vids]})")
    coarse = []
    for i in range(len(vids)):
        for j in range(i + 1, len(vids)):
            va, vb = vids[i], vids[j]
            rmse, (ia, ib) = screen_pair(cands[va], cands[vb])
            coarse.append({"pair": (va, vb), "rmse": rmse,
                           "s1": tuple(cands[va]["xy"][ia]),
                           "s2": tuple(cands[vb]["xy"][ib])})
            print(f"  粗筛 V{va}+V{vb}: RMSE={rmse:.5f} m")
    coarse.sort(key=lambda c: c["rmse"])

    # --- 前 12 组精化 + 精确行程核算
    geods = {}
    refined = []
    for c in coarse[:12]:
        va, vb = c["pair"]
        if va not in geods:
            geods[va] = geodesic_field(road, VEHICLE_INIT[va], xs)
        if vb not in geods:
            geods[vb] = geodesic_field(road, VEHICLE_INIT[vb], xs)
        sa, sb, feas = refine(c["s1"], c["s2"], geods[va], geods[vb], xs, va, vb)
        if not feas:
            refined.append({"pair": (va, vb), "ok": False})
            print(f"  精化 V{va}+V{vb}: 违反可观测域, 舍弃")
            continue
        la, _ = shortest_path(VEHICLE_INIT[va], sa)
        lb, _ = shortest_path(VEHICLE_INIT[vb], sb)
        ok = (la is not None and la <= BUDGET_M + 1e-9
              and lb is not None and lb <= BUDGET_M + 1e-9)
        met = metrics_full(sa, sb)
        refined.append({"pair": (va, vb), "ok": ok, "s1": sa, "s2": sb,
                        "L1": la, "L2": lb, "metrics": met,
                        "gamma": crossing_angle(sa, sb),
                        "B": math.dist(sa, sb)})
        print(f"  精化 V{va}+V{vb}: RMSE={met['rmse']:.5f} m, "
              f"L={la if la else -1:.2f}/{lb if lb else -1:.2f} m, "
              f"gamma={refined[-1]['gamma']:.2f}°, B={refined[-1]['B']:.3f} m"
              f"{'' if ok else '  [行程超预算, 舍弃]'}")

    ok_list = [r for r in refined if r.get("ok")]
    if not ok_list:
        raise RuntimeError("无满足全部约束的精化解")

    # --- 最优解族: 下界等号 iff d1=d2=r80 且 gamma=90° -> 族参数 theta(站点角),
    #     两站角 (theta, theta+90) 均落在 r80 可行弧内; 族上 RMSE 恒等于下界,
    #     故按 ①RMSE(同值) ②总行程最小 ③车号序 破平, 直接在族上做行程最优化.
    from q1_observable_domain import feasible_arcs, _circle_point
    arcs, _ = feasible_arcs(T2T, r80)

    def on_arc(a):
        return any(lo - 1e-9 <= a <= hi + 1e-9 for lo, hi in arcs)

    feas_angles = [round(float(a), 3) for a in np.arange(0.0, 360.0, 0.05)
                   if on_arc(float(a))]
    idx = {a: i for i, a in enumerate(feas_angles)}
    trav = {}
    for v in sorted(VEHICLE_INIT):
        row = []
        for a in feas_angles:
            L, _ = shortest_path(VEHICLE_INIT[v], _circle_point(T2T, r80, a))
            row.append(math.inf if L is None else float(L))
        trav[v] = row
    print(f"\n族上行程搜索: r80 可行弧格点 {len(feas_angles)} 个 (0.05° 步长)")

    fam = []
    vs = sorted(VEHICLE_INIT)
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            va_, vb_ = vs[i], vs[j]
            bf = None
            for a in feas_angles:
                a90 = round((a + 90.0) % 360.0, 3)
                if a90 not in idx:
                    continue
                ia, ib = idx[a], idx[a90]
                for (mv, ms), (nv, ns) in (((va_, ia), (vb_, ib)),
                                           ((va_, ib), (vb_, ia))):
                    La, Lb = trav[mv][ms], trav[nv][ns]
                    if La <= BUDGET_M and Lb <= BUDGET_M:
                        cand = (La + Lb, a, mv, nv)
                        if bf is None or cand < bf:
                            bf = cand
            if bf is not None:
                fam.append({"pair": (va_, vb_), "total_m": bf[0],
                            "theta_deg": bf[1], "vehA": bf[2], "vehB": bf[3]})
                print(f"  族搜索 V{va_}+V{vb_}: 总行程 {bf[0]:.2f} m @ theta={bf[1]:.2f}°")
            else:
                fam.append({"pair": (va_, vb_), "total_m": math.inf})
                print(f"  族搜索 V{va_}+V{vb_}: 族内无双双可达成员")
    fam_ok = [f for f in fam if math.isfinite(f["total_m"])]
    fam_ok.sort(key=lambda f: (f["total_m"], f["pair"]))
    bf = fam_ok[0]
    va, vb = bf["vehA"], bf["vehB"]
    th = bf["theta_deg"]
    sa = _circle_point(T2T, r80, th)
    sb = _circle_point(T2T, r80, (th + 90.0) % 360.0)
    La, pa = shortest_path(VEHICLE_INIT[va], sa)
    Lb, pb = shortest_path(VEHICLE_INIT[vb], sb)
    met = metrics_full(sa, sb)
    best = {"pair": (va, vb), "s1": sa, "s2": sb, "L1": La, "L2": Lb,
            "metrics": met, "gamma": crossing_angle(sa, sb),
            "B": math.dist(sa, sb), "theta_deg": th}

    # --- 最优性核对
    lb_rmse = SIGMA_B * math.sqrt(2) * r80
    gap = best["metrics"]["rmse"] - lb_rmse
    print(f"\n族核对: 精化阶段最优 RMSE "
          f"{min(r['metrics']['rmse'] for r in ok_list):.5f} m, "
          f"族上 RMSE {met['rmse']:.5f} m, 下界 {lb_rmse:.5f} m")

    # --- sigma 无关性: C ∝ sigma^2 -> 下界与解族(等号组)不随 sigma 变, RMSE* ∝ sigma
    sig_check = {}
    rmse_min = min(r["metrics"]["rmse"] for r in ok_list)
    tie_group = sorted({tuple(r["pair"]) for r in ok_list
                        if r["metrics"]["rmse"] <= rmse_min + 1e-4})
    for s in (0.5, 1.0, 2.0):
        sig_check[f"sigma={s}°"] = {
            "optimal_tie_group": [list(p) for p in tie_group],
            "rmse_at_chosen": round(met["rmse"] * s, 5),
            "note": "argmin 不随 sigma 变化 (C ∝ sigma²)"}

    # --- 五项输出
    va, vb = best["pair"]
    sa, sb = best["s1"], best["s2"]
    summary = {
        "sigma_b_deg": 1.0, "budget_m": BUDGET_M,
        "r80_m": round(r80, 4),
        "analytic_lower_bound_rmse_m": round(lb_rmse, 5),
        "optimality_gap_m": round(gap, 6),
        "reachability_counts": {f"V{k}": v for k, v in reach_info.items()},
        "coarse_screening": [{"pair": list(c["pair"]),
                              "rmse": round(c["rmse"], 5)} for c in coarse],
        "refined": [{"pair": list(r["pair"]), "accepted": r.get("ok", False),
                     "rmse_m": round(r["metrics"]["rmse"], 5),
                     "gamma_deg": round(r["gamma"], 2),
                     "baseline_m": round(r["B"], 3),
                     "travel_m": [round(r["L1"], 3), round(r["L2"], 3)],
                     "travel_s": [round(r["L1"] / 3.0, 3), round(r["L2"] / 3.0, 3)]}
                    for r in refined if r.get("ok")],
        "chosen": {
            "vehicle_pair": [va, vb],
            "station_A_veh": va, "station_A": [round(sa[0], 4), round(sa[1], 4)],
            "station_B_veh": vb, "station_B": [round(sb[0], 4), round(sb[1], 4)],
            "baseline_m": round(best["B"], 4),
            "crossing_angle_deg": round(best["gamma"], 3),
            "travel_time_s": [round(best["L1"] / 3.0, 3), round(best["L2"] / 3.0, 3)],
            "rmse_m": round(best["metrics"]["rmse"], 5),
            "detC_m4": best["metrics"]["detC"],
            "condC": best["metrics"]["condC"],
            "ellipse_1sigma_axes_m": [round(a, 5) for a in best["metrics"]["axes_1sigma"]],
            "elevation_deg": [round(elevation_angle(sa, TARGETS["T2"]), 3),
                              round(elevation_angle(sb, TARGETS["T2"]), 3)],
        },
        "optimal_family": {
            "note": "RMSE 下界等号解族: 两站均在 r80 圆且交角 90°, 由旋转角 theta 参数化; "
                    "族上精度恒等, 以总行程最小破平",
            "feasible_arc_deg": [[round(a0, 2), round(a1, 2)] for a0, a1 in arcs],
            "theta_grid_deg": 0.05,
            "pair_minima": [{"pair": list(f["pair"]),
                             "total_travel_m": round(f["total_m"], 3),
                             "theta_deg": round(f.get("theta_deg", float("nan")), 2)}
                            for f in fam],
            "chosen_theta_deg": round(th, 2),
        },
        "sigma_independence": sig_check,
        "mc_validation": {},
    }

    # --- MC 验证
    mc = mc_rmse(sa, sb, n=200000)
    summary["mc_validation"] = {"n": 200000, "rmse_mc_m": round(mc, 5),
                                "rmse_analytic_m": round(best["metrics"]["rmse"], 5),
                                "rel_diff_pct": round(abs(mc - best["metrics"]["rmse"])
                                                      / best["metrics"]["rmse"] * 100, 2)}
    print(f"\n=== 最优车-站对 (解族 theta={th:.2f}°, 族上总行程破平) ===")
    print(f"车对: V{va}+V{vb}")
    print(f"站点A(V{va}) = ({sa[0]:.4f}, {sa[1]:.4f}), 站点B(V{vb}) = ({sb[0]:.4f}, {sb[1]:.4f})")
    print(f"基线 B = {best['B']:.4f} m, 交角 gamma = {best['gamma']:.3f}°")
    print(f"行程: V{va} {best['L1']:.3f} m = {best['L1']/3:.2f} s; "
          f"V{vb} {best['L2']:.3f} m = {best['L2']/3:.2f} s")
    print(f"RMSE = sqrt(tr C) = {best['metrics']['rmse']:.5f} m "
          f"(解析下界 {lb_rmse:.5f} m, 差 {gap*100:.2f}%)")
    print(f"det C = {best['metrics']['detC']:.3e} m^4, cond(C) = {best['metrics']['condC']:.3f}")
    print(f"MC 验证: {summary['mc_validation']['rmse_mc_m']:.5f} m "
          f"(相对差 {summary['mc_validation']['rel_diff_pct']}%)")

    # --- 图 1: 可达性面板
    chosen_map = {va: sa, vb: sb}
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 8.6))
    for ax, vid in zip(axes.flat, sorted(VEHICLE_INIT)):
        dist = geods.get(vid)
        if dist is None:
            dist = geodesic_field(road, VEHICLE_INIT[vid], xs)
        plot_reachability(ax, xs, region, dist, vid,
                          chosen_map if vid in (va, vb) else {}, BUDGET_M)
        n_ok = reach_info[vid]
        ax.set_title(ax.get_title() + f"  [{n_ok} cells]", fontsize=9)
    axes.flat[5].axis("off")
    fig.suptitle("Q2(1): per-vehicle geodesic-reachable ∩ T2-observable stations "
                 f"(budget 90 m, grid {RES} m)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "q2_1_reachability.png", dpi=300)
    plt.close(fig)

    # --- 图 2: MC 散点 + 解析椭圆
    rng = np.random.default_rng(2026)
    p1, _, _ = bearing_terms(sa); p2, _, _ = bearing_terms(sb)
    eps = rng.normal(0.0, SIGMA_B, (200000, 2))
    f1, f2 = p1 + eps[:, 0], p2 + eps[:, 1]
    v1x, v1y = np.cos(f1), np.sin(f1)
    v2x, v2y = np.cos(f2), np.sin(f2)
    cr = v1x * v2y - v1y * v2x
    okm = np.abs(cr) > 1e-12
    t1 = ((sb[0] - sa[0]) * v2y - (sb[1] - sa[1]) * v2x) / np.where(okm, cr, 1.0)
    px, py = sa[0] + t1 * v1x, sa[1] + t1 * v1y
    err2 = (px - T2[0]) ** 2 + (py - T2[1]) ** 2
    sub = rng.choice(np.flatnonzero(okm), min(4000, int(okm.sum())), replace=False)

    fig, ax = plt.subplots(figsize=(7.4, 7.0))
    ax.plot(px[sub], py[sub], ".", color="#bcd7ec", ms=2, alpha=0.5,
            label="MC estimates (4000 shown)")
    C = SIGMA_B**2 * np.linalg.inv(
        sum(np.outer(u, u) / d**2
            for u, d in ((bearing_terms(sa)[1], math.dist(sa, T2)),
                         (bearing_terms(sb)[1], math.dist(sb, T2)))))
    ev, evec = np.linalg.eigh(C)
    ang = math.degrees(math.atan2(evec[1, 1], evec[0, 1]))
    for k, col in ((1, "crimson"), (2, "darkorange")):
        ax.add_patch(Ellipse((T2[0], T2[1]), 2 * k * math.sqrt(ev[1]),
                             2 * k * math.sqrt(ev[0]), angle=ang,
                             fill=False, color=col, lw=1.4,
                             label=f"{k}$\\sigma$ analytic ellipse"))
    ax.plot(T2[0], T2[1], "^", color="k", ms=10, label="T2 projection")
    ax.plot(*sa, "*", color="royalblue", ms=14, label=f"station A (V{va})")
    ax.plot(*sb, "*", color="seagreen", ms=14, label=f"station B (V{vb})")
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
    ax.set_title(f"Q2(1) MC validation at chosen pair: RMSE_MC="
                 f"{summary['mc_validation']['rmse_mc_m']:.4f} m vs analytic "
                 f"{summary['mc_validation']['rmse_analytic_m']:.4f} m", fontsize=10)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xlim(T2[0] - 2.0, T2[0] + 2.0); ax.set_ylim(T2[1] - 2.0, T2[1] + 2.0)
    fig.tight_layout()
    fig.savefig(OUT / "q2_1_mc_validation.png", dpi=300)
    plt.close(fig)

    # --- 图 3: gamma 退化扫线
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    plot_gamma_sweep(ax, r80, SIGMA_B, best["gamma"])
    fig.tight_layout()
    fig.savefig(OUT / "q2_1_gamma_degeneracy.png", dpi=300)
    plt.close(fig)

    (OUT / "q2_1_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n输出已写入 {OUT}")


if __name__ == "__main__":
    main()
