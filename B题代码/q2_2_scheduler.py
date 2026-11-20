# -*- coding: utf-8 -*-
"""Q2(2): 三警报双站联合观测的 min-max 调度(结构枚举+站点精化) + 自定义基线启发式.

运行: python -X utf8 -B q2_2_scheduler.py
输出: q2_output/ 下 q2_2_gantt.png / q2_2_map.png / q2_2_summary.json, 控制台打印方案.

口径(与框架 §4.3 一致):
- 每警报恰两辆不同车联合观测 5 s, 警报开始 = 两车 ready 较大值(GC11);
- 同一警报两站方位角须不同: 实现取 gamma=|wrap180(pan1-pan2)| >= 1° 守门, 报告实际值;
- 站点 ∈ 各警报可观测域(q1 region_mask 0.5 m 网格采样 + r80 可行弧 1° 点);
- 行驶 = shortest_path 绕障最短路 / SPEED(有限图上界性质见交接文档 §5);
- 结构空间 = 五车全用、恰一车复用: 选复用车 5 × 其有序警报对 3x2 × 其余四车分槽 4x3
  (第三警报两车可互换去重) = 360 种; 题面未要求五车全用, 此为受限搜索(论文须声明);
- 求解: 360 结构贪心选站粗筛 -> 前 12 结构 best-response 精化(复用车两槽联合 + 单槽);
  结果为可行上界, 另报离散化独立下界(每警报不计车辆共享的最早可能开始, 任意可行调度
  的 max start >= 该下界), gap 据实报告;
- 基线启发式(题目要求自定义): 逐警报独立 minimax 选站 + T1->T2->T3 先到先得指派 +
  复用车按边际延迟定序(两顺序仿真取优), 与优化方案同仿真器对比.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from geometry import (BUILDINGS, FIELD_SIZE, SPEED, TARGETS, VEHICLE_INIT,
                      elev_min_radius, shortest_path, is_station_feasible)
from simulator import GIMBAL_RATE, Task, simulate
from q1_observable_domain import _circle_point, feasible_arcs, region_mask

OUT = Path(__file__).with_name("q2_output")
OUT.mkdir(exist_ok=True)
ALERTS = ["T1", "T2", "T3"]
AZ_MIN = 1.0          # 同警报两站方位角差守门阈值(度)
N_STRUCT_TOP = 12     # 精化结构数
CAND_CAP = 600        # 每警报区域网格候选上限
TOP_STATIC = 60       # 复用车第二槽静态候选带宽


def wrap180(a):
    a = (a + 180.0) % 360.0 - 180.0
    return 180.0 if a == -180.0 else a


# ---------------------------------------------------------------- 候选站生成
def build_alert_data(name):
    T = TARGETS[name]
    xs, m = region_mask(T, 0.5)
    idx = np.argwhere(m)
    if len(idx) > CAND_CAP:
        stride = max(1, int(round(math.sqrt(len(idx) / CAND_CAP))))
        idx = idx[(idx[:, 0] % stride == 0) & (idx[:, 1] % stride == 0)]
    pts = [[float(xs[j]), float(xs[i])] for i, j in idx]
    r80 = elev_min_radius(T)
    arcs, _ = feasible_arcs(T, r80)
    for lo, hi in arcs:
        n = max(2, int(round(hi - lo)))
        for k in range(n + 1):
            p = _circle_point(T, r80, lo + (hi - lo) * k / n)
            pts.append([float(p[0]), float(p[1])])
    seen, cands = set(), []
    for p in pts:
        key = (round(p[0] * 4) / 4, round(p[1] * 4) / 4)
        if key in seen:
            continue
        seen.add(key)
        if is_station_feasible(tuple(p), T):
            cands.append(p)
    arr = np.array(cands)
    dx, dy = T.x - arr[:, 0], T.y - arr[:, 1]
    pan = np.degrees(np.arctan2(dy, dx)) % 360.0
    tilt = np.degrees(np.arctan2(T.H - 2.0, np.hypot(dx, dy)))
    slew = np.maximum(np.minimum(pan, 360.0 - pan), tilt) / GIMBAL_RATE
    return {"name": name, "T": T, "arr": arr, "pan": pan, "tilt": tilt,
            "slew": slew, "n": len(cands)}


def travel_from(pos, arr):
    out = np.full(len(arr), np.inf)
    for i, s in enumerate(arr):
        r = shortest_path(tuple(map(float, pos)), (float(s[0]), float(s[1])))
        if r is not None:
            out[i] = r[0] / SPEED
    return out


def travel_scalar(pos, s):
    r = shortest_path(tuple(map(float, pos)), (float(s[0]), float(s[1])))
    return None if r is None else r[0] / SPEED


# ---------------------------------------------------------------- 结构枚举
def gen_structures():
    out = []
    for k in range(1, 6):
        for a1 in ALERTS:
            for a2 in ALERTS:
                if a2 == a1:
                    continue
                others = [v for v in range(1, 6) if v != k]
                for p1 in others:
                    for p2 in others:
                        if p2 == p1:
                            continue
                        rest = [v for v in others if v not in (p1, p2)]
                        out.append({"k": k, "a1": a1, "a2": a2,
                                    "p1": p1, "p2": p2, "p3": rest[0], "p4": rest[1]})
    assert len(out) == 360, len(out)
    return out


def slots_of(st):
    a3 = [a for a in ALERTS if a not in (st["a1"], st["a2"])][0]
    return [(st["k"], st["a1"]), (st["p1"], st["a1"]),
            (st["k"], st["a2"]), (st["p2"], st["a2"]),
            (st["p3"], a3), (st["p4"], a3)], a3


def tasks_from(st, plan):
    a3 = [a for a in ALERTS if a not in (st["a1"], st["a2"])][0]
    return [Task(st["k"], st["a1"], plan[(st["k"], st["a1"])]),
            Task(st["p1"], st["a1"], plan[(st["p1"], st["a1"])]),
            Task(st["k"], st["a2"], plan[(st["k"], st["a2"])]),
            Task(st["p2"], st["a2"], plan[(st["p2"], st["a2"])]),
            Task(st["p3"], a3, plan[(st["p3"], a3)]),
            Task(st["p4"], a3, plan[(st["p4"], a3)])]


def alert_gamma(ad, s_a, s_b):
    p1 = math.degrees(math.atan2(ad["T"].y - s_a[1], ad["T"].x - s_a[0])) % 360.0
    p2 = math.degrees(math.atan2(ad["T"].y - s_b[1], ad["T"].x - s_b[0])) % 360.0
    return abs(wrap180(p1 - p2))


# ---------------------------------------------------------------- 选站与修方位
def dyn_order(ad, pos, idxs, k_best):
    """候选子集按 (绕障行程+slew) 升序取前 k_best 个索引."""
    best, tr = [], {}
    for i in idxs:
        t = travel_scalar(pos, ad["arr"][i])
        if t is not None:
            tr[i] = t + ad["slew"][i]
    for i in sorted(tr, key=tr.get)[:k_best]:
        best.append(i)
    return best, tr


def greedy_plan(st, AD, ready, order):
    plan = {}
    slot_list, _ = slots_of(st)
    for v, a in slot_list:
        plan[(v, a)] = tuple(map(float, AD[a]["arr"][order[(v, a)][0]]))
    s1 = plan[(st["k"], st["a1"])]
    a2 = st["a2"]
    idxs = list(order[(st["k"], a2)][:TOP_STATIC])
    tr = {}
    for i in idxs:
        t = travel_scalar(s1, AD[a2]["arr"][i])
        if t is not None:
            tr[i] = t + AD[a2]["slew"][i]
    if tr:
        i2 = min(tr, key=tr.get)
        plan[(st["k"], a2)] = tuple(map(float, AD[a2]["arr"][i2]))
    repair_all(st, plan, AD, ready, order)
    return plan


def slot_ready_proxy(v, a, st, plan, AD, ready, order):
    """修复排序用的 ready 代理: 复用车第二槽按当前第一站动态算, 其余用静态表."""
    if v == st["k"] and a == st["a2"]:
        t = travel_scalar(plan[(st["k"], st["a1"])], plan[(v, a)])
        if t is None:
            return float("inf")
        return t + AD[a]["slew"][
            int(np.argmin(np.linalg.norm(AD[a]["arr"] - plan[(v, a)], axis=1)))]
    i = int(np.argmin(np.linalg.norm(AD[a]["arr"] - plan[(v, a)], axis=1)))
    return float(ready[(v, a)][i])


def repair_all(st, plan, AD, ready, order):
    """逐警报检查方位角差; 不足 AZ_MIN 时置换 ready 代理较小的可换槽."""
    for a in ALERTS:
        vs = [v for v, aa in plan if aa == a]
        if len(vs) < 2:
            continue
        s_a, s_b = plan[(vs[0], a)], plan[(vs[1], a)]
        if alert_gamma(AD[a], s_a, s_b) >= AZ_MIN:
            continue
        best = None
        for v_slot in vs:
            other = tuple(map(float, plan[(vs[0], a) if v_slot == vs[1] else (vs[1], a)]))
            idxs = list(order[(v_slot, a)][:40])
            for i in idxs:
                s = tuple(map(float, AD[a]["arr"][i]))
                if alert_gamma(AD[a], s, other) < AZ_MIN:
                    continue
                old = plan[(v_slot, a)]
                plan[(v_slot, a)] = s
                r_new = slot_ready_proxy(v_slot, a, st, plan, AD, ready, order)
                plan[(v_slot, a)] = old
                r_other = slot_ready_proxy(
                    vs[0] if v_slot == vs[1] else vs[1], a, st, plan, AD, ready, order)
                key = (max(r_new, r_other), r_new)
                if best is None or key < best[0]:
                    best = (key, v_slot, s)
        if best is not None:
            plan[(best[1], a)] = best[2]


# ---------------------------------------------------------------- 精化
def refine(st, plan, AD, ready, order, rounds=3):
    plan = {kk: tuple(map(float, vv)) for kk, vv in plan.items()}
    tasks = tasks_from(st, plan)
    best_res = simulate(tasks)
    best_ms = best_res.max_start

    def try_plan(p):
        nonlocal best_res, best_ms
        for a in ALERTS:
            vs = [v for v, aa in p if aa == a]
            if len(vs) == 2 and alert_gamma(AD[a], p[(vs[0], a)], p[(vs[1], a)]) < AZ_MIN:
                return False
        try:
            r = simulate(tasks_from(st, p))
        except Exception:
            return False
        if r.max_start < best_ms - 1e-6:
            best_ms, best_res = r.max_start, r
            return True
        return False

    def cand_union(v, a, prev_pos, k_dyn=25):
        arr = AD[a]["arr"]
        static = list(order[(v, a)][:100])
        d = np.linalg.norm(arr - np.asarray(prev_pos), axis=1)
        near = list(np.argsort(d)[:100])
        idxs, seen = [], set()
        for i in static + near:
            if i not in seen:
                seen.add(i)
                idxs.append(int(i))
        scored, tr = [], travel_from(prev_pos, arr[idxs])
        for j, i in enumerate(idxs):
            if np.isfinite(tr[j]):
                scored.append((tr[j] + AD[a]["slew"][i], i))
        scored.sort()
        cur = plan[(v, a)]
        ci = int(np.argmin(np.linalg.norm(arr - cur, axis=1)))
        return [i for _, i in scored[:k_dyn]] + ([ci] if ci not in
                                                 [i for _, i in scored[:k_dyn]] else [])

    slot_list, _ = slots_of(st)
    for rd in range(rounds):
        improved_any = False
        if rd == 0:   # 复用车两槽联合 best-response
            k, a1, a2 = st["k"], st["a1"], st["a2"]
            arr1 = AD[a1]["arr"]
            d1 = np.linalg.norm(arr1 - np.asarray(VEHICLE_INIT[k]), axis=1)
            idx1 = [int(i) for i in np.argsort(d1)[:100]]
            s1_scored = []
            tr1 = travel_from(VEHICLE_INIT[k], arr1[idx1])
            for j, i in enumerate(idx1):
                if np.isfinite(tr1[j]):
                    s1_scored.append((tr1[j] + AD[a1]["slew"][i], i))
            s1_scored.sort()
            for _, i1 in s1_scored[:15]:
                s1 = tuple(map(float, arr1[i1]))
                idx2 = cand_union(k, a2, s1, k_dyn=8)
                for i2 in idx2:
                    p2 = dict(plan)
                    p2[(k, a1)] = s1
                    p2[(k, a2)] = tuple(map(float, AD[a2]["arr"][i2]))
                    if try_plan(p2):
                        plan = p2
                        improved_any = True
        for v, a in slot_list:
            prev = VEHICLE_INIT[v]
            if v == st["k"] and a == st["a2"]:
                prev = plan[(st["k"], st["a1"])]
            for i in cand_union(v, a, prev):
                p2 = dict(plan)
                p2[(v, a)] = tuple(map(float, AD[a]["arr"][i]))
                if try_plan(p2):
                    plan = p2
                    improved_any = True
        if not improved_any:
            break
    return plan, best_res, best_ms


# ---------------------------------------------------------------- 基线启发式
def baseline_plan(AD, ready, order):
    plan, remaining, assigned = {}, set(range(1, 6)), {}
    forced = None
    for a in ALERTS:
        avail = sorted(remaining)
        if len(avail) >= 2:
            scored = sorted(avail, key=lambda v: float(ready[(v, a)][order[(v, a)][0]]))
            v1, v2 = scored[0], scored[1]
            remaining -= {v1, v2}
            assigned[v1] = a
            assigned[v2] = a
            plan[(v1, a)] = tuple(map(float, AD[a]["arr"][order[(v1, a)][0]]))
            plan[(v2, a)] = tuple(map(float, AD[a]["arr"][order[(v2, a)][0]]))
        else:
            forced = a
    assert forced is not None
    v_last = sorted(remaining)[0]
    best = None
    for v_d, a_first in assigned.items():
        s_first = plan[(v_d, a_first)]
        idxs = cand_top(AD[forced], s_first, 100)[:25]
        for i in idxs:
            t = travel_scalar(s_first, AD[forced]["arr"][i])
            if t is None:
                continue
            key = t + AD[forced]["slew"][i]
            if best is None or key < best[0]:
                best = (key, v_d, tuple(map(float, AD[forced]["arr"][i])))
    _, v_d, s_d = best
    plan[(v_d, forced)] = s_d
    plan[(v_last, forced)] = tuple(map(float, AD[forced]["arr"][order[(v_last, forced)][0]]))
    a_first = assigned[v_d]
    repair_all({"k": v_d, "a1": a_first, "a2": forced}, plan, AD, ready, order)

    def build(order_first):
        tasks = []
        for a in ALERTS:
            vs = [v for v, aa in plan if aa == a]
            for v in vs:
                tasks.append(Task(v, a, plan[(v, a)]))
        # 复用车两任务顺序重排
        td = [t for t in tasks if t.vehicle == v_d]
        assert len(td) == 2
        rest = [t for t in tasks if t.vehicle != v_d]
        if order_first == forced:
            td.sort(key=lambda t: 0 if t.alert == forced else 1)
        else:
            td.sort(key=lambda t: 0 if t.alert != forced else 1)
        return rest + td

    cands = []
    for of in (forced, a_first):
        try:
            r = simulate(build(of))
            cands.append((r.max_start, of, r))
        except Exception:
            continue
    cands.sort(key=lambda c: c[0])
    ms, of, res = cands[0]
    return plan, v_d, a_first, of, res


def cand_top(ad, pos, cap):
    d = np.linalg.norm(ad["arr"] - np.asarray(pos), axis=1)
    return [int(i) for i in np.argsort(d)[:cap]]


# ---------------------------------------------------------------- 下界
def lower_bound(AD, ready, order):
    per = {}
    for a in ALERTS:
        best_pair = float("inf")
        for i1 in range(5):
            for i2 in range(i1 + 1, 5):
                v1, v2 = i1 + 1, i2 + 1
                j1 = order[(v1, a)][0]
                s1 = AD[a]["arr"][j1]
                r1 = float(ready[(v1, a)][j1])
                r2_best = float("inf")
                for j2 in order[(v2, a)][:40]:
                    s2 = AD[a]["arr"][j2]
                    if alert_gamma(AD[a], tuple(map(float, s1)),
                                   tuple(map(float, s2))) >= AZ_MIN:
                        r2_best = float(ready[(v2, a)][j2])
                        break
                best_pair = min(best_pair, max(r1, r2_best))
        per[a] = best_pair
    return max(per.values()), per


# ---------------------------------------------------------------- 图件
KIND_COLOR = {"arrive": "#5dade2", "slew": "#f5b041", "wait": "#bdc3c7",
              "obs": "#58d68d", "reset": "#af7ac5"}


def plot_gantt(res, title, fn):
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    vids = sorted(res.events)
    for row, v in enumerate(vids):
        for e in res.events[v]:
            if e["kind"] == "depart":
                continue
            ax.barh(row, e["t_end"] - e["t_start"], left=e["t_start"], height=0.55,
                    color=KIND_COLOR.get(e["kind"], "#cccccc"), edgecolor="k",
                    linewidth=0.4)
            if e["kind"] == "obs":
                ax.text((e["t_start"] + e["t_end"]) / 2, row, e["alert"],
                        ha="center", va="center", fontsize=8, fontweight="bold")
    starts = sorted(res.alert_start.items(), key=lambda kv: kv[1])
    for i, (a, s) in enumerate(starts):
        ax.axvline(s, color="crimson", ls="--", lw=1.0)
        nxt = starts[i + 1][1] if i + 1 < len(starts) else 1e9
        if nxt - s < 0.8:
            ax.text(s - 0.08, -0.62, f"{a} {s:.2f}s", color="crimson",
                    fontsize=8, ha="right", va="center")
        else:
            ax.text(s + 0.08, -0.62, f"{a} {s:.2f}s", color="crimson",
                    fontsize=8, ha="left", va="center")
    ax.set_ylim(-1.1, len(vids) - 0.5)
    ax.set_yticks(range(len(vids)))
    ax.set_yticklabels([f"V{v}" for v in vids])
    ax.set_xlabel("time (s)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in KIND_COLOR.values()]
    ax.legend(handles, list(KIND_COLOR.keys()), fontsize=8, ncol=5,
              loc="upper right")
    ax.set_title(f"{title}  max start = {res.max_start:.3f} s, "
                 f"makespan = {res.makespan:.3f} s", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / fn, dpi=300)
    plt.close(fig)


def plot_map(st, plan, fn):
    fig, ax = plt.subplots(figsize=(7.6, 7.4))
    for bi, b in enumerate(BUILDINGS, start=1):
        ax.add_patch(Rectangle((b.x0, b.y0), b.x1 - b.x0, b.y1 - b.y0,
                               facecolor="#9e9e9e", edgecolor="k", zorder=3))
        ax.text((b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2, f"B{bi}", ha="center",
                va="center", fontsize=8, zorder=4)
    colors = {"T1": "crimson", "T2": "purple", "T3": "seagreen"}
    for a in ALERTS:
        T = TARGETS[a]
        ax.plot(T.x, T.y, "^", color=colors[a], ms=10, zorder=6)
        ax.annotate(f"{a} ({T.x:.1f},{T.y:.1f},{T.H:.0f})", (T.x, T.y),
                    textcoords="offset points", xytext=(5, 5), fontsize=8, zorder=6)
        for (v, aa), s in plan.items():
            if aa != a:
                continue
            ax.plot(*s, "*", color=colors[a], ms=13, zorder=7)
            ax.annotate(f"V{v}", s, textcoords="offset points", xytext=(4, -10),
                        fontsize=8, color=colors[a], zorder=7)
    for v, p in VEHICLE_INIT.items():
        ax.plot(*p, "s", color="darkorange", ms=6, zorder=5)
        ax.annotate(f"V{v}", p, textcoords="offset points", xytext=(4, 4),
                    fontsize=8, zorder=5)
    for (v, a), s in plan.items():
        src = VEHICLE_INIT[v]
        if v == st["k"] and a == st["a2"]:
            src = plan[(st["k"], st["a1"])]
        r = shortest_path(tuple(map(float, src)), tuple(map(float, s)))
        if r is not None:
            pts = np.array(r[1])
            ax.plot(pts[:, 0], pts[:, 1], "-", color="#7f8c8d", lw=0.9,
                    alpha=0.8, zorder=2)
    ax.set_xlim(0, FIELD_SIZE); ax.set_ylim(0, FIELD_SIZE)
    ax.set_aspect("equal"); ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
    ax.set_title("Q2(2) optimal assignment: stations, routes and doubling vehicle",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / fn, dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------- route_audit 互证
def audit_route(pts):
    fn = OUT / "audit_tmp.json"
    fn.write_text(json.dumps([[float(p[0]), float(p[1])] for p in pts]),
                  encoding="utf-8")
    r = subprocess.run([sys.executable, "-X", "utf8", "-B", "route_audit.py", str(fn)],
                       capture_output=True, text=True,
                       cwd=Path(__file__).parent)
    try:
        return json.loads(r.stdout)
    except Exception:
        return {"error": (r.stderr or r.stdout)[:300]}


# ---------------------------------------------------------------- main
def main():
    print("== 生成候选站与静态代价表 ==")
    AD = {a: build_alert_data(a) for a in ALERTS}
    for a in ALERTS:
        print(f"  {a}: 候选站 {AD[a]['n']} 个")
    ready, order = {}, {}
    for v in range(1, 6):
        for a in ALERTS:
            tr = travel_from(VEHICLE_INIT[v], AD[a]["arr"])
            ready[(v, a)] = tr + AD[a]["slew"]
            order[(v, a)] = np.argsort(ready[(v, a)])

    print("== 阶段 A: 360 结构贪心粗筛 ==")
    structs = gen_structures()
    stageA = []
    for n, st in enumerate(structs):
        try:
            plan = greedy_plan(st, AD, ready, order)
            res = simulate(tasks_from(st, plan))
            stageA.append({"st": st, "plan": plan, "res": res})
        except Exception:
            stageA.append({"st": st, "plan": None, "res": None})
        if (n + 1) % 60 == 0:
            fin = [x for x in stageA if x["res"] is not None]
            print(f"  {n+1}/360, 可行 {len(fin)}, 当前最优 "
                  f"{min(x['res'].max_start for x in fin):.3f} s")
    finA = [x for x in stageA if x["res"] is not None]
    finA.sort(key=lambda x: x["res"].max_start)
    print(f"  粗筛完成: 可行 {len(finA)}/360, 最优 {finA[0]['res'].max_start:.3f} s")

    print(f"== 阶段 B: 前 {N_STRUCT_TOP} 结构 best-response 精化 ==")
    best = None
    for x in finA[:N_STRUCT_TOP]:
        plan, res, ms = refine(x["st"], x["plan"], AD, ready, order)
        print(f"  V{x['st']['k']}复用 {x['st']['a1']}->{x['st']['a2']}: "
              f"{x['res'].max_start:.3f} -> {ms:.3f} s")
        if best is None or ms < best[0]:
            best = (ms, x["st"], plan, res)
    ms, st, plan, res = best

    # --- 合规清单
    compliance = {"stations_feasible": True, "two_vehicles_per_alert": True,
                  "azimuth_gaps_deg": {}}
    for a in ALERTS:
        vs = [v for v, aa in plan if aa == a]
        compliance["two_vehicles_per_alert"] &= (len(vs) == 2 and vs[0] != vs[1])
        for v in vs:
            compliance["stations_feasible"] &= bool(
                is_station_feasible(plan[(v, a)], TARGETS[a]))
        compliance["azimuth_gaps_deg"][a] = round(
            alert_gamma(AD[a], plan[(vs[0], a)], plan[(vs[1], a)]), 2)

    a3 = [a for a in ALERTS if a not in (st["a1"], st["a2"])][0]
    print("\n=== Q2(2) 最优方案 ===")
    print(f"复用车 V{st['k']}: {st['a1']} -> {st['a2']}; "
          f"{st['a1']}: V{st['k']}+V{st['p1']}, {st['a2']}: V{st['k']}+V{st['p2']}, "
          f"{a3}: V{st['p3']}+V{st['p4']}")
    for a in ALERTS:
        vs = [v for v, aa in plan if aa == a]
        ss = [plan[(v, a)] for v in vs]
        print(f"  {a}: V{vs[0]}@({ss[0][0]:.3f},{ss[0][1]:.3f}) + "
              f"V{vs[1]}@({ss[1][0]:.3f},{ss[1][1]:.3f}), "
              f"start={res.alert_start[a]:.3f}s, gamma={compliance['azimuth_gaps_deg'][a]:.1f}°")
    print(f"max start = {res.max_start:.3f} s, makespan = {res.makespan:.3f} s")

    # --- route_audit 互证(复用车两段)
    src2 = plan[(st["k"], st["a1"])]
    legs = {}
    for tag, (src, dst) in (("leg1", (VEHICLE_INIT[st["k"]], plan[(st["k"], st["a1"])])),
                            ("leg2", (src2, plan[(st["k"], st["a2"])]))):
        r = shortest_path(tuple(map(float, src)), tuple(map(float, dst)))
        legs[tag] = {"sim_seconds": r[0] / SPEED if r else None,
                     "audit": audit_route(r[1]) if r else None}

    # --- 基线
    print("== 基线启发式 ==")
    bplan, v_d, a_first, of, bres = baseline_plan(AD, ready, order)
    print(f"  基线 max start = {bres.max_start:.3f} s (复用车 V{v_d}, 先执行 {of})")

    # --- 下界与 gap
    lb, lb_per = lower_bound(AD, ready, order)
    print(f"  离散化独立下界 = {lb:.3f} s, gap = {res.max_start - lb:.3f} s")

    # --- 灵敏度 ±10%
    import simulator as sim_mod
    sens = []
    for tag, (g, sp) in {"gimbal-10%": (GIMBAL_RATE * 0.9, SPEED),
                         "gimbal+10%": (GIMBAL_RATE * 1.1, SPEED),
                         "speed-10%": (GIMBAL_RATE, SPEED * 0.9),
                         "speed+10%": (GIMBAL_RATE, SPEED * 1.1)}.items():
        og, os_ = sim_mod.GIMBAL_RATE, sim_mod.SPEED
        try:
            sim_mod.GIMBAL_RATE, sim_mod.SPEED = g, sp
            r = simulate(tasks_from(st, plan))
            sens.append({"case": tag, "max_start_s": round(r.max_start, 3),
                         "makespan_s": round(r.makespan, 3)})
        finally:
            sim_mod.GIMBAL_RATE, sim_mod.SPEED = og, os_
    print(f"  灵敏度: {[(s['case'], s['max_start_s']) for s in sens]}")

    # --- 时间轴与利用率
    util = {}
    for v in sorted(res.events):
        busy = sum(e["t_end"] - e["t_start"] for e in res.events[v]
                   if e["kind"] in ("arrive", "slew", "obs", "reset"))
        util[f"V{v}"] = {"busy_s": round(busy, 3),
                         "idle_s": round(res.makespan - busy, 3)}

    summary = {
        "objective": "min max_a start_a (GC11: 警报开始=两车ready较大值, 联合观测5s)",
        "structure_space": "5x3x2x4!/2=360 (五车全用、恰一复用, 受限搜索)",
        "candidates_per_alert": {a: AD[a]["n"] for a in ALERTS},
        "azimuth_guard_deg": AZ_MIN,
        "stageA_feasible": len(finA),
        "stageA_top10": [{"double_veh": x["st"]["k"],
                          "ordered_alerts": [x["st"]["a1"], x["st"]["a2"]],
                          "max_start_s": round(x["res"].max_start, 3)}
                         for x in finA[:10]],
        "best": {
            "structure": {"double_veh": st["k"],
                          "ordered_alerts": [st["a1"], st["a2"]],
                          "a1_partner": st["p1"], "a2_partner": st["p2"],
                          "a3_vehicles": [st["p3"], st["p4"]]},
            "stations": {f"V{v}@{a}": [round(s[0], 4), round(s[1], 4)]
                         for (v, a), s in plan.items()},
            "alert_start_s": {a: round(res.alert_start[a], 3) for a in ALERTS},
            "alert_end_s": {a: round(res.alert_end[a], 3) for a in ALERTS},
            "azimuth_gap_deg": compliance["azimuth_gaps_deg"],
            "max_start_s": round(res.max_start, 3),
            "makespan_s": round(res.makespan, 3),
            "vehicle_utilization": util,
            "timelines": {f"V{v}": [{ "kind": e["kind"], "alert": e["alert"],
                                      "t_start": round(e["t_start"], 3),
                                      "t_end": round(e["t_end"], 3),
                                      "detail": e["detail"]}
                                    for e in res.events[v]]
                          for v in sorted(res.events)},
        },
        "doubling_vehicle_route_audit": legs,
        "baseline": {
            "definition": "逐警报独立minimax选站 + T1->T2->T3先到先得指派 + 复用车按边际延迟定序(两顺序取优)",
            "double_veh": v_d, "first_alert_chosen": of,
            "stations": {f"V{v}@{a}": [round(s[0], 4), round(s[1], 4)]
                         for (v, a), s in bplan.items()},
            "alert_start_s": {a: round(bres.alert_start[a], 3) for a in ALERTS},
            "max_start_s": round(bres.max_start, 3),
            "makespan_s": round(bres.makespan, 3),
        },
        "comparison": {
            "max_start_delta_s": round(bres.max_start - res.max_start, 3),
            "per_alert_delta_s": {a: round(bres.alert_start[a] - res.alert_start[a], 3)
                                  for a in ALERTS},
        },
        "lower_bound": {"per_alert_s": {a: round(v, 3) for a, v in lb_per.items()},
                        "max_lb_s": round(lb, 3),
                        "note": "每警报不计车辆共享的最早可能开始(离散候选集上取min), "
                                "任意可行调度的 max start >= 该值"},
        "gap_s": round(res.max_start - lb, 3),
        "sensitivity": sens,
        "compliance": {"stations_feasible": bool(compliance["stations_feasible"]),
                       "two_vehicles_per_alert": bool(compliance["two_vehicles_per_alert"])},
    }
    (OUT / "q2_2_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    plot_gantt(res, "Q2(2) optimal schedule", "q2_2_gantt.png")
    plot_gantt(bres, "Q2(2) baseline heuristic", "q2_2_gantt_baseline.png")
    plot_map(st, plan, "q2_2_map.png")
    print(f"\n输出已写入 {OUT}")


if __name__ == "__main__":
    main()
