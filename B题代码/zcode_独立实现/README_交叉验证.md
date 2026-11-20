# ZCode 独立实现（供交叉验证）

2026-11-21 由 ZCode 会话独立完成的 B 题四问建模与求解，与仓库内既有
`q1_observable_domain.py … q4_rescheduler.py`（QoderCN 实现）**互不依赖、零文件覆盖**。
用途：按 [B题-公共模块复核与交接](../../B题分析/B题-公共模块复核与交接.md) 的标准做双实现交叉验证。

## 文件

```
swb_common.py   公共核心：三约束可行性(向量化) / 3D slab 视线遮挡 / 云台运动学 /
                建筑R_V-邻域(圆角)可视图寻路 / 不动点时间轴仿真 / 调度器
problem1.py     Q1  可观测域 + 最优站        → results/problem1_results.json, fig_p1_region.png
problem2_1.py   Q2(1) 交会精度模型 + T2 双站  → results/problem2_1_results.json, fig_p2_1_accuracy.png
problem2_2.py   Q2(2) 分配路由 + 启发式对比   → results/problem2_2_results.json, 2 张甘特图
problem3.py     Q3  误差传播 + MC 验证        → results/problem3_results.json, fig_p3_error.png
problem4.py     Q4  45 场景重调度 KEEP/RES    → results/problem4_results.json, fig_p4_comparison.png
README.md       完整建模文档（公式/算法/全部关键数字/论文写作要点）
```

复现：`python problem1.py && python problem2_1.py && python problem2_2.py && python problem3.py && python problem4.py`
（Python ≥3.10，numpy/scipy/matplotlib；总运行约 8–10 分钟，Q4 占大半）

## 关键结果（论文可直接引用）

| 项 | 值 |
|---|---|
| Q1 最优站 | (80.50, 63.65)，r*=7.5821 m（=43/tan80°），仰角 80°，Q*=0.1319 m⁻¹；平局决胜：距 O 最近→方位角最小 |
| Q2(1) 最优组合 | V3+V5，站 (6.00,75.44)+(13.96,77.52)，B=8.23 m，γ=90°，σ_RMS=0.1436 m（解析下界 √2·σ·r_min）；精确 MC CEP50=0.1197 m |
| Q2(2) max 观测开始 | **10.65 s**（V3 链 T1→T3，V1/V4 于 B2 北侧站 T2，V2/V5 原位观测 T3/T1）vs 最近assign启发式 21.09 s（−49.5%） |
| Q3 | s* = 66.9 m；80 m 约束 a+80b ≤ 1.6384 m；s≈20 m 起定位误差方差占比 >50%；解析 vs 精确 MC 最大偏差 0.02pp |
| Q4 | 45 场景（NA16+FA5+CO3+组合21）；RES 43/45 场景 η=1.0；仅 V5 失效与纯 CO 场景 RES≈KEEP（诚实无收益案例） |

## 与既有实现的关键差异（交叉验证时重点核对）

1. **"原地观测"自由度**：本实现把"站址 ∈ 可观测域"完整放开（含车辆当前位置，行驶 0 s；
   远距站如 V3 对 T3 r≈86 m 亦合法）。Q2(2) 从"就近设站"直觉的 ~21 s 压到 10.65 s 的主要来源。
2. **寻路：精确 R_V-邻域（圆角矩形）判定**。交接文档第 2 节的 B1 角点反例
   （合法点 (20.8,14.8)，净距 1.697 m，曾误判不可达）在本实现中：
   `PathPlanner.dist((20.8,14.8),(20,14)) = 1.1314 m`，与独立验算直线长一致。
   判定改为"线段到建筑矩形的最小距离 < R_V 才阻挡"（`seg_min_rect_dist`，分段二次精确解，
   与暴力采样零偏差；贴边 1.5 m 允许、进入内部拒绝）。**建议把该判定回灌 `geometry.py`。**
3. **时间轴不动点**：联合观测开始 = max(两车 ready)，车 k 段出发时间依赖伙伴 ready——
   用单调不动点迭代（≤15 轮收敛）而非一次线性推算。
4. **调度器**：分配全枚举（C(5,2)³=1000 组合×链序，max_tasks=2 剪枝）→ 行驶下界预排序 →
   前 150 精确仿真 → 前 8 坐标下降局部搜索（站址 Top-10 邻域）。
5. **Q4 commit 机制**：被打断在途/在转/等待车辆改派时的部分完成进度处理：
   保留原任务→精确沿用 planned_ready；改派→按当前插值位置+指令姿态可用。
   "故障+中断"场景的状态冻结取自**失效修正后的执行轨迹**（被打破的观测不算完成）。
6. **Q3 验证**：线性化模型 vs 精确 atan 几何 MC（4×10⁵/点）与端到端交会 MC 双重验证。

## 已知限制

* 寻路返回可视图最短路（角点节点 + 直连），是真实圆盘最短路的**上界**——与既有实现的
  同类近似一致；本轮所有用到的站址对均无不可达（已审计 INF=0）。
* 候选站为极坐标网格（2°×1 m，r≤112 m）+ 局部搜索，Q2(2)/Q4 站址精度 ~1 m 量级；
  Q2(1) 另有连续 Nelder–Mead 精化至解析下界。
* Q4 的 RE.schedule 目标为 min max start；在窗口 [0,60] 内与 (η, t_last) 字典序
  在全部 45 场景实测一致。
