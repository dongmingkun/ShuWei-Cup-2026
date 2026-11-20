# -*- coding: utf-8 -*-
"""端到端冒烟演示: 3警报/6任务槽/5车(V2复用) 的一套完整方案跑通仿真器."""
import math

from geometry import TARGETS, elev_min_radius
from simulator import Task, simulate

t1, t2, t3 = TARGETS["T1"], TARGETS["T2"], TARGETS["T3"]
r1, r2, r3 = elev_min_radius(t1), elev_min_radius(t2), elev_min_radius(t3)

plan = [
    Task(1, "T3", (t3.x, t3.y - r3)),          # V1: T3 南站(直连)
    Task(2, "T3", (t3.x, t3.y + r3)),          # V2: T3 北站, 完成后再赴 T1(复用)
    Task(2, "T1", (t1.x, t1.y - r1)),
    Task(3, "T2", (t2.x - r2, t2.y)),          # V3: T2 西站(T2 投影在 B3 上方, 西侧可行)
    Task(4, "T1", (t1.x + r1, t1.y)),          # V4: T1 东站(直连)
    Task(5, "T2", (t2.x - r2 * math.sqrt(0.5), t2.y + r2 * math.sqrt(0.5))),  # V5: T2 西北站
]

if __name__ == "__main__":
    res = simulate(plan)
    for v in sorted(res.events):
        print(res.timeline_str(v))
        print()
    print("警报开始时刻:", {a: round(s, 3) for a, s in res.alert_start.items()})
    print("警报结束时刻:", {a: round(s, 3) for a, s in res.alert_end.items()})
    print(f"max 观测开始时间(目标值) = {res.max_start:.3f} s")
    print(f"makespan = {res.makespan:.3f} s")
