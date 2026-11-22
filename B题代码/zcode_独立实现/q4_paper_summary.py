# -*- coding: utf-8 -*-
"""Summarize q4_v2_results.json into the numbers needed for the paper §Q4."""
import json
import sys

d = json.load(open(r"D:\数维杯\B题代码\zcode_独立实现\results\q4_v2_results.json",
                   encoding="utf-8"))
rows = d["scenarios"]
print("baseline:", d["baseline"]["alert_start_s"])
print("cap cost:", d["slot_cap_cost"])

def sub(prefix):
    return [r for r in rows if r["scenario"].startswith(prefix)]

def rng(vals):
    v = [x for x in vals if x is not None and x < 1e9]
    return (f"[{min(v):.2f}, {max(v):.2f}]" if v else "n/a")

# ---- (a) grid ----
na = sub("(a) NA ")
tl1 = [r["P1"]["t_last"] for r in na]
tl1c = [r["P1c"]["t_last"] for r in na]
tl2 = [r["P2"]["t_last"] for r in na]
print(f"\n(a) n={len(na)}  P1 {rng(tl1)}  P1c {rng(tl1c)}  P2 {rng(tl2)}")
p1_better = [r["scenario"] for r in na if r["P1"]["t_last"] < r["P1c"]["t_last"] - 1e-6]
p1_worse = [(r["scenario"], r["P1"]["t_last"], r["P1c"]["t_last"])
            for r in na if r["P1"]["t_last"] > r["P1c"]["t_last"] + 1e-6]
print(f"(a) P1<P1c: {len(p1_better)}  P1>P1c: {len(p1_worse)} -> {p1_worse}")
easy = [r for r in na if float(r["new_alert"][0]) in (285.0, 315.0)]
print("(a) NE sector:", [(r["scenario"].split("NA ")[1], r["P1"]["t_last"]) for r in easy])
worst_cell = max(na, key=lambda r: r["P1"]["t_last"])
print("(a) worst cell:", worst_cell["scenario"], worst_cell["P1"]["t_last"],
      worst_cell["P1"]["starts"])
nice_cell = [r for r in na if r["scenario"] == "(a) NA θ=15° R=40"][0]
print("(a) θ=15 R=40 starts:", nice_cell["P1"]["starts"])
n_best_global = sum(1 for r in rows if r["P1"]["t_last"] < r["P1c"]["t_last"] - 1e-6)
n_worse_global = sum(1 for r in rows if r["P1"]["t_last"] > r["P1c"]["t_last"] + 1e-6)
print(f"GLOBAL P1 vs P1c: better {n_best_global}, worse {n_worse_global}, "
      f"equal {86 - n_best_global - n_worse_global}")

# ---- (b) ----
fb = [r for r in rows if r["scenario"].startswith("(b) FA ") and "CO" not in r["scenario"]]
print("\n(b):", [(r["scenario"].split()[-1], r["P0"]["eta"], r["P1"]["t_last"])
                 for r in fb])

# ---- (c) ----
co = sub("(c) CO ")
print("(c):", [(r["t_rec"], r["P0"]["t_last"], r["P1"]["t_last"], r["P2"]["t_last"])
               for r in co])

# ---- (a+c) ----
ac = sub("(a+c) NA ")
t10 = [r["P1"]["t_last"] for r in ac if r["t_rec"] == 10.0]
t30 = [r["P1"]["t_last"] for r in ac if r["t_rec"] == 30.0]
print(f"\n(a+c) tau=10 P1 {rng(t10)}  tau=30 P1 {rng(t30)}")
sw = sub("(a+c sweep) ")
print("(a+c sweep) θ=15:", [(r["t_rec"], r["P1"]["t_last"], r["P2"]["t_last"])
                            for r in sorted(sw, key=lambda x: x["t_rec"])])

# ---- (b+c), (a+b), (a+b+c) ----
bc = [r for r in rows if r["scenario"].startswith("(b+c) ")]
print("\n(b+c):", [(r["scenario"].split()[-2], r["P1"]["t_last"], r["P2"]["t_last"])
                   for r in bc])
ab = sub("(a+b) NA ")
print("(a+b):", [(r["scenario"].split("+FA ")[-1], r["P0"]["eta"], r["P1"]["t_last"],
                  r["P2"]["t_last"]) for r in ab])
abc = sub("(a+b+c) ")
print("(a+b+c):", [(r["scenario"].split("FA ")[-1], r["P1"]["t_last"], r["P2"]["t_last"])
                   for r in abc])

# ---- P2 deficit ----
gaps = sorted(((r["P2"]["t_last"] - r["P1"]["t_last"], r["scenario"]) for r in rows),
              reverse=True)
print("\nP2 deficit top-5:", [(round(g, 2), s) for g, s in gaps[:5]])
p2_eta = [r["P2"]["eta"] for r in rows]
print("P2 eta min:", min(p2_eta))

# ---- sources + azimuth compliance ----
srcs = {}
for r in rows:
    s = r["P1"].get("source")
    srcs[s] = srcs.get(s, 0) + 1
print("P1 sources:", srcs)
mg = [r["P1"].get("min_az_gap_deg") for r in rows if r["P1"].get("min_az_gap_deg")
      is not None]
print(f"P1 min azimuth gap over {len(mg)} scenarios: {min(mag):.2f} deg"
      if mg else "no gap data")

# ---- sensitivity ----
print("\nsens:")
for s in d["kinematic_sensitivity"]:
    nom = s["nominal"]["P1_tlast"]
    drifts = {c: (v["P1_tlast"] - nom) / nom * 100 for c, v in s.items()
              if c not in ("scenario", "nominal", "nominal_recheck")}
    worst = max(drifts.items(), key=lambda kv: abs(kv[1]))
    print(f"  {s['scenario']}: nominal {nom:.2f}, max drift {worst[0]} "
          f"{worst[1]:+.1f}%, recheck {s['nominal_recheck']['P1_tlast']:.2f}")
