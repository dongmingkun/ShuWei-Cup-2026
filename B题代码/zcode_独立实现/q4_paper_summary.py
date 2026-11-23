# -*- coding: utf-8 -*-
"""Summarize q4_v2_results.json (sigma<2m convention) into paper numbers."""
import json
import os
import sys

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "results", "q4_v2_results.json")
d = json.load(open(PATH, encoding="utf-8"))
rows = d["scenarios"]


def sub(prefix):
    return [r for r in rows if r["scenario"].startswith(prefix)]


def rng(vals):
    v = [x for x in vals if x is not None and x < 1e9]
    return f"[{min(v):.2f}, {max(v):.2f}]" if v else "n/a"


base = d["baseline"]["alert_start_s"]
print("baseline starts:", base, " max_sigma:", d["baseline"]["max_sigma_m"])
print("cap cost:", d["slot_cap_cost"])

na = sub("(a) NA ")
print(f"\n(a) n={len(na)}  P1 {rng([r['P1']['t_last'] for r in na])}  "
      f"P1c {rng([r['P1c']['t_last'] for r in na])}  "
      f"P2 {rng([r['P2']['t_last'] for r in na])}")
easies = sorted(na, key=lambda r: r["P1"]["t_last"])[:4]
hards = sorted(na, key=lambda r: -r["P1"]["t_last"])[:4]
print("  easiest:", [(r["scenario"].split("NA ")[1], round(r["P1"]["t_last"], 2),
                      r["P1"]["starts"]) for r in easies])
print("  hardest:", [(r["scenario"].split("NA ")[1], round(r["P1"]["t_last"], 2),
                      r["P1"]["starts"]) for r in hards])
rep = next(r for r in na if r["scenario"] == "(a) NA θ=15° R=40")
print("  rep cell θ=15 R=40:", round(rep["P1"]["t_last"], 2), rep["P1"]["starts"])

fb = [r for r in rows if r["scenario"].startswith("(b) FA ") and "CO" not in r
      ["scenario"]]
print("\n(b):", [(r["scenario"].split()[-1], round(r["P0"]["eta"], 2),
                  round(r["P0"]["t_last"], 2), round(r["P1"]["t_last"], 2),
                  round(r["P2"]["t_last"], 2)) for r in fb])

co = sub("(c) CO ")
print("\n(c):", [(r["t_rec"], round(r["P0"]["t_last"], 2),
                  round(r["P1"]["t_last"], 2), round(r["P2"]["t_last"], 2))
                 for r in co])

ac = sub("(a+c) NA ")
print(f"\n(a+c) τ=10 P1 {rng([r['P1']['t_last'] for r in ac if r['t_rec'] == 10.0])}"
      f"  τ=30 {rng([r['P1']['t_last'] for r in ac if r['t_rec'] == 30.0])}")
sw = sub("(a+c sweep) ")
print("(a+c sweep) θ=15:", [(r["t_rec"], round(r["P1"]["t_last"], 2),
                             round(r["P2"]["t_last"], 2))
                            for r in sorted(sw, key=lambda x: x["t_rec"])])

bc = [r for r in rows if r["scenario"].startswith("(b+c) ")]
print("\n(b+c):", [(r["scenario"].split()[-2], round(r["P0"]["eta"], 2),
                    round(r["P1"]["t_last"], 2), round(r["P2"]["t_last"], 2))
                   for r in bc])
ab = sub("(a+b) NA ")
print("(a+b):", [(r["scenario"].split("+FA ")[-1], round(r["P0"]["eta"], 2),
                  round(r["P1"]["t_last"], 2), round(r["P2"]["t_last"], 2))
                 for r in ab])
abc = sub("(a+b+c) ")
print("(a+b+c):", [(r["scenario"].split("FA ")[-1], round(r["P1"]["t_last"], 2),
                    round(r["P2"]["t_last"], 2)) for r in abc])

gaps = sorted(((r["P2"]["t_last"] - r["P1"]["t_last"], r["scenario"])
               for r in rows), reverse=True)
print("\nP2 deficit top-5:", [(round(g, 2), s) for g, s in gaps[:5]])
print("P2 eta min:", min(r["P2"]["eta"] for r in rows),
      " P1 eta min:", min(r["P1"]["eta"] for r in rows))

mg = [r["P1"]["min_gamma_deg"] for r in rows if r["P1"].get("min_gamma_deg")]
ms = [r["P1"]["max_sigma_m"] for r in rows if r["P1"].get("max_sigma_m")]
print(f"\nP1 audits over {len(mg)} plans: min gamma {min(mg):.1f} deg, "
      f"max sigma {max(ms):.3f} m")
mg2 = [r["P2"]["min_gamma_deg"] for r in rows if r["P2"].get("min_gamma_deg")]
ms2 = [r["P2"]["max_sigma_m"] for r in rows if r["P2"].get("max_sigma_m")]
print(f"P2 audits over {len(mg2)} plans: min gamma {min(mg2):.1f} deg, "
      f"max sigma {max(ms2):.3f} m")

print("\nsens:")
for s in d["kinematic_sensitivity"]:
    nom = s["nominal"]["P1_tlast"]
    dr = {c: (v["P1_tlast"] - nom) / nom * 100 for c, v in s.items()
          if c not in ("scenario", "nominal", "nominal_recheck")}
    w = max(dr.items(), key=lambda kv: abs(kv[1]))
    print(f"  {s['scenario']}: nominal {nom:.2f}, max drift {w[0]} {w[1]:+.1f}%, "
          f"recheck {s['nominal_recheck']['P1_tlast']:.2f}")
