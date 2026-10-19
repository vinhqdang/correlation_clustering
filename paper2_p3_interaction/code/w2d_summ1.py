import json, collections, statistics as st
R = json.load(open("instances_w2d/exp1_results_9_13_40.json"))
print("records", len(R))
for tag in ("rand1", "rand2", "max"):
    print("==", tag)
    by = collections.defaultdict(list)
    for r in R:
        if r["tag"] == tag: by[r["fam"]].append(r)
    for f, rs in by.items():
        ls = [r["l"] for r in rs]
        print(f"{f:12s} n={len(rs):3d} l mean {st.mean(ls):.2f} max {max(ls)} P(l=0) {sum(x==0 for x in ls)/len(ls):.2f} P(l<=1) {sum(x<=1 for x in ls)/len(ls):.2f} P(l<=2) {sum(x<=2 for x in ls)/len(ls):.2f}")
    byn = collections.defaultdict(list)
    for r in R:
        if r["tag"] == tag: byn[r["n"]].append(r["l"])
    print(" by n mean l:", {n: round(st.mean(v), 2) for n, v in sorted(byn.items())})
# C1
c1 = [r for r in R if "minRR" in r]
print("C1-tested (l>0,|R|>=1):", len(c1), " with |R|>=2:", sum(r["R"] >= 2 for r in c1), " counterexamples minRR>0:", sum(r["minRR"] > 0 for r in c1))
print("C1 tested by tag/n:", collections.Counter((r["tag"], r["n"]) for r in c1 if r["R"] >= 2))
rt = [r for r in c1 if r["minRT"] > 0]
print("forced R-touching unpacked edit (every optimum):", len(rt), "of", len(c1), collections.Counter(r["fam"] for r in rt))
print("  max minRT", max(r["minRT"] for r in c1), "; by n:", collections.Counter(r["n"] for r in rt))
# (4) endpoints
print("endpoint stats: key (l,c) -> max #endpoints, max #orbits, count")
d = collections.defaultdict(lambda: [0, 0, 0, 0])
for r in c1:
    k = (r["l"], r["c"]); x = d[k]; x[0] = max(x[0], len(r["ends"])); x[1] = max(x[1], r["end_orbits"]); x[2] += 1; x[3] = max(x[3], r["n"])
for k in sorted(d): print(k, d[k])
print("max orbits by l:", {l: max(r["end_orbits"] for r in c1 if r["l"] == l) for l in sorted({r["l"] for r in c1})})
print("max endpoints by l:", {l: max(len(r["ends"]) for r in c1 if r["l"] == l) for l in sorted({r["l"] for r in c1})})
print("max orbits / (2l):", max(r["end_orbits"] / (2 * r["l"]) for r in c1), " max endpoints by n at l=1:", {n: max([len(r['ends']) for r in c1 if r['l']==1 and r['n']==n] or [0]) for n in range(9,14)})
print("ends in R max", max(r["ends_in_R"] for r in c1))
