"""Chain composition of flip gadgets by identifying ONE clean socket vertex: s2 in I2 (not in the anchor member, no hard edge in the pattern)
of gadget g with s1 in I1 of gadget g+1.  Test by ILP: tight; relation between state(I1 of first) and state(I2 of last) with all anchors locked
versus one anchor released."""
import sys, glob, pickle, itertools, numpy as np, time
from w6c_lib import *
from w6c_compose2 import sfix, rel
def build(rec, N, s1, s2):
    mems, j, i1, i2, Rl, Rr, new, hard = rec
    n0 = 8; x0 = n0 - 1
    hedges = [p for i, p in enumerate(hard) if (j >> i) & 1]
    ident = {}
    for g in range(1, N): ident[(g, s1)] = (g - 1, s2)
    idx = {}
    def res(g, v):
        while (g, v) in ident: g, v = ident[(g, v)]
        return (g, v)
    def V(g, v):
        r = res(g, v)
        if r not in idx: idx[r] = len(idx)
        return idx[r]
    members = []; edges = set()
    for g in range(N):
        for mm in mems:
            t = tuple(V(g, v) for v in mm); members.append(t); edges |= {key(t[0], t[1]), key(t[1], t[2])}
        for (a, b) in hedges: edges.add(key(V(g, a), V(g, b)))
    n = len(idx); A = adj_from_edges(n, edges)
    xs = [V(g, x0) for g in range(N)]
    return A, members, xs, members[i1], members[(N - 1) * len(mems) + i2]
if __name__ == '__main__':
    recs = []
    for f in sorted(glob.glob('/tmp/claude-0/scratch/flip_8_5_*.pkl')): recs += pickle.load(open(f, 'rb'))
    part = int(sys.argv[1]); nparts = int(sys.argv[2]); N = int(sys.argv[3])
    cnt = 0
    for ri, rec in enumerate(recs):
        if ri % nparts != part: continue
        mems, j, i1, i2, Rl, Rr, new, hard = rec
        hedges = [p for i, p in enumerate(hard) if (j >> i) & 1]
        hv = {v for p in hedges for v in p}
        anchor_vs = set(mems[0])
        for s1 in mems[i1]:
            for s2 in mems[i2]:
                if s1 in anchor_vs or s2 in anchor_vs or s1 in hv or s2 in hv: continue
                cnt += 1
                A, H, xs, I1f, I2l = build(rec, N, s1, s2)
                ok, why = check_packing(A, H)
                if not ok: continue
                if ilp_opt(A)[0] != len(H): continue
                allL = rel(A, H, xs, xs, I1f, I2l)
                rels = [rel(A, H, xs, [x for x in xs if x != xs[g]], I1f, I2l) for g in range(N)]
                flip = len(allL) >= 1 and all(any(p not in allL for p in R) for R in rels)
                print(ri, (s1, s2), 'allL', sorted(allL), 'rels', [sorted(R) for R in rels], 'FLIP' if flip else '', flush=True)
    print('done part', part, 'candidates', cnt, flush=True)
