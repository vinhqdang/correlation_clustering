"""Chain composition of flip gadgets found by w6c_flip4.py by IDENTIFYING the right interface member I2 of gadget g with the left interface
member I1 of gadget g+1 (vertex map: centre->centre, ends in either order).  Tests by ILP the composite relation between
state(I1 of first gadget) and state(I2 of last gadget): all anchors locked vs one anchor released."""
import sys, glob, pickle, itertools, numpy as np, time
from w6c_lib import *
ST = {0: (0, 1, 0), 1: (1, 0, 0), 2: (1, 1, 1)}
def sfix(mem, s):
    a, b, c = mem; ab, bc, ac = ST[s]; return {key(a, b): ab, key(b, c): bc, key(a, c): ac}
def build(rec, N, endsflip):
    mems, j, i1, i2, Rl, Rr, new, hard = rec
    n0 = 8; x0 = n0 - 1
    # hard-edge pattern of the gadget
    fe = set()
    for mm in mems: fe |= {key(mm[0], mm[1]), key(mm[1], mm[2])}
    hedges = [p for i, p in enumerate(hard) if (j >> i) & 1]
    I1 = mems[i1]; I2 = mems[i2]
    # identification: gadget g+1 vertex I1[k] := gadget g vertex I2[phi(k)]
    phi = {1: 1, 0: (0 if not endsflip else 2), 2: (2 if not endsflip else 0)}
    ident = {}
    for g in range(1, N):
        for k in range(3): ident[(g, I1[k])] = (g - 1, I2[phi[k]])
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
    I1f = members[i1]; I2l = members[(N - 1) * len(mems) + i2]
    return A, members, xs, I1f, I2l
def rel(A, H, xs, locked, I1f, I2l, states=(0, 1, 2)):
    n = len(A); out = set()
    for s in states:
        for t in states:
            fx = {}; fx.update(sfix(I1f, s)); fx.update(sfix(I2l, t))
            for x in locked:
                for u in range(n):
                    if u != x: fx[key(x, u)] = 0
            r = ilp_opt(A, fixed=fx)
            if r is not None and r[0] == len(H): out.add((s, t))
    return out
def test(rec, N, ef):
    A, H, xs, I1f, I2l = build(rec, N, ef)
    ok, why = check_packing(A, H)
    if not ok: return ('invalid', why)
    if ilp_opt(A)[0] != len(H): return ('notight',)
    allL = rel(A, H, xs, xs, I1f, I2l)
    rels = [rel(A, H, xs, [x for x in xs if x != xs[g]], I1f, I2l) for g in range(N)]
    return ('ok', allL, rels)
if __name__ == '__main__':
    recs = []
    for f in glob.glob(outpath('flip_8_5_*.pkl')): recs += pickle.load(open(f, 'rb'))
    # dedupe by (members, j, i1, i2)
    seen = set(); uniq = []
    for r in recs:
        k = (tuple(r[0]), r[1], r[2], r[3])
        if k in seen: continue
        seen.add(k); uniq.append(r)
    print(len(recs), 'records,', len(uniq), 'unique', flush=True)
    part = int(sys.argv[1]); nparts = int(sys.argv[2]); N = int(sys.argv[3])
    for ri, rec in enumerate(uniq):
        if ri % nparts != part: continue
        for ef in (False, True):
            r = test(rec, N, ef)
            if r[0] != 'ok': continue
            allL, rels = r[1], r[2]
            flip = len(allL) >= 1 and all(any(p not in allL for p in R) for R in rels)
            print(ri, ef, 'allL', sorted(allL), 'rels', [sorted(R) for R in rels], 'FLIP' if flip else '', flush=True)
