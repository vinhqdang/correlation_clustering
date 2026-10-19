"""Audit 7: N*(F) over ALL admissible orderings for the triangle and several for K4 / K3,3 : T irreducible, subfamilies irreducible,
pair-maximality, and whether the ordering condition head(e_1) != tail(e_0) matters."""
import sys, itertools, random
import numpy as np, networkx as nx
import a7_lib as L, a7_ml as M, a7_n5 as N5

def admissible(edges0, rest, need_gap):
    out = []
    for perm in itertools.permutations(rest):
        seq = [edges0] + list(perm)
        heads = [h for t, h in seq]
        if any(heads[i] == heads[i + 1] for i in range(len(heads) - 1)): continue
        gap = (seq[1][1] == edges0[0])
        if need_gap is None or gap == need_gap: out.append(seq)
    return out

def test_order(seq, maxF, verbose=False, K4=False):
    A, H, idx, K, z, y = N5.build(seq); n = len(A)
    pm = L.is_pair_maximal(A, H) and L.is_packing(A, H)
    cm = M.comps_of(H)
    if K4: tight = None; red_T = None
    else:
        W, sols = M.tight_solutions(A, H, cm[0])
        tight = len(sols) > 0; red_T = any(M.R_ok(A, B, None, W) for S, B in sols)
    pair_of = {}
    for i, t in enumerate(H):
        for p in L.p3p(t): pair_of[p] = i
    red = []; tested = 0
    for r in range(1, maxF + 1):
        for F in itertools.combinations(range(len(H)), r):
            Wf = sorted(set(v for i in F for v in H[i]))
            if any(pair_of[p] not in F for p in itertools.combinations(Wf, 2) if p in pair_of): continue
            tested += 1
            W_, sols_ = M.tight_solutions(A, H, list(F))
            if any(M.R_ok(A, B_, None, W_) for S_, B_ in sols_): red.append(F)
    return pm, tight, red_T, tested, red

if __name__ == "__main__":
    e0 = (0, '*')
    tri = [(0, 1), (1, 2), (2, 0)]
    res = {}
    for gap in (False, True):
        seqs = admissible(e0, tri, gap)
        for seq in seqs[:2]:
            pm, tight, red_T, tested, red = test_order(seq, 5)
            print("triangle gap=%s" % gap, seq[1:], "pm", pm, "tight", tight, "T reducible", red_T, "subfamilies", tested, "reducible", len(red))
    # K4: vertices 0..3, orientation with in-degrees <=2: 0>1,0>2,1>2,1>3,2>3,3>0  (in: 0:1,1:1,2:2,3:2)
    K4 = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (3, 0)]
    rng = random.Random(5)
    for trial in range(8):
        t = 0   # e0 = (0,*)
        rest = K4[:]; rng.shuffle(rest)
        while True:
            seq = [e0] + rest; heads = [h for _, h in seq]
            if all(heads[i] != heads[i + 1] for i in range(len(heads) - 1)): break
            rng.shuffle(rest)
        gap = seq[1][1] == 0
        pm, tight, red_T, tested, red = test_order(seq, 4, K4=True)
        print("K4 gap=%s" % gap, seq[1:], "pm", pm, "tight", tight, "T reducible", red_T, "subfamilies(|F|<=5)", tested, "reducible", len(red), red[:3])
