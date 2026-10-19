"""w6b_far.py -- search for instances with CONNECTED I(H) (path-like, w=1, c=2) in which EVERY optimum edits a hard pair whose
ends lie in members at I(H)-distance >= K (K = 2, 3), i.e. the pointer is not contained in a bag of any tree decomposition of
I(H)^(K-1).  Base graph: a path of m members M_i = (2i, 2i+1, 2i+2) (consecutive members share a vertex) + extra vertices + random
extra edges; H = the chain, extended to a pair-maximal packing by greedily adding induced P3s (these extra members may change
I(H)).  Exhaustive optimum (all set partitions, n <= 10).
Usage: python3 w6b_far.py m n K trials seed
"""
import sys, random, itertools
import numpy as np
from w6b_lib import *

def build(m, n, rng, pe):
    A = np.zeros((n, n), dtype=np.int8)
    chain = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
    for (x, y, z) in chain:
        A[x, y] = A[y, x] = A[y, z] = A[z, y] = 1
    base = 2 * m + 1
    for a, b in pair_list(n):
        if A[a, b]: continue
        # skip pairs inside a member (they are packed)
        if any(a in t and b in t for t in chain): continue
        if rng.random() < pe: A[a, b] = A[b, a] = 1
    return A, chain

def extend_pm(adj, H, rng):
    H = list(H); used = packed_set(H)
    ts = induced_p3s(adj); rng.shuffle(ts)
    for t in ts:
        p = pairs3(t)
        if not p & used:
            H.append(t); used |= p
    return H

def main():
    m, n, K, trials, seed = [int(x) for x in sys.argv[1:6]]
    rng = random.Random(seed)
    found = []
    ok_inst = 0
    for _ in range(trials):
        pe = rng.choice([0.08, 0.12, 0.18, 0.25])
        A, chain = build(m, n, rng, pe)
        if not is_packing(A, chain): continue
        H = extend_pm(A, chain, rng)
        if not (is_packing(A, H) and pair_maximal(A, H)): continue
        inst = Inst(A, H)
        # connected I(H)?
        mm = len(H); seen = {0}; st = [0]
        while st:
            u = st.pop()
            for v in inst.I[u]:
                if v not in seen: seen.add(v); st.append(v)
        if len(seen) != mm: continue
        ok_inst += 1
        inst.solve()
        opts = inst.opt_partitions()
        good = True
        for lab in opts:
            U, S = inst.edits(lab)
            if not any(inst.vertex_dist(a, b) >= K for a, b in U): good = False; break
        if good:
            found.append((A, H, int(inst.cost.min()) - len(H), len(opts)))
    print('m=%d n=%d K=%d: connected pair-maximal instances=%d, every-optimum-has-edit-at-distance>=K: %d' % (m, n, K, ok_inst, len(found)))
    found.sort(key=lambda r: (len(r[1]), r[2], r[3]))
    for A, H, l, no in found[:4]:
        edges = [(i, j) for i in range(n) for j in range(i + 1, n) if A[i, j]]
        print('  |H|=%d l=%d #opt=%d edges=%s H=%s' % (len(H), l, no, edges, H))

if __name__ == '__main__':
    main()
