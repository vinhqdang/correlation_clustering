"""w6b_family.py -- systematic search of the family  'path of m members (P_{2m+1}, members (2i,2i+1,2i+2)) + 1 or 2 extra edges
(+ optional pendant vertices)'.  For each instance (pair-maximal extension of the chain, connected I(H)) report the number of optima
and the maximum over optima of the minimum I(H)-distance spanned by an edited hard pair.  Goal: a clean family in which EVERY
optimum edits a hard pair at I(H)-distance >= K for each K (w = 1, l = 1).
Usage: python3 w6b_family.py m extra_vertices n_extra_edges
"""
import sys, itertools, random
import numpy as np
from w6b_lib import *

def extend_pm(adj, H, rng):
    H = list(H); used = packed_set(H)
    ts = induced_p3s(adj); rng.shuffle(ts)
    for t in ts:
        p = pairs3(t)
        if not p & used:
            H.append(t); used |= p
    return H

def main():
    m = int(sys.argv[1]); xv = int(sys.argv[2]); ne = int(sys.argv[3])
    n = 2 * m + 1 + xv
    chain = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
    base = [(2 * i, 2 * i + 1) for i in range(m)] + [(2 * i + 1, 2 * i + 2) for i in range(m)]
    inchain = lambda a, b: any(a in t and b in t for t in chain)
    cand = [(a, b) for a, b in pair_list(n) if not inchain(a, b)]
    rng = random.Random(1)
    best = {}
    for extra in itertools.combinations(cand, ne):
        A = mk(n, base + list(extra))
        H = extend_pm(A, chain, rng)
        if not (is_packing(A, H) and pair_maximal(A, H)): continue
        inst = Inst(A, H)
        mm = len(H); seen = {0}; st = [0]
        while st:
            u = st.pop()
            for v in inst.I[u]:
                if v not in seen: seen.add(v); st.append(v)
        if len(seen) != mm: continue
        inst.solve(); opts = inst.opt_partitions()
        l = int(inst.cost.min()) - len(H)
        worst = 10 ** 6
        for lab in opts:
            U, S = inst.edits(lab)
            d = max([inst.vertex_dist(a, b) for a, b in U] or [0])   # farthest edit of this optimum
            worst = min(worst, d)
        key = (l, worst if worst < 10 ** 5 else 'inf')
        best.setdefault(key, []).append((extra, len(H), len(opts)))
    for key in sorted(best, key=lambda k: (k[0], str(k[1]))):
        v = best[key]
        print('l=%d every-optimum-has-edit-at-I(H)-distance>=%s : %d instances; e.g. extra=%s |H|=%d #opt=%d' % (key[0], key[1], len(v), v[0][0], v[0][1], v[0][2]))

if __name__ == '__main__':
    main()
