"""Inspect small tight instances whose zero-excess solution set is not closed under any majority: is it 2-decomposable, and where are the 3-cores of member states?"""
import itertools, numpy as np
from w6c_lib import *
from w6c_majority import states
EX = [([(0, 1), (0, 6), (1, 3), (1, 5), (1, 6), (2, 4), (2, 7), (3, 5), (3, 6), (3, 7), (5, 6), (5, 7)], [(0, 1, 3), (4, 2, 7), (6, 3, 7), (1, 5, 7), (0, 6, 5)], 8),
      ([(0, 3), (0, 5), (0, 8), (1, 7), (2, 4), (3, 5), (3, 6), (3, 7), (3, 8), (4, 5), (4, 6), (4, 8), (5, 6), (5, 8)], [(0, 3, 6), (5, 3, 7), (2, 4, 6), (0, 5, 4), (6, 5, 8), (3, 8, 4)], 9),
      ([(0, 1), (0, 4), (0, 6), (0, 7), (0, 8), (1, 3), (1, 4), (1, 7), (1, 8), (2, 6), (3, 8), (5, 7), (5, 8), (7, 8)], [(4, 0, 7), (0, 1, 3), (4, 1, 8), (0, 6, 2), (1, 7, 5), (0, 8, 5), (3, 8, 7)], 9)]
names = 'LRM'
for E, H, n in EX:
    A = adj_from_edges(n, E)
    ok, why = check_packing(A, H)
    opt, lab, same, cost, P = brute_opt(A)
    Sig = sorted({states(lab[i], H) for i in np.nonzero(cost == opt)[0]})
    p = len(H)
    print('instance n', n, 'members', H, 'valid', ok, 'OPT', opt, '|H|', len(H), '|Sigma|', len(Sig), 'pair-maximal', is_pair_maximal(A, H))
    for s in Sig: print('   ', ''.join(names[v] for v in s))
    # 2-decomposable closure
    proj = {(i, j): {(s[i], s[j]) for s in Sig} for i in range(p) for j in range(i + 1, p)}
    S2 = [t for t in itertools.product(range(3), repeat=p) if all((t[i], t[j]) in proj[(i, j)] for i in range(p) for j in range(i + 1, p))]
    extra = [t for t in S2 if t not in set(Sig)]
    print('   pairwise-consistent but not solutions:', [''.join(names[v] for v in t) for t in extra][:10], 'count', len(extra))
    # clean vertices (only soft edges, in one member as an end)
    used = set()
    for t in H: used |= pairs_of(t)
    inm = {}
    for i, m in enumerate(H):
        for pos, v in enumerate(m): inm.setdefault(v, []).append((i, pos))
    clean = [v for v in range(n) if all(key(v, u) in used for u in range(n) if A[v][u]) and v in inm and all(pos in (0, 2) for _, pos in inm[v])]
    print('   clean anchor-capable vertices:', clean, [(v, inm[v]) for v in clean])
