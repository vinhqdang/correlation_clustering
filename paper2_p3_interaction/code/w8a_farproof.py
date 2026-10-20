"""Verify the explicit repack of Prop FAR for m = 3..14: H* = {P_0..P_{m-2}, T=(1,2m-1,2m-2)}, S* = edits of {{0,1},{2,3},...,{2m-2,2m-1},{2m}}."""
from w8a_lib import *
for m in range(3, 15):
    n = 2 * m + 1
    g = G8(n, [(i, i + 1) for i in range(2 * m)] + [(1, 2 * m - 1)])
    H0 = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
    Hs = H0[:m - 1] + [(1, 2 * m - 1, 2 * m - 2)]
    part = [[2 * i, 2 * i + 1] for i in range(m)] + [[2 * m]]
    S = edits_of_partition(g, part)
    he = hard_edits(Hs, S)
    I = Igraph(Hs)
    ok = (valid_packing(g, Hs) and pair_maximal(g, Hs) and len(S) == m + 1 and he == [(2 * m - 1, 2 * m)] and
          sorted(d for _, d in I.degree()) == [2] * m and nx.is_connected(I))
    cover = set(v for t in Hs for v in t)
    print('m=%2d valid+pairmax=%s |H*|=%d |S*|=%d hard edits=%s (2m=%d uncovered: %s) I(H*) = cycle: %s  tw<=%d  [nu=%d via ILP, OPT=%d]' %
          (m, valid_packing(g, Hs) and pair_maximal(g, Hs), len(Hs), len(S), he, 2 * m, 2 * m not in cover, ok, tw_exact(I), nu_mod(g)[0], opt_ilp(g)[0]))
