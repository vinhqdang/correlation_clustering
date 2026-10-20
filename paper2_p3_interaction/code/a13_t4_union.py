"""a13_t4_union.py: the l>0 remark.  k disjoint copies of G_N at budget |H|+(k-1): the minimal certificate T^1 u ... u T^k has size k(2N+4)
(F'(T)-|H| = k > l = k-1, and F'(T-e)-|H| = k-1 for every e).  Also checks pair-maximality / validity of the union and that
a single copy at l=1 has NO certificate inside T_N (F'(T_N)-|H| = 1 is not > 1)."""
import sys, itertools
from a13_lib import *
def union(insts):
    off = 0; n = 0; E = set(); mem = []; T = []
    for I in insts:
        sh = lambda x: x + off
        E |= set(frozenset(map(sh, tuple(e))) for e in I['E'])
        mem += [tuple(map(sh, m)) for m in I['members']]
        T += [frozenset(map(sh, tuple(t))) for t in I['T']]
        off += I['n']
    return dict(n=off, E=E, members=mem, T=T)
for N, k in [(1, 2), (2, 2), (1, 3), (3, 2)]:
    inst = union([build_GN(N) for _ in range(k)])
    ok, c = validate(inst); H = len(inst['members']); T = inst['T']
    l = k - 1
    FT = Fprime(inst, T) - H
    sub = [Fprime(inst, [t for t in T if t != e]) - H for e in T]
    print('N=%d copies=%d n=%d valid=%s c=%d l=%d |T|=%d  F\'(T)-|H|=%d (>l: %s)  F\'(T-e)-|H| in %s (all<=l: %s)' % (N, k, inst['n'], ok, c, l, len(T), FT, FT > l, sorted(set(sub)), all(s <= l for s in sub)), flush=True)
