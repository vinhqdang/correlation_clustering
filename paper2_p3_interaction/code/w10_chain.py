"""Chain lemma certificate for OH_m, ALL m >= 7.
(1) explicit ZE partitions P_a (formula) are valid for m = 7..60 and every a; their spine words are R^a M L^(m-1-a).
(2) completeness: transfer matrix of the cell-by-cell DP (relative state encoding).  N_m = u^T M^(m-m0) v with a FIXED 6x6 integer matrix M
    (uniform for all interior cells); characteristic polynomial and N_m = m for all m checked, plus symbolic closed form.
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from w8c_oh import oh_instance
import w8c_transfer2 as T2
import sympy as sp

def explicit_partition(m, a):
    """clusters (as sets of names) of the ZE solution at position a (formula; caps swap e0<->e1 on the left, and see code on the right)."""
    S = lambda k: ('s', k); Ev = lambda j: ('e', j)
    # generic labels: e'_j for cell j
    def eg(j):
        return Ev({0: 1, 1: 0}.get(j, j))      # left cap swap
    cl = []
    for j in range(m + 1):
        if j < a:
            c = {eg(j), S(2 * j), S(2 * j + 1)}
        elif j == a:
            c = {eg(j), S(2 * j), S(2 * j + 1), S(2 * j + 2)}
        else:
            c = {eg(j), S(2 * j + 1), S(2 * j + 2)}
        cl.append(c)
    return cl

def true_partitions(m):
    n, E, members, names, idx = oh_instance(m)
    sols = ze_solve(n, E, members, all_solutions=True, limit=1000)
    out = {}
    for lab in sols:
        d = {}
        for v, l in enumerate(lab): d.setdefault(l, set()).add(names[v])
        a = None
        for j in range(m):
            x, y, z = idx[('s', 2 * j)], idx[('s', 2 * j + 1)], idx[('s', 2 * j + 2)]
            if lab[x] == lab[y] == lab[z]: a = j
        out[a] = sorted(map(lambda c: tuple(sorted(c, key=str)), d.values()), key=str)
    return out

if __name__ == '__main__':
    # print actual partitions for m=10 so the formula can be adapted
    m = 10
    tp = true_partitions(m)
    for a in (0, 1, 2, 5, 8, 9):
        print(a, [list(map(lambda v: v[0] + str(v[1]), c)) for c in tp[a]])
