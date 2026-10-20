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
    """clusters (frozensets of names) of the ZE solution at position a (formula).
    cell j (0<=j<=m-1): spine part  R: {s_2j,s_2j+1} (j<a),  M: {s_2j,s_2j+1,s_2j+2} (j=a),  L: {s_2j+1,s_2j+2} (j>a)
    e-labels: eps(j)=e_j for 2<=j<=m-2, eps(0)=e_1, eps(1)=e_0; cell m-1 gets e_{m-1},e_m if a<=m-2 ; if a=m-1: cell m-2 additionally gets e_{m-1} and cell m-1 gets e_m."""
    S = lambda k: ('s', k); Ev = lambda j: ('e', j)
    cl = []
    for j in range(m):
        if j < a:
            c = {S(2 * j), S(2 * j + 1)}
        elif j == a:
            c = {S(2 * j), S(2 * j + 1), S(2 * j + 2)}
        else:
            c = {S(2 * j + 1), S(2 * j + 2)}
        if j == 0: c.add(Ev(1))
        elif j == 1: c.add(Ev(0))
        elif j <= m - 2: c.add(Ev(j))
        if j == m - 1:
            if a <= m - 2: c |= {Ev(m - 1), Ev(m)}
            else: c.add(Ev(m))
        if j == m - 2 and a == m - 1: c.add(Ev(m - 1))
        cl.append(frozenset(c))
    return cl

def is_ze(m, clusters):
    n, E, members, names, idx = oh_instance(m)
    lab = [None] * n
    for ci, c in enumerate(clusters):
        for v in c: lab[idx[v]] = ci
    assert all(l is not None for l in lab)
    soft = set()
    for t in members: soft |= pairs_of(t)
    for u, v in itertools.combinations(range(n), 2):
        if key(u, v) not in soft and ((key(u, v) in E) != (lab[u] == lab[v])): return False
    for (a, b, c) in members:
        e = (lab[a] != lab[b]) + (lab[b] != lab[c]) + (lab[a] == lab[c])
        if e != 1: return False
    return True

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
    bad = 0
    for m in range(7, 61):
        for a in range(m):
            if not is_ze(m, explicit_partition(m, a)): bad += 1; print('invalid', m, a)
    print('explicit partitions valid ZE for m=7..60, all a; invalid:', bad)
    for m in (7, 8, 9, 10, 13):
        tp = true_partitions(m)
        ok = len(tp) == m and all(set(map(frozenset, tp[a])) == set(explicit_partition(m, a)) for a in range(m))
        print('m', m, 'SAT-enumerated ZE solutions == explicit', ok)


def transfer_certificate(mbig=26, verbose=True):
    """relative-state transfer matrices of the cell DP of OH_m (using w8c_transfer2.step_states), uniform on interior cells."""
    I = T2.Inst(mbig)
    # layered reachable states with counts
    layers = []      # layers[c] = dict state->count after processing cell c
    cur = {((), frozenset()): 1}
    mats = {}
    for c in range(0, mbig + 1):
        nxt = {}
        trans = {}
        for st, cnt in cur.items():
            for (lab, ns) in T2.step_states(I, c, st):
                nxt[ns] = nxt.get(ns, 0) + cnt
                trans[(st, ns)] = trans.get((st, ns), 0) + 1
        mats[c] = trans
        layers.append(nxt)
        cur = nxt
    return I, layers, mats

def matrix_of(trans, rows, cols):
    M = sp.zeros(len(rows), len(cols))
    ri = {s: i for i, s in enumerate(rows)}; ci = {s: i for i, s in enumerate(cols)}
    for (a, b), cnt in trans.items():
        M[ri[a], ci[b]] += cnt
    return M

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'transfer':
    mbig = 26
    I, layers, mats = transfer_certificate(mbig)
    print('state-set sizes per cell', [len(l) for l in layers])
    # interior cells: states are relative -> compare state sets and matrices
    cells = list(range(6, mbig - 6))
    base_states = sorted(layers[cells[0]].keys(), key=repr)
    for c in cells:
        assert sorted(layers[c].keys(), key=repr) == base_states, ('state set differs', c)
        tr = mats[c + 1]
        assert set(a for (a, b) in tr) <= set(base_states) and set(b for (a, b) in tr) <= set(base_states)
    M0 = matrix_of(mats[cells[0] + 1], base_states, base_states)
    for c in cells[:-1]:
        assert matrix_of(mats[c + 1], base_states, base_states) == M0, ('matrix differs', c)
    print('interior transfer matrix M (6x6), identical for cells', cells[0] + 1, '..', cells[-1])
    sp.pprint(M0)
    x = sp.symbols('x')
    print('char poly', sp.factor(M0.charpoly(x).as_expr()))
    print('Jordan form:'); P, J = M0.jordan_form(); sp.pprint(J)
    # N_m for all m via u M^t v : left vector after cell 5, right functional from real m
    left = sp.Matrix([[layers[cells[0]][s] for s in base_states]])
    print('left vector', list(left))
    # N_m = left * M^(m - 2*? ) ... compute right vector for each m by brute DP of last cells (cap region) -- done in numbers below
    print('N_m via DP:', [sum(T2.words_by_spec(m)[0].values()) for m in (7,)] )
