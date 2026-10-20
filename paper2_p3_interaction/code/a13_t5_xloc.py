"""a13_t5_xloc.py: wave-8 definition of locality (X_u cap X_v of a path decomposition of I(H)) with the PARTITION variant of F.
Checks whether T_N alone, or T_N plus the pairs that are far only under X-locality, is a minimal certificate."""
import sys, itertools
from a13_lib import *
from ortools.sat.python import cp_model
def Fpart(inst, cnt, T):
    n = inst['n']; A = adjm(n, inst['E']); Ts = set(T)
    m = cp_model.CpModel(); z = {}; cost = []
    for u, v in itertools.combinations(range(n), 2):
        p = frozenset((u, v)); z[p] = m.NewBoolVar('')
        if p in cnt: cost.append(z[p] if not A[u][v] else 1 - z[p])
        if p in Ts: m.Add(z[p] == A[u][v])
    for x, y, w in itertools.combinations(range(n), 3):
        a, b, c = z[frozenset((x, y))], z[frozenset((y, w))], z[frozenset((x, w))]
        m.AddBoolOr([a.Not(), b.Not(), c]); m.AddBoolOr([a.Not(), c.Not(), b]); m.AddBoolOr([b.Not(), c.Not(), a])
    m.Minimize(sum(cost)); s = cp_model.CpSolver(); s.parameters.num_workers = 4
    st = s.Solve(m)
    return None if st == cp_model.INFEASIBLE else int(round(s.ObjectiveValue()))
def xlocal(inst):
    mem = inst['members']; N = len(mem)//3
    # path decomposition: bags {A_j,B_j,C_j}, {B_j,A_{j+1}}
    bags = []
    for j in range(N):
        bags.append(set([3*j, 3*j+1, 3*j+2]))
        if j+1 < N: bags.append(set([3*j+1, 3*(j+1)]))
    n = inst['n']
    X = {}
    for u in range(n):
        M = set(i for i, m in enumerate(mem) if u in m)
        if M: X[u] = set(t for t, b in enumerate(bags) if M <= b)
    loc = set()
    for u, v in itertools.combinations(range(n), 2):
        if u in X and v in X and X[u] & X[v]: loc.add(frozenset((u, v)))
    soft = set(p for m in mem for p in pairs_of(m))
    return loc | soft
for N in range(1, int(sys.argv[1])+1):
    inst = build_GN(N); H = len(inst['members']); T = inst['T']
    cx = xlocal(inst); cm_ = counted_pairs(inst)
    links = sorted(p for p in cm_ - cx)           # local for member-adjacency, far for X-locality
    far_x = set(frozenset(p) for p in itertools.combinations(range(inst['n']), 2)) - cx
    print('N=%d  pairs local(adjacency) but far(X): %d ; T_N subset of X-far: %s' % (N, len(links), all(t in far_x for t in T)))
    for name, TT in [('T_N', T), ('T_N+links', T+links)]:
        f = Fpart(inst, cx, TT)
        sub = [Fpart(inst, cx, [t for t in TT if t != e]) for e in TT]
        print('   %-10s |T|=%d  F_X(T)-|H|=%s  minimal=%s  (F_X(T-e)-|H| values %s)' % (name, len(TT), None if f is None else f-H, f is not None and f > H and all(s is not None and s <= H or s is None for s in sub), sorted(set(None if s is None else s-H for s in sub), key=str)), flush=True)
print('--- greedy minimalisation of T_N + links under X-locality (links deleted first) ---')
for N in range(2, int(sys.argv[1])+1):
    inst = build_GN(N); H = len(inst['members']); T = inst['T']
    cx = xlocal(inst); links = sorted(counted_pairs(inst) - cx)
    cur = list(links) + list(T)
    for e in list(cur):
        t2 = [x for x in cur if x != e]
        f = Fpart(inst, cx, t2)
        if f is None or f > H: cur = t2
    f = Fpart(inst, cx, cur)
    sub = [Fpart(inst, cx, [t for t in cur if t != e]) for e in cur]
    print('N=%d minimal X-certificate size %d (3N+3=%d, 2N+4=%d); F-|H|=%d ; minimal=%s ; links kept=%d' % (N, len(cur), 3*N+3, 2*N+4, f-H, all(s is not None and s <= H for s in sub), len([c for c in cur if c in set(links)])), flush=True)
