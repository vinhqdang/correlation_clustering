"""a13_lib.py -- independent re-implementation (written from the description in notes/wave9_corebound.md Sec. 5.1 and Sec. 1 only).
Does not import any w8*/w9* file."""
import itertools, subprocess, os
from ortools.sat.python import cp_model

HERE = os.path.dirname(os.path.abspath(__file__))

def build_GN(N):
    """vertices p_j,q_j,r_j,c_j,b_j (j=1..N), a_1 extra, a_j=b_{j-1}; pendants w_1..w_{N+2}.  Returns dict."""
    names = []
    idx = {}
    def v(nm):
        if nm not in idx:
            idx[nm] = len(names); names.append(nm)
        return idx[nm]
    def a(j): return v('a1') if j == 1 else v('b%d' % (j-1))
    members = []   # (end1, centre, end2) , G-edges end1-centre, centre-end2, non-edge end1-end2
    for j in range(1, N+1):
        p, q, r, c, b = v('p%d'%j), v('q%d'%j), v('r%d'%j), v('c%d'%j), v('b%d'%j)
        members.append((q, p, a(j)))   # A_j
        members.append((q, r, b))      # B_j
        members.append((p, c, r))      # C_j
    E = set()
    def e(x, y): E.add(frozenset((x, y)))
    for (x, y, z) in members:
        e(x, y); e(y, z)
    for j in range(1, N+1):
        e(v('c%d'%j), v('q%d'%j))                   # hard edge c_j q_j
    T = []
    w = {}
    for j in range(1, N+1):
        w[j] = v('w%d'%j); e(w[j], v('p%d'%j))
        T += [frozenset((w[j], v('p%d'%j))), frozenset((w[j], v('r%d'%j)))]
    w[N+1] = v('w%d'%(N+1)); e(w[N+1], v('a1'))
    T += [frozenset((w[N+1], v('a1'))), frozenset((w[N+1], v('p1')))]
    w[N+2] = v('w%d'%(N+2)); e(w[N+2], v('b%d'%N))
    T += [frozenset((w[N+2], v('b%d'%N))), frozenset((w[N+2], v('r%d'%N)))]
    return dict(n=len(names), names=names, idx=idx, E=E, members=members, T=T)

def adjm(n, E):
    A = [[0]*n for _ in range(n)]
    for e in E:
        x, y = tuple(e); A[x][y] = A[y][x] = 1
    return A

def pairs_of(m):
    a, b, c = m
    return [frozenset((a, b)), frozenset((b, c)), frozenset((a, c))]

def validate(inst, verbose=True):
    n, E, members = inst['n'], inst['E'], inst['members']
    A = adjm(n, E)
    ok = True
    # (1) every member an induced P3 with centre = middle entry
    for (a, b, c) in members:
        if not (A[a][b] and A[b][c] and not A[a][c]): ok = False; print('not induced P3', (a, b, c))
    # (2) modification-disjoint: pair sets disjoint
    seen = {}
    for i, m in enumerate(members):
        for p in pairs_of(m):
            if p in seen: ok = False; print('pair shared', p, i, seen[p])
            seen[p] = i
    # (3) pair-maximal: no induced P3 of G whose three pairs are all outside soft pairs
    soft = set(seen)
    bad = []
    for x, y, z in itertools.combinations(range(n), 3):
        ps = [frozenset((x, y)), frozenset((y, z)), frozenset((x, z))]
        s = sum(A[u][v] for u, v in [(x, y), (y, z), (x, z)])
        if s == 2 and not any(p in soft for p in ps): bad.append((x, y, z))
    if bad: ok = False; print('not pair-maximal; free induced P3s:', bad[:5])
    # (4) c = max number of members through a vertex
    cnt = [0]*n
    for m in members:
        for x in m: cnt[x] += 1
    return ok, max(cnt)

def IH(inst):
    m = inst['members']
    adj = {i: set() for i in range(len(m))}
    for i, j in itertools.combinations(range(len(m)), 2):
        if set(m[i]) & set(m[j]): adj[i].add(j); adj[j].add(i)
    return adj

def local_pairs(inst):
    """pair {u,v} local iff u,v in V(H) and some member through u equals or is I(H)-adjacent to some member through v"""
    n, members = inst['n'], inst['members']
    adj = IH(inst)
    M = [[i for i, m in enumerate(members) if u in m] for u in range(n)]
    loc = set()
    for u, v in itertools.combinations(range(n), 2):
        if not M[u] or not M[v]: continue
        if any(i == j or j in adj[i] for i in M[u] for j in M[v]): loc.add(frozenset((u, v)))
    return loc

def counted_pairs(inst):
    soft = set(p for m in inst['members'] for p in pairs_of(m))
    return soft | local_pairs(inst)

def Fprime(inst, T, budget_hint=None, extra_counted=None, return_sol=False):
    """csp-variant: min #disagreements over counted pairs s.t. every triple all of whose far pairs lie in T is not a P3.
    Far pairs in T are fixed to their G status; far pairs not in T are free & uncosted; triples with such a pair are dropped."""
    n, E = inst['n'], inst['E']
    A = adjm(n, E)
    cnt = counted_pairs(inst) if extra_counted is None else extra_counted
    Tset = set(T)
    mdl = cp_model.CpModel()
    z = {}
    cost = []
    for u, v in itertools.combinations(range(n), 2):
        p = frozenset((u, v))
        z[p] = mdl.NewBoolVar('z%d_%d' % (u, v))
        if p in cnt:
            cost.append(z[p] if not A[u][v] else 1 - z[p])
        elif p in Tset:
            mdl.Add(z[p] == A[u][v])
    for x, y, zz in itertools.combinations(range(n), 3):
        ps = [frozenset((x, y)), frozenset((y, zz)), frozenset((x, zz))]
        if any((p not in cnt) and (p not in Tset) for p in ps): continue
        a, b, c = [z[p] for p in ps]
        mdl.AddBoolOr([a.Not(), b.Not(), c])
        mdl.AddBoolOr([a.Not(), c.Not(), b])
        mdl.AddBoolOr([b.Not(), c.Not(), a])
    mdl.Minimize(sum(cost))
    s = cp_model.CpSolver(); s.parameters.num_workers = 4; s.parameters.max_time_in_seconds = 600
    st = s.Solve(mdl)
    if st == cp_model.INFEASIBLE: return None
    assert st == cp_model.OPTIMAL, st
    return int(round(s.ObjectiveValue()))

def dp_partition(inst, counted, T, n_override=None):
    """exact partition DP in C (n<=22).  counted: set of pairs with cost; T: iterable of pairs (edge->must join, non-edge->must separate)"""
    n, E = inst['n'], inst['E']
    A = adjm(n, E)
    C = [[1 if frozenset((u, v)) in counted else 0 for v in range(n)] for u in range(n)]
    Tm = [[0]*n for _ in range(n)]
    for p in T:
        u, v = tuple(p); code = 1 if A[u][v] else 2
        Tm[u][v] = Tm[v][u] = code
    txt = '%d\n' % n + '\n'.join(' '.join(map(str, r)) for r in A) + '\n' + '\n'.join(' '.join(map(str, r)) for r in C) + '\n' + '\n'.join(' '.join(map(str, r)) for r in Tm) + '\n'
    out = subprocess.run([os.path.join(HERE, 'a13_dp')], input=txt, capture_output=True, text=True).stdout
    return int(out.strip())

def opt_exact(inst):
    n = inst['n']
    allp = set(frozenset(p) for p in itertools.combinations(range(n), 2))
    return dp_partition(inst, allp, [])
