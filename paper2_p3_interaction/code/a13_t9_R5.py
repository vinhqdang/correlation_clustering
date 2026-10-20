"""a13_t9_R5.py: realise Rk in a real packing instance (outside members own the 'dropped' pairs) and re-check the R5 argument for r<=7 by brute force."""
import itertools
from a13_lib import *
Rk = {(0,0,0),(0,0,2),(0,1,1),(0,2,0),(1,0,1),(1,1,0),(1,1,1),(1,2,1),(2,0,0),(2,2,0)}
A_, B_, C_ = (2,1,0), (4,3,1), (5,0,3)           # (end, centre, end) on vertices 0..5
dropped = [(2,4),(2,3),(4,0),(5,1)]; nonedge = [(2,5),(4,5)]
E = set()
for (a,b,c) in (A_,B_,C_): E.add(frozenset((a,b))); E.add(frozenset((b,c)))
mem = [A_, B_, C_]; nv = 6
for (s, t) in dropped:                          # outside member D=(s,x,t) with new centre x, owns the pair st (a non-edge ac pair)
    x = nv; nv += 1; mem.append((s, x, t)); E.add(frozenset((s, x))); E.add(frozenset((x, t)))
inst = dict(n=nv, E=E, members=mem, T=[])
ok, c = validate(inst)
print('realising instance: n=%d members=%d  valid(induced, pair-disjoint, [pair-maximal flag below])=%s  c=%d' % (nv, len(mem), ok, c))
A = adjm(nv, E)
for (s, t) in nonedge: assert not A[s][t]
for (s, t) in dropped: assert not A[s][t]
soft = {}
for i, m in enumerate(mem):
    for p in pairs_of(m): soft[p] = i
rel = set()
for st in itertools.product(range(3), repeat=3):
    adj = {}
    for i in range(3):
        a, b, c_ = mem[i]; s = st[i]
        adj[frozenset((a,b))] = 0 if s == 0 else 1; adj[frozenset((b,c_))] = 0 if s == 1 else 1; adj[frozenset((a,c_))] = 1 if s == 2 else 0
    good = True
    for tr in itertools.combinations(range(6), 3):
        ps = [frozenset(p) for p in itertools.combinations(tr, 2)]
        if any(p in soft and soft[p] >= 3 for p in ps): continue     # triple contains a pair owned by an outside member: belongs to another constraint
        tot = sum(adj[p] if p in adj else A[tuple(p)[0]][tuple(p)[1]] for p in ps)
        if tot == 2: good = False; break
    if good: rel.add(st)
print('relation on (A,B,C) =', sorted(rel)); print('equals Rk:', rel == Rk)
# R5 argument: Boolean restriction and the union-closure proof, brute force over all conservative f on {0,1}^r for r=3,4 (f(x) in {x_i}, idempotent)
R5 = {t for t in itertools.product((0,1), repeat=3) if sum(t) != 1}
for r in (3, 4):
    pts = list(itertools.product((0,1), repeat=r)); cnt = 0
    # enumerate f via backtracking on constraints is heavy; use SAT-free search: f is a function pts->{0,1} with f(x)=x_i forced when x constant
    import itertools as it
    free = [p for p in pts if 0 < sum(p) < r]
    found = False
    for bits in it.product((0,1), repeat=len(free)):
        f = dict(zip(free, bits)); f[tuple([0]*r)] = 0; f[tuple([1]*r)] = 1
        if any(f[tuple(1 if j == i else 0 for j in range(r))] != 0 for i in range(r)): continue       # NU: f(e_i)=0
        if any(f[tuple(0 if j == i else 1 for j in range(r))] != 1 for i in range(r)): continue       # NU: f(1-e_i)=1
        okp = True
        for rows in it.product(sorted(R5), repeat=r):
            if tuple(f[tuple(rows[i][c] for i in range(r))] for c in range(3)) not in R5: okp = False; break
        if okp: found = True; break
    print('arity r=%d: Boolean conservative NU preserving R5 exists: %s' % (r, found))
