"""Audit: C1, C3, C4 example, forced multi-edits, grid."""
import itertools, random
import numpy as np
from audit_lib import *

def opt_of(n, E):
    S, c = opt_and_costs(n, E); return int(c.min()), S, c

# C3(a): stars
for m in (1, 2, 3, 4):
    n = 2 * m + 1
    edges = [(0, i) for i in range(1, n)]
    E = adj_mask(n, edges)
    H = [(2 * i + 1, 0, 2 * i + 2) for i in range(m)]
    check_packing(n, E, H)
    o, _, _ = opt_of(n, E)
    I = {(i, j) for i in range(m) for j in range(i + 1, m) if set(H[i]) & set(H[j])}
    print('star m=%d OPT=%d expected 2m-1=%d  l=%d  I(H) edges=%d (K_m has %d)' % (m, o, 2 * m - 1, o - m, len(I), m * (m - 1) // 2))

# C3(b) 3x3 grid, row P3s vertex-disjoint; I(H) edgeless; tw(grid3)=3 (known)
n = 9
edges = []
for r in range(3):
    for c in range(3):
        v = 3 * r + c
        if c < 2: edges.append((v, v + 1))
        if r < 2: edges.append((v, v + 3))
E = adj_mask(n, edges)
H = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
check_packing(n, E, H)
print('grid3 packing ok; OPT', opt_of(n, E)[0], 'packing size 3')

# C4 example of fpt_attempt.md
n = 4
E = adj_mask(n, [(0, 1), (0, 2), (1, 3), (2, 3)])
H = [(0, 1, 3)]
pk = check_packing(n, E, H)
o, S, c = opt_of(n, E)
print('C4: OPT', o, '|H|', len(H))
P = pidx(n)
for s in S[c == o]:
    s = int(s)
    print('  optimum', [(u, v) for u in range(n) for v in range(u + 1, n) if (s >> P[(u, v)]) & 1],
          'unpacked part', [(u, v) for u in range(n) for v in range(u + 1, n) if (s >> P[(u, v)]) & 1 and not (pk >> P[(u, v)]) & 1])
print('  l=0 CSP satisfiable?', any(is_cluster(n, E ^ (1 << P[q])) for q in [(0, 1), (1, 3), (0, 3)]))
print('  H maximal?', all(pairmask(n, q) & pk for q in induced_p3s(n, E)))

# C1: arbitrary CE instances map to CEaMP with maximal vertex-disjoint packing, l = k - |H|
rng = random.Random(5); ok = 0
for _ in range(300):
    n = rng.randint(4, 8)
    E = rand_graph(rng, n, rng.choice(['gnp', 'cliques']))
    used = set(); H = []
    for (a, b, cc) in induced_p3s(n, E):
        if not ({a, b, cc} & used): H.append((a, b, cc)); used |= {a, b, cc}
    check_packing(n, E, H)
    I = [(i, j) for i in range(len(H)) for j in range(i + 1, len(H)) if set(H[i]) & set(H[j])]
    assert not I
    o = opt_of(n, E)[0]
    for k in range(0, 6):
        l = k - len(H)
        # CEaMP asks |S|-|H|<=l  <=> |S|<=k.  YES iff o<=k.
        assert (o - len(H) <= l) == (o <= k)
    ok += 1
print('C1 identity ok on', ok)

# minimal instances where every optimum needs a P3 with 2 / 3 edits
rng = random.Random(77)
best2 = best3 = None
for t in range(40000):
    n = rng.randint(4, 7)
    E = rand_graph(rng, n, rng.choice(['gnp', 'hub', 'cliques', 'path']))
    H = rand_packing(rng, n, E, rng.choice(['max', 'part', 'overlap']))
    if not H: continue
    pk = check_packing(n, E, H)
    o, S, c = opt_of(n, E)
    sel = S[c == o]
    mx = np.zeros(len(sel), dtype=np.int64)
    for q in H:
        mx = np.maximum(mx, popcount(sel & np.uint64(pairmask(n, q))).astype(np.int64))
    if mx.min() >= 2 and (best2 is None or (n, len(H)) < best2[0]): best2 = ((n, len(H)), n, E, H, o)
    if mx.min() >= 3 and (best3 is None or (n, len(H)) < best3[0]): best3 = ((n, len(H)), n, E, H, o)
def show(x):
    if x is None: print('  none found'); return
    _, n, E, H, o = x
    print('  n=%d edges=%s H=%s OPT=%d |H|=%d' % (n, [(u, v) for u in range(n) for v in range(u + 1, n) if has(n, E, u, v)], H, o, len(H)))
print('forced 2-edit (every optimum has a P3 with >=2 edits):'); show(best2)
print('forced 3-edit:'); show(best3)
