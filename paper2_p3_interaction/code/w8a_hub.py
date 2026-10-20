"""Hub stress test: spine-like family where the outside vertex y is adjacent to c_i, c_j of two far spine positions with the same head
(instead of {u_0,c_0}).  Is NF** (A-local covered edits + hub condition) still attainable?  joint3 ILP."""
import itertools, sys
import networkx as nx
from w8a_lib import *
from w8a_ilp3 import joint3
def build(arcs, hub_nbrs):
    N = len(arcs) - 1
    idx = {}; nv = 0
    for i in range(N + 1):
        for nm in 'abcu': idx[(nm, i)] = nv; nv += 1
    z = nv; y = nv + 1; n = nv + 2
    K = {}
    for i, (t, h) in enumerate(arcs):
        K.setdefault(t, []).append(idx[('a', i)])
        for nm in 'bcu': K.setdefault(h, []).append(idx[(nm, i)])
    E = set()
    for v, vs in K.items():
        for x, w in itertools.combinations(vs, 2): E.add(key(x, w))
    H = []
    for i in range(N + 1):
        a, b, c, u = (idx[(nm, i)] for nm in 'abcu')
        E.add(key(a, b)); un = idx[('u', i + 1)] if i < N else z
        E.add(key(b, un)); H.append((a, b, c)); H.append((u, b, un))
    for (nm, i) in hub_nbrs: E.add(key(y, idx[(nm, i)]))
    return n, E, H, idx, y
configs = {
 'y~c_1,c_3 (far, same head)': ([(3, '*'), (2, 1), (1, 2), (3, 1)], [('c', 1), ('c', 3)]),
 'y~u_1,c_1,c_3': ([(3, '*'), (2, 1), (1, 2), (3, 1)], [('u', 1), ('c', 1), ('c', 3)]),
 'y~c_1,u_3': ([(3, '*'), (2, 1), (1, 2), (3, 1)], [('c', 1), ('u', 3)]),
 'y~u_0,c_0 (original trigger)': ([(3, '*'), (2, 1), (1, 2), (3, 1)], [('u', 0), ('c', 0)]),
}
for name, (arcs, hn) in configs.items():
    n, E, H, idx, y = build(arcs, hn)
    g = G8(n, list(E))
    ok = valid_packing(g, H); pm = pair_maximal(g, H)
    opt, _ = opt_ilp(g); nu, _ = nu_mod(g)
    v = joint3(g, opt, nu, 2, tl=120)
    print('%-32s n=%d |H|=%d valid=%s pairmax=%s OPT=%d nu=%d  NF** violations=%s' % (name, n, len(H), ok, pm, opt, nu, v), flush=True)
