"""a13_t10_conn.py: a CONNECTED family for l=1 (so component-wise preprocessing cannot remove it): chain G_{2N} (blocks 1..2N, joined at b_N=a_{N+1})
plus one extra pendant w_m adjacent to b_N only, far pairs w_m b_N, w_m r_N, w_m p_{N+1}.  T = T_{2N} + these 3 pairs.
Checks validity (induced, pair-disjoint, pair-maximal, connected), F'(T)-|H| = 2 > l = 1, and minimality at l=1."""
import sys, itertools, networkx as nx
from a13_lib import *
for N in (1, 2, 3, 4):
    M = 2*N
    inst = build_GN(M); idx = inst['idx']
    wm = inst['n']; n = wm + 1
    E = set(inst['E']); E.add(frozenset((wm, idx['b%d' % N])))
    T = list(inst['T']) + [frozenset((wm, idx['b%d' % N])), frozenset((wm, idx['r%d' % N])), frozenset((wm, idx['p%d' % (N+1)]))]
    inst2 = dict(inst, n=n, E=E, T=T)
    ok, c = validate(inst2); H = len(inst2['members'])
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(tuple(e) for e in E)
    far = set(frozenset(p) for p in itertools.combinations(range(n), 2)) - counted_pairs(inst2)
    assert all(t in far for t in T)
    FT = Fprime(inst2, T) - H
    sub = [Fprime(inst2, [t for t in T if t != e]) - H for e in T]
    print('N=%d (2N=%d blocks) n=%d valid=%s c=%d connected(G)=%s |T|=%d F\'(T)-|H|=%d  max_e F\'(T-e)-|H|=%d  => minimal certificate for l=1: %s' % (N, M, n, ok, c, nx.is_connected(G), len(T), FT, max(sub), FT > 1 and max(sub) <= 1), flush=True)
