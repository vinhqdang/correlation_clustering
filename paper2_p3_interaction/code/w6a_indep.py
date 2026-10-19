"""Independent re-check (partition based, no tight_choices): for the whole-only chains, every proper connected subfamily F is irreducible:
for EVERY optimal partition of G[W_F], some outside vertex y has N(y) cap W_F neither empty nor a cluster.  Whole T has an optimal partition satisfying R."""
import itertools
from w6a_lib import *
ex = {
 3: (7, [(0,1),(0,5),(1,2),(1,3),(3,4),(3,5),(3,6),(4,6)], [(0,1,2),(1,3,4),(5,3,6)]),
 4: (9, [(0,1),(0,7),(1,2),(1,3),(1,7),(3,4),(3,5),(5,6),(5,7),(5,8),(6,8)], [(0,1,2),(1,3,4),(3,5,6),(7,5,8)]),
 5: (11, [(0,1),(0,5),(0,8),(1,2),(1,3),(1,9),(2,9),(3,4),(3,5),(3,6),(4,6),(5,7),(5,8),(7,9),(7,10)], [(0,1,2),(1,3,4),(5,3,6),(7,5,8),(9,7,10)]),
}
def reducible(adj, H, F):
    W = sorted({v for i in F for v in H[i]}); n = len(adj)
    if any(len(set(t) & set(W)) >= 2 for j, t in enumerate(H) if j not in F): return None   # not pair-closed
    idx = {v: i for i, v in enumerate(W)}
    sub = from_edges(len(W), {(idx[a], idx[b]) for a in W for b in W if a < b and adj[a] >> b & 1})
    parts, o = all_opt_partitions(sub)
    for part in parts:
        cl = [frozenset(W[i] for i in range(len(W)) if c >> i & 1) for c in part]
        if all(set(adj_nb) in ([set()] + [set(c) for c in cl]) for adj_nb in [ {x for x in W if adj[y] >> x & 1} for y in range(n) if y not in W]): return True
    return False
for m, (n, E, H) in ex.items():
    adj = from_edges(n, E)
    assert is_pair_maximal(adj, H) and opt(adj) == len(H)
    allF = [F for k in range(1, m) for F in itertools.combinations(range(m), k)]
    from w6a_rule import conn_subfamilies
    conn = [F for F in conn_subfamilies(H, m - 1)]
    res = [reducible(adj, H, F) for F in conn]
    whole = reducible(adj, H, tuple(range(m)))
    print("m", m, "proper connected subfamilies", len(conn), "reducible:", sum(1 for r in res if r), "not pair-closed:", sum(1 for r in res if r is None), "whole reducible:", whole)
