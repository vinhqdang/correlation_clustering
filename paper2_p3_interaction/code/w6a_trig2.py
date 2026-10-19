"""Family N*(F): N(F') for F' = F + pendant edge e*=(v->v*), plus one trigger vertex y hard-adjacent to u_{e*} and c_{e*} only.
Claims checked: pair-maximal; explicit partition costs |H|+1 (so YES with l=1); (small n) OPT = |H|+1; T=I(H)-component of all members is tight
but Rule-6.2 irreducible as a whole; no connected pair-closed subfamily of size <= smax is reducible; Gamma' contains F as a minor."""
import sys, itertools, random, subprocess
from w6a_lib import *
from w6a_rule import *
from w6a_fam import build
from w6a_trig import copt
def make(edges, k_star=0):
    # edges: orientation of F; add pendant edge from head of edges[k_star] ... use vertex 'P' new sink
    ed = list(edges) + [(edges[0][0], "P")]
    names = {}
    ed2 = []
    for t, h in ed:
        for x in (t, h):
            if x not in names: names[x] = len(names)
        ed2.append((names[t], names[h]))
    # ensure star edge is last; cyclic order as given
    n, adj, H, cl, sigma = build(ed2)
    N = len(ed2); estar = N - 1
    idx = lambda nm, e: 4 * e + "abcu".index(nm)
    y = n; n2 = n + 1
    adj = [a for a in adj] + [0]
    for v in (idx("u", estar), idx("c", estar)):
        adj[y] |= 1 << v; adj[v] |= 1 << y
    return n2, adj, H, cl + [[]], y, estar, ed2
def run(edges, smax, do_opt):
    n, adj, H, cl, y, estar, ed2 = make(edges)
    pm = is_pair_maximal(adj, H)
    # explicit partition: cliques; y joins clique of v*
    lab = {}
    for i, c in enumerate(cl[:-1]):
        for v in c: lab[v] = i
    idx = lambda nm, e: 4 * e + "abcu".index(nm)
    lab[y] = lab[idx("u", estar)]
    cost = sum(1 for a in range(n) for b in range(a + 1, n) if (adj[a] >> b & 1) != (lab[a] == lab[b]))
    T = list(range(len(H))); W = members_vertices(H, T)
    ch = tight_choices(adj, H, T, W)
    tight = bool(ch); whole_red = any(R_ok(adj, S, W) for S in ch)
    red = 0
    for F in conn_subfamilies(H, smax):
        W2 = members_vertices(H, F)
        if pair_closed(H, F, W2) and any(R_ok(adj, S, W2) for S in tight_choices(adj, H, F, W2)): red += 1
    o = copt(adj) if do_opt else None
    return dict(n=n, H=len(H), pairmax=pm, cost_explicit=cost, OPT=o, T_tight=tight, T_whole_reducible=whole_red, nred_small=red)
if __name__ == "__main__":
    tri = [(0, 1), (1, 2), (2, 0)]
    K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    print("triangle+pendant", run(tri, 3, True), flush=True)
    print("K4+pendant", run(K4, 2, False), flush=True)
