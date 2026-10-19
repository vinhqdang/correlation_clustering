"""Family N(F): 'clique-necklace with expander leaves'.  F = digraph (edges oriented), cyclic edge order sigma.
Vertices per edge e: a_e,b_e,c_e,u_e.  Cliques K_v+ = {a_e: tail(e)=v} u {b_e,c_e,u_e: head(e)=v}.
Members P_e=(a_e,b_e,c_e), X_e=(u_e,b_e,u_sigma(e)).  Edges: all pairs inside cliques, a_e b_e, b_e u_sigma(e).
Checks: induced+pair-disjoint, pair-maximal, cost({K_v+})=|H| (so OPT=|H|, l=0), every connected pair-closed subfamily of size<=smax
is Rule-6.2 irreducible (no tight S_F with R), whole component reducible, hard degrees, quotient by cliques contains F as minor."""
import sys, itertools, random
from w6a_lib import *
from w6a_rule import *
def build(edges):   # edges: list of (tail, head) in cyclic order
    N = len(edges); idx = {}; nv = 0
    for e in range(N):
        for nm in "abcu": idx[(nm, e)] = nv; nv += 1
    sigma = []
    for e in range(N):
        for k in range(1, N):
            f = (e + k) % N
            if edges[f][1] != edges[e][1]: sigma.append(f); break
        else: return None
    verts = {}
    for e, (t, h) in enumerate(edges):
        verts.setdefault(t, []).append(idx[("a", e)])
        for nm in "bcu": verts.setdefault(h, []).append(idx[(nm, e)])
    E = set()
    for v, vs in verts.items():
        for x, y in itertools.combinations(vs, 2): E.add(key(x, y))
    H = []
    for e in range(N):
        a, b, c, u = (idx[(nm, e)] for nm in "abcu")
        E.add(key(a, b)); E.add(key(b, idx[("u", sigma[e])]))
        H.append((a, b, c)); H.append((u, b, idx[("u", sigma[e])]))
    return nv, from_edges(nv, E), H, list(verts.values()), sigma
def check(edges, smax):
    r = build(edges)
    if r is None: return None
    n, adj, H, cl, sigma = r
    ok_ind = all((adj[x] >> y & 1) and (adj[y] >> z & 1) and not (adj[x] >> z & 1) for x, y, z in H)
    pd = len(soft_pairs(H)) == 3 * len(H)
    pm = is_pair_maximal(adj, H)
    # cost of clique partition
    lab = {}
    for i, c in enumerate(cl):
        for v in c: lab[v] = i
    cost = sum(1 for x in range(n) for y in range(x + 1, n) if (adj[x] >> y & 1) != (lab[x] == lab[y]))
    hd = max(hard_deg(adj, H)); comp, nc = hard_components(adj, H)
    sizes = max(sum(1 for v in range(n) if comp[v] == c) for c in range(nc))
    red = []
    for F in conn_subfamilies(H, smax):
        W = members_vertices(H, F)
        if not pair_closed(H, F, W): continue
        if any(R_ok(adj, S, W) for S in tight_choices(adj, H, F, W)): red.append(F)
    T = list(range(len(H))); Wt = members_vertices(H, T)
    whole = any(R_ok(adj, S, Wt) for S in tight_choices(adj, H, T, Wt)) if len(H) <= 14 else None
    return dict(n=n, members=len(H), induced=ok_ind, pairdisj=pd, pairmax=pm, cost=cost, hardc=2 * len(edges), maxhd=hd, maxcomp=sizes, reducible_small=red[:3], nred=len(red), whole=whole, w_I=None)
if __name__ == "__main__":
    # necklace-like small: cycle-necklace of m cliques needs F with a cyclic structure: use K4 (6 edges), K33 (9 edges)
    K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    K33 = [(i, 3 + j) for i in range(3) for j in range(3)]
    tri = [(0, 1), (1, 2), (2, 0)]
    for name, ed in [("triangle", tri), ("K4", K4), ("K33", K33)]:
        rng = random.Random(1)
        for t in range(3):
            e2 = ed[:]; rng.shuffle(e2)
            if t == 0: e2 = ed[:]
            print(name, t, check(e2, 3 if len(ed) < 8 else 2), flush=True)
