"""Path version N*_path(F): I(H) is a caterpillar (tree, w=1).  F' = F + pendant edge e*=(v->P).  Order e_0=e*, then the other edges of F in an order with
consecutive heads different; X_{e_i}=(u_i,b_i,u_{i+1}), last X=(u_last,b_last,z) with z a fresh vertex.  Trigger y hard-adjacent to u_{e*},c_{e*}.
Checks as in w6a_trig2 plus: I(H) acyclic."""
import sys, itertools, random
from w6a_lib import *
from w6a_rule import *
from w6a_trig import copt
from w6a_csp import tight_choices_dfs
def build_path(edges):
    ed = [(edges[0][0], "P")] + list(edges)
    # order the rest so consecutive heads differ
    rng = random.Random(7)
    for attempt in range(1000):
        rest = list(edges); rng.shuffle(rest); order = [ed[0]]
        ok = True
        while rest:
            for i, e in enumerate(rest):
                if e[1] != order[-1][1]: order.append(rest.pop(i)); break
            else: ok = False; break
        if ok: break
    else: return None
    names = {}
    ed2 = []
    for t, h in order:
        for x in (t, h):
            if x not in names: names[x] = len(names)
        ed2.append((names[t], names[h]))
    N = len(ed2); idx = {}; nv = 0
    for e in range(N):
        for nm in "abcu": idx[(nm, e)] = nv; nv += 1
    z = nv; nv += 1; y = nv; nv += 1
    verts = {}
    for e, (t, h) in enumerate(ed2):
        verts.setdefault(t, []).append(idx[("a", e)])
        for nm in "bcu": verts.setdefault(h, []).append(idx[(nm, e)])
    E = set()
    for v, vs in verts.items():
        for p, q in itertools.combinations(vs, 2): E.add(key(p, q))
    H = []
    for e in range(N):
        a, b, c, u = (idx[(nm, e)] for nm in "abcu")
        nxt = idx[("u", e + 1)] if e + 1 < N else z
        E.add(key(a, b)); E.add(key(b, nxt)); H.append((a, b, c)); H.append((u, b, nxt))
    E.add(key(y, idx[("u", 0)])); E.add(key(y, idx[("c", 0)]))
    return nv, from_edges(nv, E), H, list(verts.values()), y, idx, ed2
def is_forest_I(H):
    m = len(H); Ie = Iedges(H); par = list(range(m))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for i, j in Ie:
        a, b = f(i), f(j)
        if a == b: return False
        par[a] = b
    return True
def run(edges, smax, do_opt):
    r = build_path(edges)
    if r is None: return None
    n, adj, H, cl, y, idx, ed2 = r
    pm = is_pair_maximal(adj, H)
    ind = all((adj[x] >> yy & 1) and (adj[yy] >> zz & 1) and not (adj[x] >> zz & 1) for x, yy, zz in H)
    pdj = len(soft_pairs(H)) == 3 * len(H)
    lab = {}
    for i, c in enumerate(cl):
        for v in c: lab[v] = i
    z = n - 2
    lab[z] = len(cl); lab[y] = lab[idx[("u", 0)]]
    cost = sum(1 for a in range(n) for b in range(a + 1, n) if (adj[a] >> b & 1) != (lab[a] == lab[b]))
    T = list(range(len(H))); W = members_vertices(H, T)
    ch = tight_choices_dfs(adj, H, T, W); whole = any(R_ok(adj, S, W) for S in ch)
    red = 0
    for F in conn_subfamilies(H, smax):
        W2 = members_vertices(H, F)
        if pair_closed(H, F, W2) and any(R_ok(adj, S, W2) for S in tight_choices(adj, H, F, W2)): red += 1
    o = copt(adj) if do_opt else None
    comps = I_components(H)
    return dict(n=n, H=len(H), induced=ind, pairdisj=pdj, pairmax=pm, forest=is_forest_I(H), ncomp=len(comps), cost_explicit=cost, OPT=o, T_tight=bool(ch), T_whole_reducible=whole, nred_small=red, maxhd=max(hard_deg(adj, H)))
if __name__ == "__main__":
    tri = [(0, 1), (1, 2), (2, 0)]
    K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    print("triangle+pendant", run(tri, 4, True), flush=True)
    print("K4+pendant", run(K4, 3, False), flush=True)
    K33 = [(i, 3 + j) for i in range(3) for j in range(3)]
    print("K33+pendant", run(K33, 3, False), flush=True)
    # cubic graphs: cube, Petersen
    cube = [(a, a ^ (1 << b)) for a in range(8) for b in range(3) if a < a ^ (1 << b)]
    print("cube+pendant", run(cube, 3, False), flush=True)
    pet = [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    print("Petersen+pendant", run(pet, 2, False), flush=True)
