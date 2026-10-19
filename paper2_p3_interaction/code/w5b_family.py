"""Wave 5 B: 'far default links' family.  F = cubic graph; every edge e=(u->v) of F becomes a P3 (a_e,b_e,c_e),
a_e in clique K_u, b_e,c_e in clique K_v, one cross edge a_e b_e.  Claims checked by computer:
 H = {P3s} is a pair-disjoint induced packing (vertex-disjoint, so w = tw(I(H)) = 0, c = 1), pair-maximal,
 the partition {K_v} costs |H| so OPT = |H| (l = 0, YES), all hard degrees are <= 5, and the graph
 Gamma = soft pairs + hard edges has F as a minor (contract the connected hard components K_v), so tw(Gamma) >= tw(F)."""
import sys, itertools, random
import networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in
sys.path.insert(0, '.')

def mmd_plus(G):
    """minor-min-width lower bound on treewidth (contract a min-degree vertex with a min-degree neighbour)"""
    G = nx.Graph(G); lb = 0
    while G.number_of_nodes() > 1:
        v = min(G.nodes, key=G.degree); lb = max(lb, G.degree(v))
        if G.degree(v) == 0: G.remove_node(v); continue
        u = min(G.neighbors(v), key=G.degree); G = nx.contracted_nodes(G, u, v, self_loops=False)
    return lb

def build(F):
    idx = {}; P3 = []; node_members = {v: [] for v in F}
    cnt = 0
    for (u, v) in F.edges():
        a, b, c = cnt, cnt + 1, cnt + 2; cnt += 3; P3.append((a, b, c))
        node_members[u].append(a); node_members[v] += [b, c]
    n = cnt; E = set()
    for v, mem in node_members.items():
        for x, y in itertools.combinations(mem, 2): E.add((min(x, y), max(x, y)))
    for a, b, c in P3: E.add((min(a, b), max(a, b)))
    return n, E, P3, node_members

def check(F, name):
    n, E, P3, K = build(F)
    adj = lambda x, y: (min(x, y), max(x, y)) in E
    soft = set()
    for t in P3:
        for x, y in itertools.combinations(t, 2): soft.add((min(x, y), max(x, y)))
    ok_ind = all(adj(a, b) and adj(b, c) and not adj(a, c) for a, b, c in P3)
    ok_disj = len(soft) == 3 * len(P3)
    # pair-maximal: every induced P3 of G has a soft pair
    pm = True
    for b in range(n):
        nb = [x for x in range(n) if x != b and adj(x, b)]
        for a, c in itertools.combinations(nb, 2):
            if not adj(a, c):
                ps = [(min(a, b), max(a, b)), (min(b, c), max(b, c)), (min(a, c), max(a, c))]
                if not any(p in soft for p in ps): pm = False
    # BFK Rule 6.2 conditions for every packed P3 with S = {ab}: clusters {a},{b,c}
    nbr = lambda x: set(y for y in range(n) if y != x and adj(x, y))
    rule = True
    for a, b, c in P3:
        W = {a, b, c}
        if (nbr(b) - W) != (nbr(c) - W): rule = False
        if (nbr(a) & nbr(b)) - W or (nbr(a) & nbr(c)) - W: rule = False
    print('   BFK Rule 6.2 (S={ab}) applicable to every P3:', rule)
    lab = {}
    for v, mem in K.items():
        for x in mem: lab[x] = v
    cost = sum(1 for x, y in itertools.combinations(range(n), 2) if (lab[x] == lab[y]) != adj(x, y))
    hard_deg = [sum(1 for y in range(n) if y != x and adj(x, y) and (min(x, y), max(x, y)) not in soft) for x in range(n)]
    Gam = nx.Graph(); Gam.add_nodes_from(range(n))
    Gam.add_edges_from(soft); Gam.add_edges_from(p for p in E if p not in soft)
    ub = treewidth_min_fill_in(Gam)[0]
    print(f"{name}: |F|={F.number_of_nodes()} nodes, n={n}, |H|={len(P3)}, induced={ok_ind}, pair-disjoint={ok_disj}, pair-maximal={pm}, "
          f"cost(K_v partition)={cost}, maxhard_deg={max(hard_deg)}, tw(F)>={mmd_plus(F)}, tw(Gamma)<={ub}, tw(Gamma)>=(minor F)")
    return ok_ind and ok_disj and pm and cost == len(P3)

if __name__ == '__main__':
    ok = True
    ok &= check(nx.complete_graph(4), 'K4')
    ok &= check(nx.complete_bipartite_graph(3, 3), 'K33')
    ok &= check(nx.petersen_graph(), 'Petersen')
    for m, sd in ((12, 1), (20, 2), (30, 3), (40, 4)):
        ok &= check(nx.random_regular_graph(3, m, seed=sd), f'cubic{m}')
    print('ALL OK' if ok else 'FAIL')
