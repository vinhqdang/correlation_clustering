"""a13_t7_L.py: Theorem L.  (a) explicit instance where the triple scope needs radius 3 in I(H) (6-path of members), (b) random radius check,
(c) width of the primal graph of the leaf-test CSP (variables: members, local-pair bits; one constraint per triple with a non-far pair) on G_N for growing N."""
import itertools, random, networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in
from a13_lib import *
def dist_IH(inst):
    adj = IH(inst); G = nx.Graph(); G.add_nodes_from(adj)
    for i in adj:
        for j in adj[i]: G.add_edge(i, j)
    return dict(nx.all_pairs_shortest_path_length(G)), G
def scope_radius(inst, D, cnt):
    n, mem = inst['n'], inst['members']
    M = [[i for i, m in enumerate(mem) if u in m] for u in range(n)]
    worst = 0; rows = []
    for tr in itertools.combinations(range(n), 3):
        ps = [frozenset(p) for p in itertools.combinations(tr, 2)]
        nf = [p for p in ps if p in cnt]
        if not nf: continue
        W = set()
        for p in nf:
            for u in p: W |= set(M[u])
        r = min(max(D[c].get(w, 99) for w in W) for c in range(len(mem)))
        worst = max(worst, r)
    return worst
# (a) explicit 6-path: members A',A,B,B',C,C' ; shared vertices u(A',A) z1(A,B) v(B,B') z2(B',C) w(C,C')
names = {}
def V(s): 
    names.setdefault(s, len(names)); return names[s]
mem = [(V('x1'), V('u'), V('y1')), (V('u'), V('a'), V('z1')), (V('z1'), V('b'), V('v')), (V('v'), V('d'), V('z2')), (V('z2'), V('e'), V('w')), (V('w'), V('f'), V('y2'))]
# member (end,centre,end): chain shares an END vertex each time; edges end-centre-end
E = set()
for (a, b, c) in mem: E.add(frozenset((a, b))); E.add(frozenset((b, c)))
inst = dict(n=len(names), E=E, members=mem, T=[])
print('six-path instance valid:', validate(inst))
D, G = dist_IH(inst); cnt = counted_pairs(inst)
u, v, w = V('u'), V('v'), V('w')
print('pairs uv,vw local:', frozenset((u, v)) in cnt, frozenset((v, w)) in cnt, ' uw local:', frozenset((u, w)) in cnt)
Wm = set(i for i, m in enumerate(mem) if (set(m) & {u, v, w}))
print('W members', sorted(Wm), 'radius needed =', min(max(D[c][x] for x in Wm) for c in range(6)), ' eccentricity profile', [max(D[c][x] for x in Wm) for c in range(6)])
# (b) random check of the lemma
random.seed(3); worst = 0
for it in range(40):
    m = random.randint(5, 8); verts = []; mem = []
    nxt = [0]
    def newv(): nxt[0] += 1; return nxt[0]-1
    # random tree-like packing: each new member attaches to a random old member by sharing one vertex
    mem.append((newv(), newv(), newv()))
    for k in range(m-1):
        old = random.choice(mem); sh = random.choice(old)
        others = [newv(), newv()]; tri = [sh] + others; random.shuffle(tri); mem.append(tuple(tri))
    E = set()
    for (a, b, c) in mem: E.add(frozenset((a, b))); E.add(frozenset((b, c)))
    inst = dict(n=nxt[0], E=E, members=mem, T=[])
    if nxt[0] > 17: continue
    # ensure members induced (add nothing); skip if not valid
    ok, c = validate(inst)
    if not ok: continue
    D, G = dist_IH(inst); worst = max(worst, scope_radius(inst, D, counted_pairs(inst)))
print('random tree-like instances: max radius needed =', worst)
# (c) CSP primal graph width on G_N
for N in (3, 6, 10, 14):
    inst = build_GN(N); n = inst['n']; cnt = counted_pairs(inst); soft = {}
    for i, m in enumerate(inst['members']):
        for p in pairs_of(m): soft[p] = i
    PG = nx.Graph()
    for tr in itertools.combinations(range(n), 3):
        ps = [frozenset(p) for p in itertools.combinations(tr, 2)]
        if not any(p in cnt for p in ps): continue
        vs = [('m', soft[p]) if p in soft else ('b', p) for p in ps if p in cnt]
        for a, b in itertools.combinations(vs, 2): PG.add_edge(a, b)
        for a in vs: PG.add_node(a)
    tw, _ = treewidth_min_fill_in(PG)
    print('G_%d: n=%d, CSP vars=%d, primal-graph width (min-fill) = %d' % (N, n, PG.number_of_nodes(), tw), flush=True)
