"""Check: (T, N[X_t]) is a tree decomposition of I(H)^2 with bag size <= (w+1)(3w+1), w = width of the given TD of I(H)
(elimination-order TDs, so w >= tw; the claim is for any TD).  Also max degree <= 3(c-1) and c <= w+1.
Members = random 3-subsets of [N] with pairwise <=1 common vertex (pair-disjoint)."""
import random, itertools
def rand_members(rng, N, m):
    H = []
    for _ in range(m * 5):
        t = frozenset(rng.sample(range(N), 3))
        if all(len(t & s) <= 1 for s in H): H.append(t)
        if len(H) == m: break
    return H
def elim_td(adj):
    """min-degree elimination -> (bags list, tree edges)"""
    nodes = set(adj); g = {v: set(adj[v]) for v in adj}; order = []; bags = {}
    while g:
        v = min(g, key=lambda x: (len(g[x]), x)); nb = set(g[v]); bags[v] = {v} | nb; order.append(v)
        for a in nb: g[a] |= nb - {a}; g[a].discard(v)
        del g[v]
    pos = {v: i for i, v in enumerate(order)}; edges = []
    for v in order:
        rest = bags[v] - {v}
        if rest: edges.append((v, min(rest, key=lambda x: pos[x])))
    return bags, edges
def is_td(adj, bags, edges):
    # tree-ness not rechecked (elimination tree is a forest; add dummy links irrelevant for the running-intersection test per component)
    import collections
    tadj = collections.defaultdict(set)
    for a, b in edges: tadj[a].add(b); tadj[b].add(a)
    for v in adj:
        S = {t for t in bags if v in bags[t]}
        if not S: return False
        seen = {next(iter(S))}; st = list(seen)
        while st:
            x = st.pop()
            for y in tadj[x]:
                if y in S and y not in seen: seen.add(y); st.append(y)
        if seen != S: return False
    for u in adj:
        for v in adj[u]:
            if not any(u in b and v in b for b in bags.values()): return False
    return True
rng = random.Random(5); worst = 0; checked = 0; tight = 0
for it in range(4000):
    N = rng.randint(6, 16); m = rng.randint(2, 9)
    H = rand_members(rng, N, m)
    if len(H) < 2: continue
    adj = {i: {j for j in range(len(H)) if j != i and H[i] & H[j]} for i in range(len(H))}
    bags, edges = elim_td(adj); w = max(len(b) for b in bags.values()) - 1
    c = max(sum(1 for t in H if v in t) for v in range(N))
    D = max(len(a) for a in adj.values())
    assert c <= w + 1 or w == 0 and c <= 1 or True
    if D > 3 * (c - 1): print('DEG VIOLATION', H)
    if c > w + 1: print('c>w+1', H, w, c)
    sq = {i: set() for i in adj}
    for i in adj:
        for j in adj[i]: sq[i].add(j); sq[i] |= (adj[j] - {i})
    nb = {t: set().union(*[({x} | adj[x]) for x in bags[t]]) for t in bags}
    if not is_td(sq, nb, edges): print('NOT A TD', H); continue
    mx = max(len(b) for b in nb.values()); bound = (w + 1) * (3 * w + 1)
    if mx > bound: print('SIZE VIOLATION', mx, bound)
    worst = max(worst, mx / bound); checked += 1; tight += mx == bound
print('checked', checked, 'max size/bound', worst, 'tight', tight)
