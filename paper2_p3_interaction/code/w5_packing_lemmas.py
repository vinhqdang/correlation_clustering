"""Wave 5: brute-force test of the statements of Section 'Structure of solutions with respect to a packing'.
Own code (imports nothing from the other scripts).  Exhaustive over all set partitions for n <= 8.
Statements tested: MAX (a)-(d) + tight star, HE, HUB+ (nu <= l, alpha <= 2l+c+1, corollary |X| <= 2l+1),
SM (i), (ii), (iii sharp), pieces inside clusters, soft-graph treewidth bound, D = 2c.
Usage: python3 w5_packing_lemmas.py [n_graphs] [seed]
"""
import itertools, random, sys
from functools import lru_cache

def partitions(n):
    def rec(i, labels, k):
        if i == n:
            yield tuple(labels); return
        for c in range(k + 1):
            labels.append(c)
            yield from rec(i + 1, labels, max(k, c + 1))
            labels.pop()
    yield from rec(0, [], 0)

def popcount(x): return bin(x).count("1")

class Ctx:
    def __init__(self, n):
        self.n = n
        self.pairs = list(itertools.combinations(range(n), 2))
        self.idx = {p: i for i, p in enumerate(self.pairs)}
        self.parts = []
        for lab in partitions(n):
            m = 0
            for i, (a, b) in enumerate(self.pairs):
                if lab[a] == lab[b]: m |= 1 << i
            self.parts.append((lab, m))
    def pid(self, a, b): return self.idx[(a, b) if a < b else (b, a)]

def treewidth(nodes, adj):
    """exact treewidth by subset DP (nodes: list, adj: dict node->set). returns 0 for edgeless/empty."""
    nodes = list(nodes); k = len(nodes)
    if k == 0: return 0
    pos = {v: i for i, v in enumerate(nodes)}
    nb = [0] * k
    for v in nodes:
        for u in adj[v]: nb[pos[v]] |= 1 << pos[u]
    full = (1 << k) - 1
    @lru_cache(None)
    def q(S, v):  # |{w not in S, w != v : path v..w through S}|
        seen = 1 << v; stack = [v]; reach = 0
        while stack:
            x = stack.pop()
            for y in range(k):
                if nb[x] >> y & 1 and not seen >> y & 1:
                    seen |= 1 << y
                    if S >> y & 1: stack.append(y)
                    else: reach |= 1 << y
        return popcount(reach)
    dp = {0: -1}
    for S in range(1, full + 1):
        best = 10 ** 9
        for v in range(k):
            if S >> v & 1:
                T = S ^ (1 << v)
                best = min(best, max(dp[T], q(T, v)))
        dp[S] = best
    return max(dp[full], 0)

def induced_p3s(n, E):
    out = []
    for b in range(n):
        for a, c in itertools.combinations([x for x in range(n) if x != b], 2):
            if (min(a, b), max(a, b)) in E and (min(c, b), max(c, b)) in E and (a, c) not in E:
                out.append((a, b, c))
    return out

def pairs_of(P):
    a, b, c = P
    return [tuple(sorted(x)) for x in ((a, b), (b, c), (a, c))]

def random_graph(n, rng):
    t = rng.random()
    E = set()
    if t < 0.4:
        p = rng.choice([0.2, 0.35, 0.5, 0.7])
        for pr in itertools.combinations(range(n), 2):
            if rng.random() < p: E.add(pr)
    elif t < 0.8:
        k = rng.randint(1, 3); lab = [rng.randrange(k) for _ in range(n)]
        for pr in itertools.combinations(range(n), 2):
            e = lab[pr[0]] == lab[pr[1]]
            if rng.random() < 0.15: e = not e
            if e: E.add(pr)
    else:  # hub / star-like
        h = rng.randrange(n)
        for v in range(n):
            if v != h and rng.random() < 0.8: E.add(tuple(sorted((h, v))))
        for pr in itertools.combinations(range(n), 2):
            if h not in pr and rng.random() < 0.1: E.add(pr)
    return E

def greedy_packing(n, E, used, rng, stop=1.0):
    """random greedy pair-disjoint packing extending nothing; 'used' = set of pairs already packed (mutated)"""
    P3 = induced_p3s(n, E); rng.shuffle(P3); res = []
    for P in P3:
        if rng.random() > stop: continue
        ps = pairs_of(P)
        if all(p not in used for p in ps):
            res.append(P); used.update(ps)
    return res

def inter_graph(H):
    adj = {i: set() for i in range(len(H))}
    for i, j in itertools.combinations(range(len(H)), 2):
        if set(H[i]) & set(H[j]): adj[i].add(j); adj[j].add(i)
    return adj

def max_cliques_free_alpha(vs, E):
    """independence number of G[vs] by brute force"""
    vs = list(vs); best = 0
    for r in range(len(vs), 0, -1):
        for sub in itertools.combinations(vs, r):
            if all(tuple(sorted(p)) not in E for p in itertools.combinations(sub, 2)):
                return r
    return 0

def max_matching(vs, edges):
    vs = list(vs); edges = [e for e in edges]
    @lru_cache(None)
    def rec(i, usedmask):
        if i == len(edges): return 0
        a, b = edges[i]
        best = rec(i + 1, usedmask)
        if not (usedmask >> a & 1) and not (usedmask >> b & 1):
            best = max(best, 1 + rec(i + 1, usedmask | 1 << a | 1 << b))
        return best
    return rec(0, 0)

def check(ctx, E, H, stats, tag, ell_extra=3):
    n = ctx.n
    Emask = 0
    for p in E: Emask |= 1 << ctx.idx[p]
    packed = set()
    for P in H: packed.update(pairs_of(P))
    pmask = 0
    for p in packed: pmask |= 1 << ctx.idx[p]
    hsz = len(H)
    c = max([sum(1 for P in H if v in P) for v in range(n)] + [0])
    Q = {v: set() for v in range(n)}
    for (a, b) in packed: Q[a].add(b); Q[b].add(a)
    D = max(len(Q[v]) for v in range(n))
    assert D == 2 * c, ("D=2c", D, c)
    # costs
    best = 10 ** 9
    sols = []
    for lab, same in ctx.parts:
        S = same ^ Emask
        cst = popcount(S)
        if cst < best: best = cst
        if cst <= hsz + ell_extra: sols.append((lab, S, cst))
    OPT = best
    assert OPT >= hsz, "OPT>=|H|"
    # ---- soft graph treewidth and tw(I(H))
    w = treewidth(range(len(H)), inter_graph(H)) if H else 0
    Wadj = {v: set(Q[v]) for v in range(n)}
    twW = treewidth(range(n), Wadj)
    assert twW <= (3 * (w + 1) - 1 if H else 0), ("tw(W)", twW, w)
    stats["tw(W)"] += 1
    # ---- HE (pair-maximal H only)
    pairmax = all(any(p in packed for p in pairs_of(P)) for P in induced_p3s(n, E))
    if pairmax:
        for b in range(n):
            hardnb = [a for a in range(n) if a != b and tuple(sorted((a, b))) in E and tuple(sorted((a, b))) not in packed]
            for a, cc in itertools.combinations(hardnb, 2):
                assert tuple(sorted((a, cc))) in E or tuple(sorted((a, cc))) in packed, "HE"
                stats["HE"] += 1
    # ---- HUB+ : every vertex u, every l >= OPT-|H| (l up to ell_extra)
    for u in range(n):
        Nh = [x for x in range(n) if x != u and tuple(sorted((u, x))) in E and tuple(sorted((u, x))) not in packed]
        gam = [(x, y) for x, y in itertools.combinations(Nh, 2)
               if (x, y) not in E and (x, y) not in packed]
        nu = max_matching(Nh, tuple(gam)) if False else None
        # matching on indices
        ix = {x: i for i, x in enumerate(Nh)}
        nu = max_matching(range(len(Nh)), tuple((ix[x], ix[y]) for x, y in gam))
        al = max_cliques_free_alpha(Nh, E)
        # clique of Gamma (hard non-adjacent pairwise): largest
        omega = 0
        for r in range(len(Nh), 0, -1):
            if any(all(pr in set(gam) or (pr[1], pr[0]) in set(gam) for pr in itertools.combinations(sub, 2))
                   for sub in itertools.combinations(Nh, r)):
                omega = r; break
        for l in range(max(0, OPT - hsz), ell_extra + 1):
            assert nu <= l, ("HUB nu", nu, l)
            assert 2 * nu <= 2 * l
            assert al <= 2 * l + c + 1, ("HUB alpha", al, l, c)
            assert omega <= 2 * l + 1, ("HUB(a)", omega, l)
            stats["HUB"] += 1
    # ---- SM on every solution with cost <= |H| + ell_extra
    for lab, S, cst in sols:
        U = S & ~pmask
        u = popcount(U)
        assert u <= cst - hsz, "|U|<=|S|-|H|"
        F = U & Emask
        I = U & ~Emask
        # kept hard edges
        par = list(range(n))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for (a, b) in E:
            i = ctx.idx[(a, b)]
            if (a, b) not in packed and not (F >> i & 1):
                par[find(a)] = find(b)
        piece = [find(v) for v in range(n)]
        pieces = {}
        for v in range(n): pieces.setdefault(piece[v], []).append(v)
        clusters = {}
        for v in range(n): clusters.setdefault(lab[v], []).append(v)
        for r_, vs in pieces.items():
            assert len({lab[v] for v in vs}) == 1, "piece inside cluster"
        for cl in clusters.values():
            ps = {}
            for v in cl: ps.setdefault(piece[v], []).append(v)
            plist = list(ps.values())
            if len(plist) < 2: continue
            # (i)
            for A, B in itertools.combinations(plist, 2):
                nI = 0
                for a in A:
                    for b in B:
                        pr = (min(a, b), max(a, b))
                        assert pr in packed or (I >> ctx.idx[pr] & 1), "SM(i)"
                        if I >> ctx.idx[pr] & 1: nI += 1
                a_, b_ = max(len(A), len(B)), min(len(A), len(B))
                assert b_ * (a_ - D) <= nI <= u, "SM(ii) product"
                assert a_ <= D + u, "SM(ii) size"
                stats["SM(ii)"] += 1
            bound = max(D + u + 1, 2 * D, 2 * u)
            assert len(cl) <= bound, ("SM(iii)", len(cl), D, u)
            loose = max(2 * D + (cst - hsz), 2 * (cst - hsz))
            assert len(cl) <= loose, "loose bound implied"
            stats["SM(iii)"] += 1
            stats["SM(iii)_slack_min"] = min(stats.get("SM(iii)_slack_min", 99), bound - len(cl))
    return OPT

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rng = random.Random(seed)
    stats = {k: 0 for k in ["graphs", "MAX", "tw(W)", "HE", "HUB", "SM(ii)", "SM(iii)", "MAXtw_eq"]}
    ctxs = {n: Ctx(n) for n in (5, 6, 7, 8)}
    for it in range(N):
        n = rng.choice([5, 6, 7, 7, 8]) if it % 3 else rng.choice([6, 7])
        ctx = ctxs[n]
        E = random_graph(n, rng)
        stats["graphs"] += 1
        # H: partial random packing (any packing)
        used = set()
        H = greedy_packing(n, E, used, rng, stop=rng.choice([0.3, 0.6, 1.0]))
        OPT = check(ctx, E, H, stats, "H")
        # MAX
        used2 = set(used)
        M = greedy_packing(n, E, used2, rng, 1.0)  # maximal: every P3 with all pairs unused was tried
        HM = H + M
        allp = set()
        for P in HM: allp.update(pairs_of(P))
        assert len(allp) == 3 * len(HM), "H u M pair-disjoint"
        assert OPT >= len(H) + len(M), "MAX (a)"
        assert all(any(p in allp for p in pairs_of(P)) for P in induced_p3s(n, E)), "MAX (c) pair-maximal"
        wH = treewidth(range(len(H)), inter_graph(H)) if H else 0
        wHM = treewidth(range(len(HM)), inter_graph(HM)) if HM else 0
        assert wHM <= wH + len(M), ("MAX (d)", wH, wHM, len(M))
        cH = max([sum(1 for P in H if v in P) for v in range(n)] + [0])
        cHM = max([sum(1 for P in HM if v in P) for v in range(n)] + [0])
        assert cHM <= cH + len(M)
        assert cHM <= wHM + 1 or not HM
        if wHM == wH + len(M) and M: stats["MAXtw_eq"] += 1
        for l in range(0, 4):
            yes1 = OPT <= len(H) + l
            yes2 = len(M) <= l and OPT <= len(HM) + (l - len(M))
            assert yes1 == yes2, "MAX (b)"
        stats["MAX"] += 1
        # all statements again for H u M (primed instance), incl. HE since pair-maximal
        check(ctx, E, HM, stats, "HM")
    # tight star example: p,q,a_i,b_i around centre v
    for m in range(1, 6):
        n = 3 + 2 * m
        v = 0
        E = {(0, x) for x in range(1, n)}
        H = [(1, 0, 2)]
        M = [(2 * i + 1, 0, 2 * i + 2) for i in range(1, m + 1)]
        HM = H + M
        used = set()
        for P in HM: used.update(pairs_of(P))
        assert all(any(p in used for p in pairs_of(P)) for P in induced_p3s(n, E)), "star pair-maximal"
        assert treewidth(range(1), inter_graph(H)) == 0
        assert treewidth(range(len(HM)), inter_graph(HM)) == m
    print("star example: tw(I(H)) = 0, tw(I(H u M)) = m verified for m = 1..5")
    print(stats)
    print("ALL OK")

main()
