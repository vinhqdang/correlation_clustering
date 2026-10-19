"""Tests of the lemmas in notes/fpt_attempt.md (brute force, n <= 7).
 T1  Theorem A: bucket-elimination DP over x in {nonempty subsets}^H, U enumerated
     (<= l unpacked pairs) agrees with brute-force OPT <= |H|+l.
 T2  Lemma 2 (viability): in every solution of cost <= |H|+l, each unpacked edit uv
     has |D(uv)| <= |U|-1.
 T3  Lemma 3 (maximality): extending H to maximal H' and l -> l-(|H'|-|H|) is equivalent;
     for maximal H each induced P3 has >= 2 vertices in V(H); non-adjacent pairs inside N(r),
     r outside V(H), are packed pairs.
 T4  Lemma 4 (locality of the interaction): a triple containing an unpacked edited pair uv and
     a packed pair has its third vertex z in Q(u) u Q(v).
 T5  search: forced unpacked edits with an endpoint outside V(H), and forced edits with both
     endpoints outside V(H) (locality counterexample data).
"""
import itertools, random, sys
from check_csp_view import pairs, induced_p3s, is_cluster, partitions

def fs(*a): return frozenset(a)
def rand_graph(rng, n):
    p = rng.choice([0.25, 0.4, 0.55])
    return {fs(a, b) for a, b in itertools.combinations(range(n), 2) if rng.random() < p}
def rand_packing(rng, n, E, maximal):
    cand = induced_p3s(n, E); rng.shuffle(cand)
    used, H = set(), []
    for P in cand:
        ps = pairs(P)
        if not any(q in used for q in ps):
            if maximal or rng.random() < .6:
                H.append(P); used.update(ps)
    return H
def sol_sets(n, E, budget):
    """all edit sets of partitions with cost <= budget (as sets of frozenset pairs)"""
    out = set()
    for lab in partitions(n):
        S = frozenset(fs(a, b) for a, b in itertools.combinations(range(n), 2)
                      if (lab[a] == lab[b]) != (fs(a, b) in E))
        if len(S) <= budget: out.add(S)
    return out
def opt(n, E):
    best = 10**9
    for lab in partitions(n):
        c = sum((lab[a] == lab[b]) != (fs(a, b) in E) for a, b in itertools.combinations(range(n), 2))
        best = min(best, c)
    return best

# ---------- T1: bucket elimination on I(H) -------------
def dp_feasible(n, E, H, U, budget):
    """min over x in nonempty-subsets^H of sum |x_i| (+|U|) s.t. E xor S is cluster; via bucket elimination"""
    p = len(H); owner = {}
    for i, P in enumerate(H):
        for j, q in enumerate(pairs(P)): owner[q] = (i, j)
    F = []  # factors: (scope, table dict)
    for i in range(p): F.append(((i,), {(m,): bin(m).count('1') for m in range(1, 8)}))
    for tri in itertools.combinations(range(n), 3):
        ps = [fs(a, b) for a, b in itertools.combinations(tri, 2)]
        sc = sorted({owner[q][0] for q in ps if q in owner})
        tab = {}
        for vals in itertools.product(range(1, 8), repeat=len(sc)):
            asg = dict(zip(sc, vals)); e = []
            for q in ps:
                adj = q in E
                if q in owner and (asg[owner[q][0]] >> owner[q][1]) & 1: adj = not adj
                if q in U: adj = not adj
                e.append(adj)
            tab[vals] = 0 if sum(e) != 2 else None
        F.append((tuple(sc), tab))
    # eliminate with min-degree order
    vars_ = set(range(p)); INF = 10**9
    while vars_:
        def deg(v): return len({u for s, _ in F if v in s for u in s} - {v})
        v = min(vars_, key=deg); vars_.remove(v)
        rel = [f for f in F if v in f[0]]; F = [f for f in F if v not in f[0]]
        sc = tuple(sorted({u for s, _ in rel for u in s} - {v}))
        tab = {}
        for vals in itertools.product(range(1, 8), repeat=len(sc)):
            asg = dict(zip(sc, vals)); best = INF
            for m in range(1, 8):
                asg[v] = m; c = 0
                for s, t in rel:
                    r = t[tuple(asg[u] for u in s)]
                    if r is None: c = INF; break
                    c += r
                best = min(best, c)
            tab[vals] = None if best >= INF else best
        F.append((sc, tab))
    tot = 0
    for s, t in F:
        r = t[()]
        if r is None: return None
        tot += r
    return tot + len(U)

def T1(trials, seed):
    rng = random.Random(seed); cnt = 0
    for _ in range(trials):
        n = rng.randint(4, 6); E = rand_graph(rng, n); H = rand_packing(rng, n, E, False)
        if not H or len(H) > 3: continue
        o = opt(n, E); unp = [fs(a, b) for a, b in itertools.combinations(range(n), 2)
                              if all(fs(a, b) not in pairs(P) for P in H)]
        for l in range(0, 3):
            bf = o <= len(H) + l
            ok = False
            for k in range(0, l + 1):
                for U in itertools.combinations(unp, k):
                    r = dp_feasible(n, E, H, set(U), len(H) + l)
                    if r is not None and r <= len(H) + l: ok = True; break
                if ok: break
            assert ok == bf, ("T1", n, sorted(map(sorted, E)), H, l)
        cnt += 1
    print("T1 ok on", cnt)

# ---------- T2, T4 -------------
def Q(H, u):
    return {x for P in H if u in P for x in P if x != u}
def T2_T4(trials, seed):
    rng = random.Random(seed); cnt = 0
    for _ in range(trials):
        n = rng.randint(4, 7); E = rand_graph(rng, n); H = rand_packing(rng, n, E, rng.random() < .5)
        if not H: continue
        packed = {q for P in H for q in pairs(P)}
        o = opt(n, E)
        for l in range(0, 3):
            if o > len(H) + l: continue
            for S in sol_sets(n, E, len(H) + l):
                U = [q for q in S if q not in packed]
                assert all(any(q in S for q in pairs(P)) for P in H)
                for uv in U:
                    u, v = tuple(uv)
                    D = []
                    for z in range(n):
                        if z in (u, v) or z in Q(H, u) | Q(H, v): continue
                        Ef = E ^ {uv}
                        a, b, c = fs(u, v) in Ef, fs(u, z) in Ef, fs(v, z) in Ef
                        if a + b + c == 2: D.append(z)
                    assert len(D) <= len(U) - 1, ("T2", n, sorted(map(sorted, E)), H, S)
                    # T4
                    for z in range(n):
                        if z in (u, v): continue
                        if any(fs(*t) in packed for t in [(u, z), (v, z)]):
                            assert z in Q(H, u) | Q(H, v)
                cnt += 1
    print("T2/T4 ok on", cnt, "solutions")

# ---------- T3 -------------
def T3(trials, seed):
    rng = random.Random(seed); cnt = 0
    for _ in range(trials):
        n = rng.randint(4, 7); E = rand_graph(rng, n); H0 = rand_packing(rng, n, E, False)
        if not H0: continue
        H1 = list(H0); used = {q for P in H1 for q in pairs(P)}
        for P in induced_p3s(n, E):
            if not any(q in used for q in pairs(P)): H1.append(P); used.update(pairs(P))
        o = opt(n, E); t = len(H1) - len(H0)
        for l in range(0, 4):
            assert (o <= len(H0) + l) == (o <= len(H1) + l - t)
        Vh = {x for P in H1 for x in P}
        for P in induced_p3s(n, E):
            assert sum(x in Vh for x in P) >= 2 and any(q in used for q in pairs(P))
        for r in set(range(n)) - Vh:
            N = [x for x in range(n) if fs(r, x) in E]
            for a, b in itertools.combinations(N, 2):
                if fs(a, b) not in E: assert fs(a, b) in used
        cnt += 1
    print("T3 ok on", cnt)

# ---------- T5 -------------
def T5(trials, seed):
    rng = random.Random(seed); forcedR = forcedRR = 0; exR = exRR = None
    for _ in range(trials):
        n = rng.randint(5, 7); E = rand_graph(rng, n); H = rand_packing(rng, n, E, True)
        if not H: continue
        packed = {q for P in H for q in pairs(P)}; Vh = {x for P in H for x in P}
        o = opt(n, E)
        if o <= len(H): continue
        sols = [S for S in sol_sets(n, E, o) if len(S) == o]
        def cnt(S, f): return sum(1 for q in S if q not in packed and f(q))
        if all(cnt(S, lambda q: not q <= Vh) > 0 for S in sols):
            forcedR += 1; exR = exR or (n, sorted(map(sorted, E)), H)
        if all(cnt(S, lambda q: not (q & Vh)) > 0 for S in sols):
            forcedRR += 1; exRR = exRR or (n, sorted(map(sorted, E)), H)
    print("T5: instances where every optimum has an unpacked edit with an endpoint outside V(H):", forcedR, exR)
    print("T5: ... with BOTH endpoints outside V(H):", forcedRR, exRR)

if __name__ == "__main__":
    T1(150, 1); T2_T4(150, 2); T3(400, 3); T5(400, 4)
