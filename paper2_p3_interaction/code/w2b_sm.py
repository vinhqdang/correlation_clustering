"""Wave-2 checks. (SM) small-merge lemma, (MAX) maximalisation lemma, (HE) hard-edge lemma.
Brute force over all set partitions, n <= 9.
SM:  for every solution S (cost <= |H|+l, l<=2), U = S minus packed pairs, F = deleted hard edges,
     I = inserted hard pairs, pieces = components of (E minus packed minus F).
     (i)  two vertices of one cluster in different pieces form a packed pair or a pair of I;
     (ii) for pieces A,B of one cluster: min(|A|,|B|)*(max(|A|,|B|)-Ds) <= |I cap AxB|, hence max<=Ds+|U|;
     (iii) a cluster with >=2 pieces has size <= max(2Ds+|U|+... ) -> checked as  s*(s-M-Ds) <= 2|U| with M=Ds+|U|.
MAX: greedy maximal extension M of all-hard induced P3s: |M|<=l if OPT<=|H|+l, and equivalence of instances.
HE:  if H is maximal, hard edges ab, bc imply ac in E or packed.
"""
import itertools, random, sys
from check_csp_view import pairs, induced_p3s, partitions
fs = frozenset
def gen(rng, n):
    # cluster graph with random clique sizes, then perturb a few pairs
    lab = [rng.randint(0, max(1, n // 3)) for _ in range(n)]
    E = {fs((a, b)) for a, b in itertools.combinations(range(n), 2) if lab[a] == lab[b]}
    for _ in range(rng.randint(1, 4)):
        E ^= {fs(rng.sample(range(n), 2))}
    return E
def rpack(rng, n, E, maximal=True, extra_ok=True):
    cand = induced_p3s(n, E); rng.shuffle(cand); used = set(); H = []
    for P in cand:
        if not any(q in used for q in pairs(P)) and (maximal or rng.random() < .6):
            H.append(P); used.update(pairs(P))
    return H
def cost_edits(lab, n, E):
    return {fs((a, b)) for a, b in itertools.combinations(range(n), 2) if (lab[a] == lab[b]) != (fs((a, b)) in E)}
def sm_test(trials, seed, nmax):
    rng = random.Random(seed); nsol = nmulti = ninst = 0; maxpiece = 0
    for _ in range(trials):
        n = rng.randint(5, nmax); E = gen(rng, n); H = rpack(rng, n, E, rng.random() < .5)
        if not H: continue
        packed = {q for P in H for q in pairs(P)}
        soft = {x: set() for x in range(n)}
        for q in packed:
            a, b = tuple(q); soft[a].add(b); soft[b].add(a)
        Ds = max(len(s) for s in soft.values())
        best = None; sols = []
        for lab in partitions(n):
            S = cost_edits(lab, n, E)
            if len(S) <= len(H) + 2: sols.append((lab, S))
        if not sols: continue
        ninst += 1
        for lab, S in sols:
            U = {q for q in S if q not in packed}
            assert len(U) <= len(S) - len(H)
            F = {q for q in U if q in E}; I = U - F
            # pieces
            par = list(range(n))
            def f(x):
                while par[x] != x: par[x] = par[par[x]]; x = par[x]
                return x
            for q in E:
                if q not in packed and q not in F:
                    a, b = tuple(q); par[f(a)] = f(b)
            piece = [f(x) for x in range(n)]
            nsol += 1
            for c in set(lab):
                C = [x for x in range(n) if lab[x] == c]
                ps = {}
                for x in C: ps.setdefault(piece[x], []).append(x)
                for x, y in itertools.combinations(C, 2):
                    if piece[x] != piece[y]:
                        assert fs((x, y)) in packed or fs((x, y)) in I, ("SM(i)", n, E, H, lab)
                lu = len(U)
                for (p1, A), (p2, B) in itertools.combinations(ps.items(), 2):
                    ins = sum(fs((a, b)) in I for a in A for b in B)
                    lo, hi = min(len(A), len(B)), max(len(A), len(B))
                    assert lo * (hi - Ds) <= ins, ("SM(ii)", n, E, H, lab)
                    assert hi <= Ds + lu
                if len(ps) >= 2:
                    nmulti += 1; s = len(C); M = Ds + lu
                    assert s * (s - M - Ds) <= 2 * lu, ("SM(iii)", n, E, H, lab)
                    maxpiece = max(maxpiece, max(len(v) for v in ps.values()))
    print(f"SM ok: instances {ninst}, solutions {nsol}, multi-piece clusters {nmulti}, max piece in a multi-piece cluster {maxpiece}")
def opt(n, E):
    return min(len(cost_edits(l, n, E)) for l in partitions(n))
def max_test(trials, seed):
    rng = random.Random(seed); cnt = viol = 0
    for _ in range(trials):
        n = rng.randint(4, 8); E = gen(rng, n); H = rpack(rng, n, E, False)
        packed = {q for P in H for q in pairs(P)}
        M = []; used = set(packed)
        for P in induced_p3s(n, E):
            if not any(q in used for q in pairs(P)): M.append(P); used.update(pairs(P))
        o = opt(n, E)
        for l in range(0, 4):
            lhs = o <= len(H) + l
            rhs = (len(M) <= l) and (o <= len(H) + len(M) + (l - len(M)))
            assert lhs == rhs, ("MAX", n, E, H, M, l)
        # HE on H' = H+M
        pk = used
        for b in range(n):
            for a, c in itertools.combinations(range(n), 2):
                if b in (a, c): continue
                ab, bc, ac = fs((a, b)), fs((b, c)), fs((a, c))
                if ab in E and bc in E and ab not in pk and bc not in pk:
                    assert ac in E or ac in pk, "HE"
        cnt += 1
    print("MAX and HE ok on", cnt)
if __name__ == "__main__":
    sm_test(int(sys.argv[1]), 11, int(sys.argv[2])); max_test(int(sys.argv[3]), 12)
