"""Sanity check of the CSP view of Cluster Editing above a modification-disjoint
P3 packing at l = 0 (brute force on small random graphs).

Claims checked
  (A) every triple of vertices has its packed pairs in packed P3s that pairwise
      share a vertex (so a constraint only links P3s that intersect);
  (B) 'some set S with exactly one edit per packed P3 and no other edit makes G a
      cluster graph' agrees between (i) enumerating the 3^|H| local choices and
      (ii) the optimum over all partitions being equal to |H|.
"""
import itertools, random, sys

def pairs(P):  # P = (a, b, c) induced path a-b-c: edges ab, bc, non-edge ac
    a, b, c = P
    return [frozenset((a, b)), frozenset((b, c)), frozenset((a, c))]

def induced_p3s(n, E):
    out = []
    for b in range(n):
        for a, c in itertools.combinations(range(n), 2):
            if a != b and c != b and frozenset((a, b)) in E and frozenset((b, c)) in E and frozenset((a, c)) not in E:
                out.append((a, b, c))
    return out

def packing(n, E, rng):
    cand = induced_p3s(n, E); rng.shuffle(cand)
    used, H = set(), []
    for P in cand:
        ps = pairs(P)
        if not any(p in used for p in ps):
            H.append(P); used.update(ps)
    return H

def is_cluster(n, E):
    for a, b, c in itertools.permutations(range(n), 3):
        if frozenset((a, b)) in E and frozenset((b, c)) in E and frozenset((a, c)) not in E:
            return False
    return True

def partitions(n):
    def rec(i, labels, m):
        if i == n:
            yield list(labels); return
        for l in range(m + 1):
            labels.append(l); yield from rec(i + 1, labels, max(m, l + 1)); labels.pop()
    yield from rec(0, [], 0)

def opt(n, E):
    best = 10**9
    for lab in partitions(n):
        cost = sum((lab[a] == lab[b]) != (frozenset((a, b)) in E) for a, b in itertools.combinations(range(n), 2))
        best = min(best, cost)
    return best

def main(trials=400, seed=1):
    rng = random.Random(seed); checked = hit = 0
    for _ in range(trials):
        n = rng.randint(4, 7); p = rng.choice([0.25, 0.4, 0.55])
        E = {frozenset(e) for e in itertools.combinations(range(n), 2) if rng.random() < p}
        H = packing(n, E, rng)
        if not H: continue
        owner = {q: i for i, P in enumerate(H) for q in pairs(P)}
        for tri in itertools.combinations(range(n), 3):
            ps = [frozenset(q) for q in itertools.combinations(tri, 2)]
            owners = sorted({owner[q] for q in ps if q in owner})
            for i, j in itertools.combinations(owners, 2):
                assert set(H[i]) & set(H[j]), ("(A) fails", H, tri)
        sat = False
        for ch in itertools.product(range(3), repeat=len(H)):
            F = set(E)
            for P, c in zip(H, ch):
                F ^= {pairs(P)[c]}
            if is_cluster(n, F): sat = True; break
        assert sat == (opt(n, E) == len(H)), ("(B) fails", n, sorted(map(sorted, E)), H)
        checked += 1; hit += sat
    print(f"checked {checked} instances, {hit} with OPT = |H|; (A) and (B) hold")

if __name__ == "__main__":
    main(*(int(x) for x in sys.argv[1:]))
