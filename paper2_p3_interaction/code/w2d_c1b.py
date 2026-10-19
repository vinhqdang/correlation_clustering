"""Targeted C1 search with large R (core + attached cliques / overlapping cliques / complete multipartite minus noise).
Tests: for every random maximal H (3 per graph) with |R|>=2, l>0: min #unpacked R-R edits over optima. usage: NMIN NMAX PER"""
import sys, json, random, itertools, collections
from multiprocessing import Pool
from w2d_common import *

def gen(fam, rng, n):
    if fam == "core_cliques":
        k = rng.randint(4, 6); E = set(gnp(rng, k, rng.choice([.4, .6]))); nxt = k
        while nxt < n:
            s = min(n - nxt, rng.randint(2, 4)); B = list(range(nxt, nxt + s)); nxt += s
            E |= set(itertools.combinations(B, 2))
            att = rng.sample(range(k), rng.randint(1, 3))
            for b in B:
                for a in att:
                    if rng.random() < .85: E.add((a, b))
        return frozenset(E)
    if fam == "overlap_cliques":
        k = rng.randint(2, 3); sets = [set(rng.sample(range(n), rng.randint(3, n // 2 + 1))) for _ in range(k)]
        E = set()
        for S in sets: E |= set(itertools.combinations(sorted(S), 2))
        for q in itertools.combinations(range(n), 2):
            if rng.random() < .05: E ^= {q}
        return frozenset(E)
    if fam == "multipartite":
        k = rng.randint(2, 4); part = [rng.randrange(k) for _ in range(n)]
        E = {q for q in itertools.combinations(range(n), 2) if part[q[0]] != part[q[1]]}
        for q in itertools.combinations(range(n), 2):
            if rng.random() < .06: E ^= {q}
        return frozenset(E)
    if fam == "star_cliques":
        hubs = rng.randint(1, 2); E = set(); nxt = hubs
        while nxt < n:
            s = min(n - nxt, rng.randint(2, 4)); B = list(range(nxt, nxt + s)); nxt += s
            E |= set(itertools.combinations(B, 2))
            for h in rng.sample(range(hubs), rng.randint(1, hubs)):
                for b in B:
                    if rng.random() < .8: E.add((h, b))
        return frozenset(E)

def run(args):
    fam, n, seed = args
    rng = random.Random(f"c1b-{fam}-{n}-{seed}")
    E = gen(fam, rng, n); O = opt(n, E); ps = list(itertools.combinations(range(n), 2)); idx = {p: i for i, p in enumerate(ps)}
    out = []
    for t in range(3):
        H = rand_maximal_packing(rng, n, E); l = O - len(H); V = vh(H); R = [v for v in range(n) if v not in V]
        if l <= 0 or len(R) < 2: continue
        PP = packed_pairs(H)
        RR = np.array([1.0 if (p not in PP and p[0] not in V and p[1] not in V) else 0 for p in ps])
        S = solve(n, E, w=RR, cap=O)[1]; bad = int(sum(RR[idx[p]] for p in S))
        out.append(dict(fam=fam, n=n, seed=seed, t=t, l=l, R=len(R), minRR=bad))
        if bad: json.dump(dict(n=n, E=sorted(E), H=H, opt=O, l=l, note="C1 counterexample"), open(f"instances_w2d/c1b_counter_{fam}_{n}_{seed}_{t}.json", "w"))
    return out

if __name__ == "__main__":
    nmin, nmax, per = map(int, sys.argv[1:4])
    jobs = [(f, n, s) for n in range(nmin, nmax + 1) for f in ("core_cliques", "overlap_cliques", "multipartite", "star_cliques") for s in range(per)]
    with Pool(4) as p: res = [r for rs in p.imap_unordered(run, jobs, chunksize=4) for r in rs]
    json.dump(res, open(f"instances_w2d/c1b_results_{nmin}_{nmax}_{per}.json", "w"))
    print("tested", len(res), "counterexamples", sum(r["minRR"] > 0 for r in res), "R>=3:", sum(r["R"] >= 3 for r in res), "R>=5:", sum(r["R"] >= 5 for r in res))
    print(collections.Counter((r["fam"], r["n"]) for r in res))
