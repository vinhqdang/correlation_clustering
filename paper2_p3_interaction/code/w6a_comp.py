"""Structure of hard components in optimal partitions + tw sandwich (Gamma' = W u hard graph; Q = Gamma'/hard components).
Checks on all optimal partitions of random pair-maximal instances (n <= 10)."""
import sys, random, itertools
from w6a_lib import *
def treewidth(n, adjl):
    # exact tw by DP over subsets: TW(S)=min_v max(TW(S-v), |Q(S-v,v)|)
    N = 1 << n
    def Q(S, v):
        seen = 1 << v; st = [v]; cnt = 0; nb = 0
        stack = [v]; vis = {v}
        while stack:
            x = stack.pop()
            for u in range(n):
                if adjl[x] >> u & 1 and u not in vis:
                    if S >> u & 1: vis.add(u); stack.append(u)
                    else: vis.add(u); nb += 1
        return nb
    tw = {0: -1}
    for S in range(1, N):
        best = 99
        x = S
        while x:
            v = (x & -x).bit_length() - 1; x &= x - 1
            T = S ^ (1 << v)
            val = max(tw[T], Q(T, v))
            if val < best: best = val
        tw[S] = best
    return tw[N - 1]
def main():
    seed, k, a, b = map(int, sys.argv[1:5]); rng = random.Random(seed)
    st = dict(inst=0, parts=0, nonexc_fail=0, split_fail=0, split_cases=0, split_bound_fail=0, uncl_nonexc=0, sand_fail=0, sand_checked=0, clean_bound_fail=0, hc_clean=0, hc_unclean=0)
    for it in range(k):
        n = rng.randint(a, b); adj = random_instance(n, rng); H = maximal_pair_packing(adj, rng)
        if not H: continue
        soft = soft_pairs(H); comp, nc = hard_components(adj, H); hd = hard_deg(adj, H)
        D = max(sum(1 for t in H if v in t) for v in range(n)) * 2
        Dv = [2 * sum(1 for t in H if v in t) for v in range(n)]
        parts, o = all_opt_partitions(adj); st["inst"] += 1
        # cleanliness of components (pair-maximal => hard components may be unclean only through soft edges)
        for c in range(nc):
            vs = [v for v in range(n) if comp[v] == c]
            unclean = any(not (adj[x] >> y & 1) and key(x, y) not in soft for x, y in itertools.combinations(vs, 2))
            st["hc_unclean" if unclean else "hc_clean"] += 1
        for part in parts[:6]:
            st["parts"] += 1
            lab = {}
            for i, c in enumerate(part):
                for v in range(n):
                    if c >> v & 1: lab[v] = i
            Sset = {(x, y) for x in range(n) for y in range(x + 1, n) if (adj[x] >> y & 1) != (lab[x] == lab[y])}
            U = {p for p in Sset if p not in soft}
            ends = {v for p in U for v in p}
            exc = {lab[v] for v in ends}
            for i, c in enumerate(part):
                vs = [v for v in range(n) if c >> v & 1]
                if i in exc: continue
                comps_in = {comp[v] for v in vs}
                # union of whole comps
                for cc in comps_in:
                    if any(comp[v] == cc and not (c >> v & 1) for v in range(n)): st["nonexc_fail"] += 1
                    # clean
                    cv = [v for v in range(n) if comp[v] == cc]
                    if any(not (adj[x] >> y & 1) and key(x, y) not in soft for x, y in itertools.combinations(cv, 2)): st["nonexc_fail"] += 1; st["uncl_nonexc"] += 1
                # cross pairs soft
                for x, y in itertools.combinations(vs, 2):
                    if comp[x] != comp[y] and key(x, y) not in soft: st["nonexc_fail"] += 1
                if len(comps_in) == 1:
                    for y in vs:
                        if len(vs) > hd[y] + Dv[y] + 1: st["clean_bound_fail"] += 1
            for cc in range(nc):
                cv = [v for v in range(n) if comp[v] == cc]
                cls = {lab[v] for v in cv}
                if len(cls) >= 2:
                    st["split_cases"] += 1
                    if not cls <= exc: st["split_fail"] += 1
                    nUK = sum(1 for p in U if comp[p[0]] == cc and comp[p[1]] == cc and adj[p[0]] >> p[1] & 1)
                    mx = 1 + 2 * D + len(U)
                    if len(cv) > (nUK + 1) * mx: st["split_bound_fail"] += 1
        # sandwich
        if n <= 9:
            G2 = [0] * n
            for x in range(n):
                for y in range(n):
                    if x != y and (adj[x] >> y & 1 or key(x, y) in soft): G2[x] |= 1 << y
            Qn = nc; Q2 = [0] * nc
            for x in range(n):
                for y in range(n):
                    if G2[x] >> y & 1 and comp[x] != comp[y]: Q2[comp[x]] |= 1 << comp[y]
            tg = treewidth(n, G2); tq = treewidth(nc, Q2); km = max(sum(1 for v in range(n) if comp[v] == c) for c in range(nc))
            st["sand_checked"] += 1
            if not (tq <= tg <= km * (tq + 1) - 1): st["sand_fail"] += 1
    print(st)
main()
