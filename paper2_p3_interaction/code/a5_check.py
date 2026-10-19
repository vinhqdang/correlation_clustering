"""Independent brute-force audit of wave5_locality claims (n<=8 random, partitions enumerated)."""
import sys, random, itertools
from a5_lib import *
def run(n, N, seed0=0):
    st = dict(inst=0, optparts=0, L1a=0, L1b=0, T1T2=0, cor2=0, cor2pairs=0, L3=0, maxpart=0, maxpart_bound=0,
              partner_gt_DD2=0, ctrl_L1=0, ctrl_parts=0, nonlocal_pairs_in_tiny_with_soft=0)
    for sd in range(seed0, seed0 + N):
        rng = random.Random(f"a5-{n}-{sd}")
        p = rng.uniform(.25, .7)
        adj = mk(n, [(i, j) for i in range(n) for j in range(i+1, n) if rng.random() < p])
        H = rand_pm(adj, rng, rng.choice([0, .3, .6]))
        if not H or not pair_maximal(adj, H): continue
        pk = packed(H); st['inst'] += 1
        Q = {v: {x for x in range(n) if frozenset((v, x)) in pk} for v in range(n)}
        D = max(len(q) for q in Q.values()); cnt = [0]*n
        for t in H:
            for v in t: cnt[v] += 1
        c = max(cnt)
        opt, opts = optima(adj)
        for lab in partitions(n) if False else opts:
            st['optparts'] += 1
            cl = {}
            for v in range(n): cl.setdefault(lab[v], []).append(v)
            # Lemma 1
            for C in cl.values():
                for v in C:
                    d = sum(adj[v][x] for x in C if x != v)
                    if 2*d < len(C) - 1: st['L1a'] += 1
                for a, b in itertools.combinations(C, 2):
                    if not adj[a][b] and not any(adj[a][z] and adj[z][b] for z in C if z not in (a, b)): st['L1b'] += 1
            pid = pieces(adj, pk, lab); U = U_set(adj, pk, lab); u = len(U)
            # Thm 1
            for (a, b) in U:
                if adj[a][b]: continue
                C = cl[lab[a]]; ok = False
                for z in C:
                    if z in (a, b) or not (adj[a][z] and adj[z][b]): continue
                    sa, sb = frozenset((a, z)) in pk, frozenset((z, b)) in pk
                    if sa and sb: ok = True
                    elif sa != sb:
                        h, o = (a, b) if not sa else (b, a)   # h is the endpoint with hard edge to z
                        if pid[h] == pid[z]: ok = True
                if not ok: st['T1T2'] += 1
            # Cor 2
            for C in cl.values():
                ps = {}
                for v in C: ps.setdefault(pid[v], []).append(v)
                for A, B in itertools.combinations(list(ps.values()), 2):
                    soft = any(frozenset((a, b)) in pk for a in A for b in B)
                    if soft:
                        # control: how often does a cross inserted pair have no common soft neighbour?
                        for a in A:
                            for b in B:
                                if frozenset((a, b)) not in pk and not (Q[a] & Q[b]): st['nonlocal_pairs_in_tiny_with_soft'] += 1
                        continue
                    st['cor2'] += 1
                    for a in A:
                        for b in B:
                            st['cor2pairs'] += 1
                            if (a, b) not in U and (b, a) not in U: st['cor2'] += 10**6
                            zs = [z for z in Q[a] & Q[b] if z in C and z not in A and z not in B]
                            if not zs: st['cor2'] += 10**6
                            # partner count
                    # partner counts
            # partner bound: for any a, number of b with common soft nbr and ab hard
            for a in range(n):
                partners = set()
                for z in Q[a]:
                    for b in Q[z]:
                        if b != a and frozenset((a, b)) not in pk: partners.add(b)
                st['maxpart'] = max(st['maxpart'], len(partners))
                if len(partners) > D*(D-2) and D >= 2: st['partner_gt_DD2'] += 1
            # Lemma 3
            for C in cl.values():
                for z, b in itertools.permutations(C, 2):
                    M = [a for a in range(n) if a not in (z, b) and frozenset((a, z)) not in pk and frozenset((a, b)) not in pk
                         and adj[a][z] and adj[a][b] == 0 or (a not in (z, b) and frozenset((a, z)) not in pk and frozenset((a, b)) not in pk and adj[a][z] == 0 and adj[a][b])]
                    if len(M) > u: st['L3'] += 1
        # control: suboptimal partitions (cost opt+1) - Lemma 1 should fail somewhere
        if inst_ctrl(adj, opt): st['ctrl_L1'] += 1
    return st
def inst_ctrl(adj, opt):
    n = len(adj)
    for lab in partitions(n):
        if cost(adj, lab) != opt + 1: continue
        cl = {}
        for v in range(n): cl.setdefault(lab[v], []).append(v)
        for C in cl.values():
            for v in C:
                if 2*sum(adj[v][x] for x in C if x != v) < len(C) - 1: return True
    return False
if __name__ == '__main__':
    n = int(sys.argv[1]); N = int(sys.argv[2])
    print(n, run(n, N))
