"""Audit tests: Obs 1, Obs 2, Lemma 0, 1, 2, 3, 4, F1, critical cliques, C3, C4, C1."""
import itertools, random, sys
import numpy as np
from audit_lib import *

def owner_map(n, H):
    P = pidx(n); own = {}
    for i, p in enumerate(H):
        a, b, c = p
        for (x, y) in ((a, b), (b, c), (a, c)):
            own[P[(x, y)]] = i
    return own

def run(seed, trials, nlo=5, nhi=9):
    rng = random.Random(seed)
    stats = dict(inst=0, obs1_equiv_yes=0, obs1_equiv_no=0, obs2_triples=0, lem2_sols=0, lem2_pairs=0,
                 lem4_pairs=0, lem4_inst=0, f1=0, crit=0, packed_multi2=0, packed_multi3=0,
                 forced_multi2_inst=0, forced_multi3_inst=0, lem0_sols=0)
    for t in range(trials):
        n = rng.randint(nlo, nhi)
        kind = rng.choice(['gnp', 'hub', 'cliques', 'path'])
        E = rand_graph(rng, n, kind)
        mode = rng.choice(['max', 'part', 'overlap', 'max'])
        H = rand_packing(rng, n, E, mode)
        if not H: continue
        pk = check_packing(n, E, H)
        Pm = pidx(n); p = len(H)
        stats['inst'] += 1
        pms = [pairmask(n, q) for q in H]
        own = owner_map(n, H)
        # ---- Obs 2
        Vsets = [set(q) for q in H]
        for tri in itertools.combinations(range(n), 3):
            u, v, w = tri
            ow = {own[Pm[x]] for x in ((u, v), (v, w), (u, w)) if Pm[x] in own}
            for i in ow:
                for j in ow:
                    if i < j: assert Vsets[i] & Vsets[j], ('Obs2 violated', H, tri)
            stats['obs2_triples'] += 1
        # ---- exhaustive partitions
        S, cost = opt_and_costs(n, E)
        opt = int(cost.min())
        assert opt >= p
        # Lemma 0 / Obs 1 on all solutions with cost <= p+3
        sel = np.nonzero(cost <= p + 3)[0]
        unpk = np.uint64(pk)
        notpk = ~unpk
        Ssel = S[sel]; csel = cost[sel]
        tot = np.zeros(len(sel), dtype=np.int64)
        for pm in pms:
            e = popcount(Ssel & np.uint64(pm)).astype(np.int64)
            assert (e >= 1).all(), 'Lemma0: some P3 unedited in a cluster-editing set'
            tot += e
        U = popcount(Ssel & notpk).astype(np.int64)
        assert (tot + U == csel).all()
        stats['lem0_sols'] += len(sel)
        # Obs 1: opt == p iff exists x in {1,2,3}^p cluster
        if p <= 9:
            exists = False
            opts = [[1 << 0, 1 << 1, 1 << 2] for _ in H]
            # explicit pair bits per P3
            bits = []
            for q in H:
                a, b, c = q
                bits.append([1 << Pm[(a, b)], 1 << Pm[(b, c)], 1 << Pm[(a, c)]])
            for ch in itertools.product(range(3), repeat=p):
                Sx = 0
                for i, c in enumerate(ch): Sx |= bits[i][c]
                if is_cluster(n, E ^ Sx): exists = True; break
            assert exists == (opt == p), ('Obs1', exists, opt, p)
            stats['obs1_equiv_yes' if exists else 'obs1_equiv_no'] += 1
            if opt == p:
                s0 = S[cost == p]
                for pm in pms: assert (popcount(s0 & np.uint64(pm)) == 1).all()
                assert (popcount(s0 & notpk) == 0).all()
        # multi-edit: do solutions of cost <= opt have packed P3 with >=2 edits? is it forced?
        optsel = np.nonzero(cost == opt)[0]
        m2 = np.zeros(len(optsel), dtype=bool); m3 = np.zeros(len(optsel), dtype=bool)
        for pm in pms:
            e = popcount(S[optsel] & np.uint64(pm))
            m2 |= e >= 2; m3 |= e >= 3
        if m2.any(): stats['packed_multi2'] += 1
        if m3.any(): stats['packed_multi3'] += 1
        if m2.all(): stats['forced_multi2_inst'] += 1
        if m3.all(): stats['forced_multi3_inst'] += 1
        # ---- Q(u), Lemma 2, 3
        Q = [set() for _ in range(n)]
        for q in H:
            for x in q:
                for y in q:
                    if x != y: Q[x].add(y)
        pkset = {k for k in range(n * (n - 1) // 2) if (pk >> k) & 1}
        sols = list(sel)
        rng.shuffle(sols)
        for si in sols[:200]:
            Sx = int(S[si])
            Uset = [k for k in range(n * (n - 1) // 2) if (Sx >> k) & 1 and k not in pkset]
            kk = {Pm[(a, b)]: (a, b) for a in range(n) for b in range(a + 1, n)}
            for k in Uset:
                u, v = kk[k]
                Ep = E ^ (1 << k)
                D = []
                for z in range(n):
                    if z in (u, v) or z in Q[u] or z in Q[v]: continue
                    t3 = [x for x in induced_p3s(n, Ep) if set(x) == {u, v, z}]
                    if t3: D.append(z)
                for z in D:
                    assert ((Sx >> Pm[(u, z)]) & 1) or ((Sx >> Pm[(v, z)]) & 1), 'Lemma2'
                assert len(D) <= len(Uset) - 1, 'Lemma2 count'
                # Lemma 3
                for z in range(n):
                    if z in (u, v): continue
                    if Pm[(u, z)] in pkset or Pm[(v, z)] in pkset:
                        assert z in Q[u] or z in Q[v]
                stats['lem2_pairs'] += 1
            stats['lem2_sols'] += 1
        # ---- Lemma 4 for maximal H
        maximal = all(any(pairmask(n, q3) & pk for _ in [0]) for q3 in induced_p3s(n, E))
        if maximal:
            Vh = set().union(*Vsets)
            R = [x for x in range(n) if x not in Vh]
            N = [{y for y in range(n) if y != x and has(n, E, x, y)} for x in range(n)]
            if len(R) >= 2: stats['lem4_inst'] += 1
            for r, s in itertools.combinations(R, 2):
                if has(n, E, r, s): assert N[r] | {r} == N[s] | {s}, 'Lemma4 edge'
                else: assert not (N[r] & N[s]), 'Lemma4 nonedge'
                stats['lem4_pairs'] += 1
            # (ii)
            comp = {}
            for x in range(n):
                if x in comp: continue
                st = [x]; comp[x] = x
                while st:
                    y = st.pop()
                    for z in N[y]:
                        if z not in comp: comp[z] = x; st.append(z)
            for r in R:
                cc = {y for y in range(n) if comp[y] == comp[r]}
                if cc & Vh: assert N[r] & Vh, 'Lemma4(ii)'
        # ---- F1: every optimum has clusters connected in G ; critical clique property
        labs = all_partitions(n)[1]
        optidx = np.nonzero(cost == opt)[0]
        Nc = [frozenset(N_ | {x}) for x, N_ in enumerate([{y for y in range(n) if y != x and has(n, E, x, y)} for x in range(n)])]
        twin_ok = False
        for oi in optidx:
            lab = labs[oi]
            for c in set(lab):
                vs = [x for x in range(n) if lab[x] == c]
                # connected in G?
                seen = {vs[0]}; st = [vs[0]]
                while st:
                    y = st.pop()
                    for z in vs:
                        if z not in seen and has(n, E, y, z): seen.add(z); st.append(z)
                assert len(seen) == len(vs), 'F1 violated'
            if all(lab[x] == lab[y] for x in range(n) for y in range(n) if Nc[x] == Nc[y]): twin_ok = True
        assert twin_ok, 'critical clique property'
        stats['f1'] += 1; stats['crit'] += 1
    return stats

if __name__ == '__main__':
    seed = int(sys.argv[1]); trials = int(sys.argv[2])
    print(run(seed, trials))
