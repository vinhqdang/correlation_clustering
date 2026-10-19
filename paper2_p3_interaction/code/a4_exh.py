"""Exhaustive audit (all set partitions) of MAX, HE, HUB, SM, soft-graph treewidth.
usage: python3 a4_exh.py seed ninst nmin nmax"""
import sys, random, collections, time, itertools
import numpy as np
import a4_lib as L

seed, ninst, nmin, nmax = map(int, sys.argv[1:5])
rng = random.Random(seed)
C = collections.Counter(); worst = collections.defaultdict(float); bad = []
t0 = time.time()
for it in range(ninst):
    n = rng.randint(nmin, nmax)
    A = L.gen_graph(rng, n)
    if rng.random() < .5: A = L.noise(A, rng, rng.randint(1, 3))
    mode = rng.choice(['maximal', 'partial', 'overlap', 'overlap_partial'])
    H = L.rand_packing(A, rng, stop=.5 if 'partial' in mode else 0.0, prefer_overlap='overlap' in mode)
    if mode == 'maximal': H = L.greedy_extend(A, H, rng)
    if mode == 'overlap': H = L.greedy_extend(A, H, rng, True)
    assert L.is_packing(A, H)
    pk = L.packed_matrix(n, H); h = len(H)
    lab, same, cost = L.all_costs(A); OPT = int(cost.min())
    C['inst'] += 1; C['maximal'] += L.is_pair_maximal(A, H)
    D = L.soft_degree(pk); c = L.cmax(n, H)
    assert D == 2 * c
    # ---------- MAX
    for rep in range(3):
        M = L.rand_M(A, H, rng); HM = H + M
        assert L.is_packing(A, HM) and L.is_pair_maximal(A, HM)           # H u M packing and pair-maximal
        assert OPT >= h + len(M)                                          # Lemma 0 for the enlarged packing
        for l in range(0, 4):
            lhs = OPT <= h + l
            rhs = len(M) <= l and OPT <= len(HM) + (l - len(M))
            if lhs != rhs: bad.append(('MAX-equiv', seed, it, l))
        if len(H) <= 11 and len(HM) <= 13:
            w = L.tw_exact(L.ih_adj(H)); w2 = L.tw_exact(L.ih_adj(HM))
            if w2 > w + len(M): bad.append(('MAX-tw', seed, it, w, w2, len(M)))
            C['maxtw_checked'] += 1; worst['tw_increase_minus_M'] = max(worst['tw_increase_minus_M'], w2 - w - len(M))
            C['tw_incr_eq_M'] += (w2 - w == len(M) and len(M) > 0)
        C['M_total'] += len(M); C['M_nonempty'] += bool(M)
        # HE on the pair-maximal packing HM
        pkM = L.packed_matrix(n, HM)
        for a, b, cc in itertools.combinations(range(n), 3):
            for x, y, z in ((a, b, cc), (b, a, cc), (a, cc, b)):    # centre y
                if A[x][y] and A[y][z] and not pkM[x, y] and not pkM[y, z]:
                    if not (A[x][z] or pkM[x, z]): bad.append(('HE', seed, it))
                    C['HE_triples'] += 1
    # HE on H itself when H is pair-maximal; and failure count when not
    for a, b, cc in itertools.combinations(range(n), 3):
        for x, y, z in ((a, b, cc), (b, a, cc), (a, cc, b)):
            if A[x][y] and A[y][z] and not pk[x, y] and not pk[y, z] and not (A[x][z] or pk[x, z]):
                if L.is_pair_maximal(A, H): bad.append(('HE-H', seed, it))
                else: C['HE_fails_nonmaximal'] += 1
    # ---------- HUB (needs YES instance with budget l)
    hs = L.hub_stats(A, pk, H)
    for l in range(0, 4):
        if OPT <= h + l:
            C['yes_l%d' % l] += 1
            for u, (X, al, nh) in enumerate(hs):
                if X > 2 * l + 1: bad.append(('HUB-X', seed, it, l, u, X))
                if al > (D + 1) * (2 * l + 1): bad.append(('HUB-alpha', seed, it, l, u))
                if al > (c + 1) * (2 * l + 1): bad.append(('HUB-alpha-sharp', seed, it, l, u, al, c))
                worst['X/(2l+1) l=%d' % l] = max(worst['X/(2l+1) l=%d' % l], X / (2 * l + 1))
                worst['alpha/((c+1)(2l+1)) l=%d' % l] = max(worst['alpha/((c+1)(2l+1)) l=%d' % l], al / ((c + 1) * (2 * l + 1)))
    # ---------- soft-graph treewidth
    if h <= 13 and h > 0:
        w = L.tw_exact(L.ih_adj(H)); tw = L.tw_exact(L.soft_adj(n, H))
        C['softtw_checked'] += 1
        if tw > max(0, 3 * (w + 1) - 1): bad.append(('softtw', seed, it, w, tw))
        C['softtw_tight'] += (tw == 3 * (w + 1) - 1)
        worst['tw(W)/(3(w+1)-1)'] = max(worst['tw(W)/(3(w+1)-1)'], tw / (3 * (w + 1) - 1))
    # ---------- SM on all solutions with cost <= h+3
    idxs = np.nonzero(cost <= h + 3)[0]
    if len(idxs) > 250: idxs = rng.sample(list(idxs), 250)
    P = L.part_table(n)[2]
    for r in idxs:
        v, st = L.sm_check(A, pk, lab[r], D)
        C['sol'] += 1; C['multi_clusters'] += st['multi']
        if st['multi']:
            worst['SM maxC multi'] = max(worst['SM maxC multi'], st['maxC'])
        for x in v: bad.append(('SM', seed, it, x))
print('seed', seed, 'time %.0fs' % (time.time() - t0))
print(dict(C)); print({k: round(v, 3) for k, v in worst.items()})
print('VIOLATIONS', len(bad), bad[:10])
