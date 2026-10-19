"""ILP (HiGHS) audit, n = 10..13, same statements.  usage: python3 a4_ilp.py seed ninst nmin nmax"""
import sys, random, collections, time, itertools
import numpy as np
import a4_lib as L

seed, ninst, nmin, nmax = map(int, sys.argv[1:5])
rng = random.Random(seed); C = collections.Counter(); worst = collections.defaultdict(float); bad = []
t0 = time.time()
def hard_components(A, pk):
    n = len(A); comp = list(range(n))
    def f(x):
        while comp[x] != x: comp[x] = comp[comp[x]]; x = comp[x]
        return x
    for i, j in L.prs(n):
        if A[i][j] and not pk[i, j]: comp[f(i)] = f(j)
    return [f(i) for i in range(n)]
for it in range(ninst):
    n = rng.randint(nmin, nmax)
    A = L.gen_graph(rng, n)
    if rng.random() < .6: A = L.noise(A, rng, rng.randint(1, 4))
    mode = rng.choice(['maximal', 'partial', 'overlap', 'overlap_partial'])
    H = L.rand_packing(A, rng, stop=.4 if 'partial' in mode else 0.0, prefer_overlap='overlap' in mode)
    if mode in ('maximal', 'overlap'): H = L.greedy_extend(A, H, rng, mode == 'overlap')
    pk = L.packed_matrix(n, H); h = len(H); D = L.soft_degree(pk); c = L.cmax(n, H)
    r = L.ilp(A); OPT = r[1]; C['inst'] += 1
    # MAX
    M = L.rand_M(A, H, rng); HM = H + M
    assert L.is_packing(A, HM) and L.is_pair_maximal(A, HM)
    if OPT < h + len(M): bad.append(('MAX-lb', seed, it))
    for l in range(4):
        if (OPT <= h + l) != (len(M) <= l and OPT <= len(HM) + l - len(M)): bad.append(('MAX-eq', seed, it))
    if len(H) <= 14 and len(HM) <= 15:
        w = L.tw_exact(L.ih_adj(H)); w2 = L.tw_exact(L.ih_adj(HM)); C['tw'] += 1
        if w2 > w + len(M): bad.append(('MAX-tw', seed, it))
        worst['tw incr - |M|'] = max(worst['tw incr - |M|'], w2 - w - len(M))
        if len(H) <= 14 and len(H) > 0:
            tw = L.tw_exact(L.soft_adj(n, H)); C['softtw'] += 1
            if tw > max(0, 3 * (w + 1) - 1): bad.append(('softtw', seed, it, w, tw))
            worst['tw(W)/(3(w+1)-1)'] = max(worst['tw(W)/(3(w+1)-1)'], tw / (3 * (w + 1) - 1))
    # HE
    pkM = L.packed_matrix(n, HM)
    for a, b, cc in itertools.combinations(range(n), 3):
        for x, y, z in ((a, b, cc), (b, a, cc), (a, cc, b)):
            if A[x][y] and A[y][z] and not pkM[x, y] and not pkM[y, z] and not (A[x][z] or pkM[x, z]): bad.append(('HE', seed, it))
    # HUB
    hs = L.hub_stats(A, pk, H)
    for l in range(4):
        if OPT <= h + l:
            C['yes_l%d' % l] += 1
            for u, (X, al, nh) in enumerate(hs):
                if X > 2 * l + 1 or al > (D + 1) * (2 * l + 1) or al > (c + 1) * (2 * l + 1): bad.append(('HUB', seed, it, l, u, X, al))
                worst['alpha/((c+1)(2l+1)) l=%d' % l] = max(worst['alpha/((c+1)(2l+1)) l=%d' % l], al / ((c + 1) * (2 * l + 1)))
    # SM : optimum, random-objective solutions with cost <= h + l, forced multi-piece maximum cluster
    hc = hard_components(A, pk); P = L.prs(n)
    for l in range(0, 4):
        lim = h + l
        if OPT > lim: continue
        sols = []
        for k in range(2):
            obj = np.array([rng.uniform(-1, 1) for _ in P]); rr = L.ilp(A, obj=obj, cost_le=lim)
            if rr: sols.append(rr)
        # forced multi-piece: a,b in different hard components, same cluster, maximise cluster of a
        cand = [(a, b) for a, b in P if hc[a] != hc[b]]
        for a, b in rng.sample(cand, min(3, len(cand))) if cand else []:
            obj = np.zeros(len(P))
            idx = L._tri(n)[1]
            for v in range(n):
                if v != a: obj[idx[L.key(a, v)]] = -1
            rr = L.ilp(A, obj=obj, cost_le=lim, fixed={(a, b): 1})
            if rr: sols.append(rr); C['forced_multi_l%d' % l] += 1
        for z, cost, idx in sols:
            lab = L.labels_from_z(n, z, idx)
            v, st = L.sm_check(A, pk, lab, D); C['sol'] += 1; C['multi'] += st['multi']
            if st['multi']: worst['maxC multi'] = max(worst['maxC multi'], st['maxC'])
            for x in v: bad.append(('SM', seed, it, l, x))
print('seed', seed, 'time %.0fs' % (time.time() - t0)); print(dict(C)); print({k: round(v, 3) for k, v in worst.items()})
print('VIOLATIONS', len(bad), bad[:10])
