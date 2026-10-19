"""Adversarial searches.  usage: python3 a4_adv.py part seed steps
part = star   : MAX treewidth tightness family (explicit)
part = hub    : hill-climb to push hard-independent set / alpha(G[N_h(u)]) toward the HUB bounds (n = 9, exhaustive OPT)
part = pieces : hill-climb to push multi-piece clusters toward SM(iii) and my sharper bound (n = 9, all solutions cost <= |H|+3)"""
import sys, random, itertools, collections
import numpy as np
import a4_lib as L

part = sys.argv[1]; seed = int(sys.argv[2]); steps = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
rng = random.Random(seed)

def star_family():
    # centre v in H-member P=(p,v,q); arms a_i - v - b_i; K_{1,2m+2}
    for m in range(1, 6):
        n = 3 + 2 * m; A = np.zeros((n, n), np.int8)
        for x in range(1, n): A[0, x] = A[x, 0] = 1
        H = [(1, 0, 2)]
        M = [(3 + 2 * i, 0, 4 + 2 * i) for i in range(m)]
        assert L.is_packing(A, H + M) and L.is_pair_maximal(A, H + M)
        # M is maximal for H?
        assert all(L.p3_pairs(t) & (L.used_pairs(H) | L.used_pairs(M)) for t in L.p3_list(A))
        w = L.tw_exact(L.ih_adj(H)); w2 = L.tw_exact(L.ih_adj(H + M))
        # soft graph tw
        ws = L.tw_exact(L.soft_adj(n, H + M)) if n <= 15 else None
        print('star m=%d: tw(I(H))=%d tw(I(H u M))=%d  (w+|M|=%d)  tw(W(H u M))=%s  bound 3(w+1)-1 for H u M=%d' % (m, w, w2, w + m, ws, 3 * (w2 + 1) - 1))

def subdivided_k4():
    # M = 4 P3s whose intersection pattern is K4, H = 6 P3 on the K4 edges (each H-P3 meets two M-P3s): I(H) edgeless.
    n = 0; vid = {}
    # K4 vertices 0..3 = M-P3s; edges e=(i,j): H-P3 E_e with two private vertices u_e^i (shared with M_i) and u_e^j, plus a middle one.
    # M_i has three vertices, one per incident edge.  H-P3 on (x, y, z): x in M_i, z in M_j, y private, need induced path x-y-z.
    edges = list(itertools.combinations(range(4), 2)); cnt = 0; Hv = []; Mv = {i: [] for i in range(4)}
    for (i, j) in edges:
        x, y, z = cnt, cnt + 1, cnt + 2; cnt += 3; Hv.append((x, y, z)); Mv[i].append(x); Mv[j].append(z)
    n = cnt; A = np.zeros((n, n), np.int8)
    for (x, y, z) in Hv: A[x, y] = A[y, x] = A[y, z] = A[z, y] = 1
    # M_i: path through its three vertices (a - b - c): edges ab, bc, a-c non-edge (pairs are between different H-P3 so hard)
    Ms = []
    for i in range(4):
        a, b, c = Mv[i]; A[a, b] = A[b, a] = A[b, c] = A[c, b] = 1; Ms.append((a, b, c))
    assert L.is_packing(A, Hv + Ms), 'not a packing'
    w = L.tw_exact(L.ih_adj(Hv)); w2 = L.tw_exact(L.ih_adj(Hv + Ms))
    print('subdivided K4: n=%d |H|=6 |M|=4 tw(I(H))=%d tw(I(H u M))=%d' % (n, w, w2))

def build_H_for_hub(A, rng):
    """pack P3s whose end pairs lie inside N(0) first (to create soft pairs among hub leaves), then extend"""
    n = len(A); N0 = {x for x in range(n) if A[0][x]}
    ts = [t for t in L.p3_list(A) if t[0] in N0 and t[2] in N0 and t[1] != 0]
    rng.shuffle(ts); H = []; used = set()
    for t in ts:
        ps = L.p3_pairs(t)
        if not (ps & used): H.append(t); used |= ps
    return L.greedy_extend(A, H, rng)

def hub_climb():
    n = 9; lab, same, P = L.part_table(n)[0], L.part_table(n)[1], L.part_table(n)[2]
    best = {}
    for l in (0, 1, 2):
        A = L.gen_graph(rng, n, 'hub'); cur = -1; bestA = None
        for s in range(steps):
            B = A.copy()
            for _ in range(rng.randint(1, 2)):
                i, j = rng.sample(range(n), 2); B[i, j] = B[j, i] = 1 - B[i, j]
            H = build_H_for_hub(B, rng); pk = L.packed_matrix(n, H)
            _, _, cost = L.all_costs(B); OPT = int(cost.min())
            if OPT > len(H) + l: continue
            c = L.cmax(n, H); D = L.soft_degree(pk)
            hs = L.hub_stats(B, pk, H)
            X = max(x[0] for x in hs); al = max(x[1] for x in hs)
            assert X <= 2 * l + 1 and al <= (c + 1) * (2 * l + 1) and al <= (D + 1) * (2 * l + 1)
            score = al / ((c + 1) * (2 * l + 1)) + 0.01 * X / (2 * l + 1) + 0.001 * al
            if score >= cur: cur = score; A = B; bestA = (B, H, X, al, c, D, OPT)
        B, H, X, al, c, D, OPT = bestA
        print('HUB climb l=%d: best |X|=%d (bound %d), alpha=%d, c=%d, D=%d, (c+1)(2l+1)=%d, (D+1)(2l+1)=%d, |H|=%d OPT=%d'
              % (l, X, 2 * l + 1, al, c, D, (c + 1) * (2 * l + 1), (D + 1) * (2 * l + 1), len(H), OPT))

def pieces_climb():
    n = 9; labs = L.part_table(n)[0]
    A = L.gen_graph(rng, n, 'fano'); H = None; best = -10 ** 9; bestinfo = None
    for s in range(steps):
        B = A.copy()
        for _ in range(rng.randint(0, 2)):
            i, j = rng.sample(range(n), 2); B[i, j] = B[j, i] = 1 - B[i, j]
        H = L.rand_packing(B, rng, stop=rng.choice([0, .3]), prefer_overlap=rng.random() < .5)
        pk = L.packed_matrix(n, H); h = len(H); D = L.soft_degree(pk)
        _, _, cost = L.all_costs(B)
        idxs = np.nonzero(cost <= h + 3)[0]
        if len(idxs) > 300: idxs = rng.sample(list(idxs), 300)
        score = -10 ** 9; info = None
        for r in idxs:
            v, st = L.sm_check(B, pk, labs[r], D)
            assert not v, v
            if st['multi']:
                u = st['u']; sharp = max(D + u + 1, 2 * D, 2 * u); orig = max(2 * D + u, 2 * u)
                sc = st['maxC'] - sharp
                if sc > score: score = sc; info = (st['maxC'], D, u, sharp, orig)
        if score >= best: best = score; A = B; bestinfo = info
    print('PIECES climb: best |C|-sharp = %d, info (|C|, Ds, |U|, sharp bound, orig bound) = %s' % (best, bestinfo))

{'star': lambda: (star_family(), subdivided_k4()), 'hub': hub_climb, 'pieces': pieces_climb}[part]()
