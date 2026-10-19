"""Exhaustive (all set partitions) audit of Twin lemma, C1 strong form, exchange formula,
N3/N4/N5, P2/P3/P4.  usage: python3 a3a_exh.py seed ninst nmin nmax [hyp]
hyp = 'pairmax' (default, H maximal w.r.t. pair-disjointness) | 'vertexmax' | 'nonmax' (random sub-packing)"""
import sys, random, itertools, collections
import numpy as np
import a3a_lib as L

def canon(lab):
    m = {}; out = []
    for x in lab:
        if x not in m: m[x] = len(m)
        out.append(m[x])
    return tuple(out)

def run(seed, ninst, nmin, nmax, hyp='pairmax'):
    rng = random.Random(seed)
    C = collections.Counter(); fails = []
    for it in range(ninst):
        n = rng.randint(nmin, nmax)
        f, A, H = L.make_instance(rng, n)
        if hyp == 'vertexmax':
            H = L.vertex_greedy_packing(A, rng)
        elif hyp == 'nonmax':
            H = H[:rng.randint(0, len(H))]
        VH, R, Q, classes = L.info(A, H)
        pm = L.is_pair_maximal(A, H)
        if hyp == 'pairmax': assert pm
        if len(R) < 2: continue
        l4 = L.lemma4_holds(A, R)
        if pm: assert l4
        C['inst'] += 1; C['pairmax' if pm else 'notpairmax'] += 1
        if not l4: C['lemma4_fail'] += 1
        costs = L.costs_all(A); OPT = int(costs.min()); opt = np.nonzero(costs == OPT)[0]
        Lab = L.rgs(n); S = L.same_matrix(n)
        pi = L.pidx(n); ed_all = np.array([A[i, j] for i, j in itertools.combinations(range(n), 2)])
        lopt = OPT - len(H); assert lopt >= 0
        C['opt_gt_H'] += lopt > 0
        # twin lemma, every optimum, all true-twin pairs of G
        tw = L.true_twin_pairs(A)
        for (u, v) in tw:
            C['twin_pairs'] += 1
            if not S[opt, pi[(u, v)]].all(): C['TWIN_FAIL'] += 1; fails.append(('twin', A.tolist(), H))
        # NF violation: R-R pairs inserted or deleted
        rr = [pi[(r, s)] for r, s in itertools.combinations(R, 2)]
        viol = S[:, rr].astype(np.int32) @ (1 - ed_all[rr]) + (1 - S[:, rr].astype(np.int32)) @ ed_all[rr]
        nf = opt[viol[opt] == 0]
        if len(nf) == 0:
            C['C1_FAIL(no NF optimum)'] += 1; fails.append(('c1', hyp, A.tolist(), H))
        if len(nf) < len(opt): C['instances_with_non_NF_optimum(tie)'] += 1
        if (viol[opt] == 0).sum() == 0: pass
        # strict C1 weak form: some optimum has no R-R edit == viol==0
        # --- exchange formula and 2-class inequality over ALL partitions
        if len(classes) >= 2 and l4:
            for X1, X2 in itertools.permutations(classes, 2):
                A1 = [z for z in range(n) if z not in X1 and A[X1[0], z]]
                x0 = X1[0]
                M1 = (Lab == Lab[:, x0:x0 + 1])
                whole1 = M1[:, X1].all(1); has2 = M1[:, X2].all(1)
                rows = np.nonzero(whole1 & has2)[0]
                if len(rows) == 0: continue
                Csz = M1[rows].sum(1); a1 = M1[rows][:, A1].sum(1) if A1 else np.zeros(len(rows), int)
                d1 = len(X1) * (len(X1) + 2 * a1 - Csz)
                # actual cost change
                if len(rows) > 400: sel = np.random.RandomState(seed + it).choice(len(rows), 400, replace=False)
                else: sel = np.arange(len(rows))
                for t in sel:
                    lab = Lab[rows[t]].astype(int).copy(); c0 = costs[rows[t]]
                    lab2 = lab.copy()
                    for x in X1: lab2[x] = 99
                    C['exchange_formula_tests'] += 1
                    if L.cost_labels(A, lab2) - c0 != d1[t]: C['EXCHANGE_FORMULA_FAIL'] += 1; fails.append(('xf', A.tolist(), H, lab.tolist(), X1, X2))
                # the inequality needs both directions
                if X1 < X2:
                    M2 = (Lab == Lab[:, X2[0]:X2[0] + 1])
                    A2 = [z for z in range(n) if z not in X2 and A[X2[0], z]]
                    a2 = M2[rows][:, A2].sum(1) if A2 else np.zeros(len(rows), int)
                    d2 = len(X2) * (len(X2) + 2 * a2 - Csz)
                    bad = ((d1 > 0) & (d2 > 0)).sum()
                    C['two_class_partitions(all, not only optimal)'] += len(rows)
                    if bad: C['INEQUALITY_FAIL'] += int(bad); fails.append(('ineq', A.tolist(), H))
        # --- move DAG from every (capped) optimum
        if len(classes) >= 2 and l4:
            cls_of = {}
            for ci, X in enumerate(classes):
                for x in X: cls_of[x] = ci
            def shared(lab):
                cnt = 0
                lc = [lab[X[0]] for X in classes]
                for a, b in itertools.combinations(range(len(classes)), 2):
                    if lc[a] == lc[b]: cnt += 1
                return cnt
            sample = opt if len(opt) <= 80 else opt[np.random.RandomState(it).choice(len(opt), 80, replace=False)]
            for oi in sample:
                start = tuple(int(x) for x in Lab[oi])
                stack = [start]; seen = {canon(start)}; nstates = 0
                while stack and nstates < 1500:
                    lab = stack.pop(); nstates += 1
                    sh = shared(lab)
                    if sh == 0:
                        C['leaf_NF'] += 1
                        if int(costs[[i for i in [0]][0]]) is None: pass
                        continue
                    found = False
                    for X in classes:
                        lx = lab[X[0]]
                        if sum(1 for X2 in classes if lab[X2[0]] == lx) < 2: continue
                        lab2 = list(lab)
                        for x in X: lab2[x] = 99
                        if L.cost_labels(A, lab2) <= OPT:
                            found = True
                            assert L.cost_labels(A, lab2) == OPT
                            if shared(lab2) >= sh: C['POTENTIAL_NOT_DECREASING'] += 1
                            cl = canon(lab2)
                            if cl not in seen: seen.add(cl); stack.append(lab2)
                    if not found: C['DEAD_END(non-NF optimum with no free move)'] += 1; fails.append(('dead', A.tolist(), H))
                C['dag_starts'] += 1
        # --- statements on NF optima and on all optima
        cls_info = []
        for X in classes:
            AX = [z for z in range(n) if z not in X and A[X[0], z]]
            cls_info.append((X, AX))
        qmax = max(len(Q[v]) for v in range(n))
        # N4 (graph only)
        if l4 and pm:
            for X, AX in cls_info:
                assert all(z in VH for z in AX)
                for z, z2 in itertools.combinations(AX, 2):
                    if not A[z, z2] and z2 not in Q[z]: C['N4_nonedge_in_A_not_packed'] += 1
                for a in VH:
                    if a in AX: continue
                    for z in AX:
                        if A[a, z] and a not in Q[z]: C['N4_edge_not_packed'] += 1
                        C['N4_checks'] += 1
            adjcl = collections.Counter()
            for X, AX in cls_info:
                for a in AX: adjcl[a] += 1
            if any(v > 1 for v in adjcl.values()): C['N4_vertex_two_classes'] += 1
        # N3 on EVERY optimum
        for X, AX in cls_info:
            if len(X) > lopt:
                C['N3_classes'] += 1
                target = set(X) | set(AX)
                for oi in opt:
                    lab = Lab[oi]; cl = {v for v in range(n) if lab[v] == lab[X[0]]}
                    if cl != target: C['N3_FAIL'] += 1; fails.append(('N3', A.tolist(), H, X)); break
        # P2, P3, P4 on NF optima
        for oi in nf:
            lab = Lab[oi]
            for X, AX in cls_info:
                cl = {v for v in range(n) if lab[v] == lab[X[0]]}
                B = cl & VH
                if B and not (B & set(AX)):
                    C['P2_FAIL'] += 1; fails.append(('P2', A.tolist(), H, X))
                C['P2_tests'] += 1
                for a in B - set(AX):
                    zs = B & set(AX)
                    if not zs: C['P3_FAIL(no z)'] += 1
                    for z in zs:
                        ok = (a in Q[z]) or ((not A[z, a]) and True)
                        if A[z, a] and a not in Q[z]: C['P3_FAIL'] += 1
                        if (a not in Q[z]) and not (not A[z, a]): C['P3_FAIL2'] += 1
                        C['P3_tests'] += 1
                if len(AX) >= lopt + qmax:
                    C['P4_tests'] += 1
                    if B - set(AX): C['P4_FAIL'] += 1; fails.append(('P4', A.tolist(), H, X))
        # N5 shrink: deleting one x from class X with |X|>=2 ; only claim if OPT(G-x)-|H| <= |X|-2
        for X, AX in cls_info:
            if len(X) >= 2:
                x = X[-1]; keep = [v for v in range(n) if v != x]
                A2 = A[np.ix_(keep, keep)]
                O2 = int(L.costs_all(A2).min())
                C['N5_tests'] += 1
                if O2 - len(H) <= len(X) - 2:
                    C['N5_applicable'] += 1
                    if O2 != OPT: C['N5_FAIL'] += 1; fails.append(('N5', A.tolist(), H, X))
                else:
                    C['N5_not_applicable'] += 1
                if O2 != OPT: C['OPT_changes_on_class_deletion(nonapplicable)'] += 1
    return C, fails

if __name__ == '__main__':
    seed, ninst, nmin, nmax = map(int, sys.argv[1:5]); hyp = sys.argv[5] if len(sys.argv) > 5 else 'pairmax'
    C, fails = run(seed, ninst, nmin, nmax, hyp)
    print('seed', seed, 'hyp', hyp, dict(sorted(C.items())))
    print('n_fail_records', len(fails))
    for f in fails[:2]: print(f)
