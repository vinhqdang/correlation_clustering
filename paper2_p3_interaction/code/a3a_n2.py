"""N2 identity |U| = #HH + sum_X |X||A_X xor B_X| on every optimum, plus generator coverage statistics (exhaustive, n<=9)."""
import random, itertools, collections, numpy as np, sys
import a3a_lib as L
rng = random.Random(int(sys.argv[1])); C = collections.Counter()
for it in range(int(sys.argv[2])):
    n = rng.randint(7, 9); f, A, H = L.make_instance(rng, n)
    VH, R, Q, cls = L.info(A, H)
    if len(R) < 2: continue
    C['inst'] += 1; C['ge2classes'] += len(cls) >= 2; C['bigclass(>=2)'] += any(len(X) >= 2 for X in cls)
    C['ge2classes_nonempty_A'] += sum(1 for X in cls if any(A[X[0], z] for z in range(n) if z not in X)) >= 2
    C['maxclass=%d' % max(len(X) for X in cls)] += 1
    costs = L.costs_all(A); OPT = costs.min(); opt = np.nonzero(costs == OPT)[0]
    C['l=%d' % (OPT - len(H))] += 1
    packed = set()
    for t in H: packed.update(L.tri_pairs(t))
    for oi in opt:
        lab = L.rgs(n)[oi]; U = 0; hh = 0
        for i, j in itertools.combinations(range(n), 2):
            if (lab[i] == lab[j]) != bool(A[i, j]) and frozenset((i, j)) not in packed:
                U += 1; hh += (i in VH and j in VH)
        rhs = hh
        for X in cls:
            AX = {z for z in range(n) if z not in X and A[X[0], z]}
            B = {v for v in VH if lab[v] == lab[X[0]]}
            rhs += len(X) * len(AX ^ B)
        C['N2_tests'] += 1
        if rhs != U: C['N2_FAIL'] += 1
print(dict(sorted(C.items())))
