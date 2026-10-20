"""a13_t2_cert.py: Theorem CE.  F'(T_N), F'(T_N - e), F'(empty) by CP-SAT (pair-variable formulation, written fresh);
partition-variant F (exact subset DP) for N<=2;  sanity check of the DP against CP-SAT on random graphs."""
import sys, random, itertools
from a13_lib import *
# sanity: DP vs CP-SAT optimum on random graphs (all pairs counted, no T)
random.seed(1)
for it in range(30):
    n = random.randint(5, 10)
    E = set(frozenset(p) for p in itertools.combinations(range(n), 2) if random.random() < .4)
    inst = dict(n=n, E=E, members=[], T=[])
    allp = set(frozenset(p) for p in itertools.combinations(range(n), 2))
    a = dp_partition(inst, allp, []); b = Fprime(inst, [], extra_counted=allp)
    assert a == b, (a, b)
print('DP == CP-SAT on 30 random graphs (OPT)')
# sanity with T: partition variant must be >= csp variant
for it in range(30):
    n = random.randint(5, 9)
    E = set(frozenset(p) for p in itertools.combinations(range(n), 2) if random.random() < .4)
    inst = dict(n=n, E=E, members=[], T=[])
    allp = [frozenset(p) for p in itertools.combinations(range(n), 2)]
    cnt = set(p for p in allp if random.random() < .6)
    T = [p for p in allp if p not in cnt and random.random() < .5]
    a = dp_partition(inst, cnt, T); b = Fprime(inst, T, extra_counted=cnt)
    if b is None: continue
    assert a >= b, (a, b)
print('F_partition >= F_csp on 30 random (counted,T) splits')
Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
for N in range(1, Nmax+1):
    inst = build_GN(N); H = len(inst['members']); T = inst['T']
    far = set(frozenset(p) for p in itertools.combinations(range(inst['n']), 2)) - counted_pairs(inst)
    assert all(t in far for t in T)
    f0 = Fprime(inst, []) - H
    fT = Fprime(inst, T) - H
    sub = [Fprime(inst, [t for t in T if t != e]) - H for e in T]
    print('N=%d |T|=%d  F\'(empty)-|H|=%d  F\'(T)-|H|=%d  max_e F\'(T-e)-|H|=%d (all==0: %s)' % (N, len(T), f0, fT, max(sub), all(s == 0 for s in sub)), flush=True)
    if N <= 2:
        cnt = counted_pairs(inst)
        pf = dp_partition(inst, cnt, T) - H
        ps = [dp_partition(inst, cnt, [t for t in T if t != e]) - H for e in T]
        print('    partition variant (exact DP): F(T)-|H|=%d, F(T-e)-|H| = %s' % (pf, ps), flush=True)
