import sys, random, itertools
from a15_lib import *
def run(k, n, S, p, label):
    t0 = time.time()
    W, meta = build_gt(k, n, S, p, thin_mode=False)
    I = Inst(W)
    sols = gt_solvable(k, n, S)
    N, size, e, ne, cid = contract(I)
    cost, lab = wopt(N, size, e, ne, verbose=False)
    ex = cost - len(I.members)
    print(label, 'S', [[sorted(x) for x in r] for r in S], 'GT-solutions', len(sols), '| n', I.n, '|H|', len(I.members), 'valid', I.valid[0], 'c', I.c, 'u_max', max(W.used.values()), '| classes', N, 'OPT excess', ex, 'k^2', k * k,
          'OK' if ((ex <= k * k) == (len(sols) > 0)) and (len(sols) == 0 or ex == k * k) else 'MISMATCH', round(time.time() - t0, 1), flush=True)
    return ex, len(sols)
if __name__ == '__main__':
    k = int(sys.argv[1]); n = int(sys.argv[2]); p = int(sys.argv[3])
    which = sys.argv[4]
    if which == 'fixed':
        Y = [[{(1, 1)}, {(1, 1)}], [{(1, 1)}, {(1, 1)}]]
        N1 = [[{(2, 1)}, {(1, 1)}], [{(1, 1)}, {(1, 1)}]]   # a11=2 > a12=1
        N2 = [[{(1, 2)}, {(1, 1)}], [{(1, 1)}, {(1, 1)}]]   # b11=2 > b21=1
        for nm, S in (('YES', Y), ('NO-a', N1), ('NO-b', N2)):
            run(k, n, S, p, nm)
    else:
        seed = int(which)
        rng = random.Random(seed)
        for _ in range(int(sys.argv[5])):
            dens = rng.choice([0.3, 0.5])
            S = [[set((a, b) for a in range(1, n + 1) for b in range(1, n + 1) if rng.random() < dens) or {(rng.randint(1, n), rng.randint(1, n))} for _ in range(k)] for _ in range(k)]
            run(k, n, S, p, 'rnd')
