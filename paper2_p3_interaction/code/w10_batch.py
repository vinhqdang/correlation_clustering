"""Batch test of the thin guest model: ZE(G - cuts(choice)) feasible  <=>  choice is a Grid-Tiling(<=) solution."""
import sys, random, itertools, time
from w10_red import *
def gt_all(k, n, S):
    return gt_solvable(k, n, S)
def run(k, n, seed, dens, nsamp):
    rng = random.Random(seed)
    S = [[set((a, b) for a in range(1, n+1) for b in range(1, n+1) if rng.random() < dens) for _ in range(k)] for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if not S[i][j]: S[i][j].add((rng.randint(1, n), rng.randint(1, n)))
    cells = [(i, j) for i in range(k) for j in range(k)]
    W, meta = build_gt(k, n, S, 2, thin_mode=True)
    I = ThinInst(W)
    sols = gt_all(k, n, S) if k*k <= 6 or True else []
    tests = []
    solkeys = list(sols)
    rng.shuffle(solkeys)
    tests += [(t, True) for t in solkeys[:nsamp]]
    # invalid: random choices that are not solutions (also "near misses": one cell changed from a solution)
    tries = 0; inval = []
    while len(inval) < nsamp and tries < 20 * nsamp:
        tries += 1
        if solkeys and rng.random() < 0.5:
            t = dict(rng.choice(solkeys)); c = rng.choice(cells); t[c] = rng.choice(sorted(S[c[0]][c[1]]))
        else:
            t = {c: rng.choice(sorted(S[c[0]][c[1]])) for c in cells}
        if t not in sols and t not in inval: inval.append(t)
    tests += [(t, False) for t in inval]
    bad = 0
    for t, exp in tests:
        lab = I.solve(drop_hard=cut_choice(meta, t))
        if (lab is not None) != exp: bad += 1; print('MISMATCH', seed, t, exp, flush=True)
    return len(tests), len(solkeys), bad
if __name__ == '__main__':
    k, n, s0, s1, dens, nsamp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
    T = B = 0; t0 = time.time(); nsol = 0
    for seed in range(s0, s1):
        a, s, b = run(k, n, seed, dens, nsamp)
        T += a; B += b; nsol += (s > 0)
    print('k', k, 'n', n, 'seeds', s0, s1, 'tests', T, 'instances with a GT solution', nsol, 'mismatches', B, round(time.time()-t0, 1))
