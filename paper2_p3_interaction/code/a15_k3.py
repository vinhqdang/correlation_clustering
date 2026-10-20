"""a15: contracted exact OPT for k=3 (lean: no validity/c computation, which is O(n |H|))."""
import sys, random
from a15_lib import *


class Lean:
    def __init__(self, W):
        self.W = W
        self.n, self.E, self.members = W.finish()


def run(k, n, S, p, label):
    t0 = time.time()
    W, meta = build_gt(k, n, S, p, thin_mode=False)
    I = Lean(W)
    sols = gt_solvable(k, n, S)
    N, size, e, ne, cid = contract(I)
    cost, lab = wopt(N, size, e, ne, verbose=False)
    ex = cost - len(I.members)
    ok = ((ex <= k * k) == (len(sols) > 0)) and (len(sols) == 0 or ex == k * k)
    print(label, 'S', [[sorted(x) for x in r] for r in S], 'GT-solutions', len(sols), '| n', I.n, '|H|', len(I.members),
          'u_max', max(W.used.values()), '| classes', N, 'OPT excess', ex, 'k^2', k * k, 'OK' if ok else 'MISMATCH', round(time.time() - t0, 1), flush=True)


if __name__ == '__main__':
    k, n, p, seed, reps = [int(x) for x in sys.argv[1:6]]
    rng = random.Random(seed)
    for r in range(reps):
        if r == 0:
            S = [[{(1, 1)} for _ in range(k)] for _ in range(k)]
        elif r == 1:
            S = [[{(1, 1)} for _ in range(k)] for _ in range(k)]
            S[1][1] = {(2, 1)}          # middle cell: a11=2 > right neighbour's a=1  (NO)
        else:
            S = [[set((a, b) for a in range(1, n + 1) for b in range(1, n + 1) if rng.random() < 0.4)
                  or {(rng.randint(1, n), rng.randint(1, n))} for _ in range(k)] for _ in range(k)]
        run(k, n, S, p, 'k3')
