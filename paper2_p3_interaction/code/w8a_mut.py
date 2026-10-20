"""Skeleton-mutation search for NF* violations (bad edits at soft distance > r) around FAR_m and spine-like skeletons, n ~ 13..19."""
import sys, random, time
from w8a_lib import *
from w8a_ilp2 import joint2
def far(m): return 2 * m + 1, [(i, i + 1) for i in range(2 * m)] + [(1, 2 * m - 1)]
def farchords(m, rng):
    n, E = far(m); E = list(E)
    for _ in range(rng.randint(0, 2)):
        a, b = rng.sample(range(n), 2); E.append((a, b))
    return n, E
def score(n, E, r):
    g = G8(n, E); gg = nx.Graph(list(g.E))
    if gg.number_of_nodes() < n or not nx.is_connected(gg) or len(g.E) < 2: return None
    opt, _ = opt_ilp(g); nu, _ = nu_mod(g)
    if opt == nu: return (0, 0)
    v = joint2(g, opt, nu, r, tl=40)
    return None if v is None else (v, opt - nu)
if __name__ == '__main__':
    seed = int(sys.argv[1]); r = int(sys.argv[2]); T = float(sys.argv[3]); rng = random.Random(seed); t0 = time.time()
    best = (-1, 0); cnt = 0; hist = {}
    while time.time() - t0 < T:
        m = rng.randint(4, 8); n, E = farchords(m, rng)
        cur = score(n, E, r)
        if cur is None: continue
        E = sorted(set(key(a, b) for a, b in E)); stale = 0
        while stale < 15 and time.time() - t0 < T:
            E2 = set(E)
            for _ in range(rng.randint(1, 2)):
                a, b = rng.sample(range(n), 2); p = key(a, b)
                if p in E2 and rng.random() < .5: E2.discard(p)
                else: E2.add(p)
            s = score(n, sorted(E2), r); cnt += 1
            if s is None: stale += 1; continue
            hist[s[0]] = hist.get(s[0], 0) + 1
            if s[0] >= cur[0]:
                if s[0] > cur[0]: stale = 0
                else: stale += 1
                cur, E = s, sorted(E2)
            else: stale += 1
            if s[0] > best[0]: best = s; print('seed', seed, 'NEW BEST', s, 'n', n, 'E', sorted(E2), flush=True)
    print('seed', seed, 'r', r, 'evaluated', cnt, 'hist of bad', hist)
