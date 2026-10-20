"""Hill-climb for graphs where every maximum packing + optimum leaves >=1 hard edit (both ends covered) at soft distance > r."""
import sys, random, time
from w8a_lib import *
from w8a_ilp3 import joint3
from w8a_ilp import rand_graph
def score(E, n, r):
    g = G8(n, E)
    gg = nx.Graph(list(g.E))
    if gg.number_of_nodes() < n or not nx.is_connected(gg) or len(g.E) < 2: return None
    opt, _ = opt_ilp(g); nu, _ = nu_mod(g)
    v = joint3(g, opt, nu, r, tl=30)
    return None if v is None else (v, opt - nu)
if __name__ == '__main__':
    seed = int(sys.argv[1]); n = int(sys.argv[2]); r = int(sys.argv[3]); T = float(sys.argv[4])
    rng = random.Random(seed); t0 = time.time(); best = None; bestE = None
    while time.time() - t0 < T:
        E = rand_graph(n, rng.choice(['gnp', 'tree+chords', 'cycle+chords', 'cubic']), rng)
        s = score(E, n, r)
        if s is None: continue
        cur = (s[0], -s[1]) if s else None   # maximise bad, tie-break: small gap
        stale = 0
        while stale < 25 and time.time() - t0 < T:
            E2 = set(key(a, b) for a, b in E)
            a, b = rng.sample(range(n), 2); p = key(a, b)
            if p in E2: E2.remove(p)
            else: E2.add(p)
            s2 = score(sorted(E2), n, r)
            if s2 is None: stale += 1; continue
            c2 = (s2[0], -s2[1])
            if c2 >= cur:
                if c2 > cur: stale = 0
                else: stale += 1
                cur, E = c2, sorted(E2)
            else: stale += 1
        if best is None or cur > best: best, bestE = cur, E; print('seed', seed, 'n', n, 'r', r, 'best (bad,-gap)', best, 'E', bestE, flush=True)
