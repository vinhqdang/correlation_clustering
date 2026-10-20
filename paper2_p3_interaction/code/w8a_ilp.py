"""Joint ILP: choose a packing H* of size >= nu_target (default max) and a partition (optimal) minimising the number of hard edits with BOTH ends in V(H*).
Returns min 'bad' count, plus witness.  scipy milp / HiGHS."""
import random, sys, collections, time, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from w8a_lib import *

def joint(g, opt, nu, size=None, mode='both'):
    n = g.n; P = g.p3s(); size = nu if size is None else size
    pid = {}
    for a in range(n):
        for b in range(a + 1, n): pid[(a, b)] = len(pid)
    m = len(pid); nt = len(P)
    # variable layout: y[m] (together), x[nt] (packing), v[n] (covered), bad[m]
    oy, ox, ov, ob = 0, m, m + nt, m + nt + n
    N = ob + m
    rows = []; lo = []; hi = []
    def add(coefs, l, h):
        r = np.zeros(N)
        for i, c in coefs: r[i] += c
        rows.append(r); lo.append(l); hi.append(h)
    for a, b, d in itertools.combinations(range(n), 3):
        ab, bd, ad = pid[(a, b)], pid[(b, d)], pid[(a, d)]
        add([(ab, 1), (bd, 1), (ad, -1)], -np.inf, 1)
        add([(ab, 1), (ad, 1), (bd, -1)], -np.inf, 1)
        add([(bd, 1), (ad, 1), (ab, -1)], -np.inf, 1)
    # cost = |E| + sum_{nonedge} y - sum_{edge} y  == opt
    r = np.zeros(N)
    for p, i in pid.items(): r[oy + i] = -1 if p in g.E else 1
    rows.append(r); lo.append(opt - len(g.E)); hi.append(opt - len(g.E))
    # packing: pairs disjoint
    cover = {}
    for i, t in enumerate(P):
        for p in mp(t): cover.setdefault(p, []).append(i)
    for p, l in cover.items():
        add([(ox + i, 1) for i in l], -np.inf, 1)
    add([(ox + i, 1) for i in range(nt)], size, np.inf)
    for i, t in enumerate(P):
        for a in t: add([(ov + a, 1), (ox + i, -1)], 0, np.inf)
    for p, i in pid.items():
        a, b = p
        # z_p = edit indicator
        if p in g.E: zc = [(oy + i, -1)]; zconst = 1
        else: zc = [(oy + i, 1)]; zconst = 0
        soft = [(ox + j, 1) for j in cover.get(p, [])]
        # bad >= z + v_a + v_b - soft - 2
        coefs = [(ob + i, 1)] + [(k, -c) for k, c in zc] + [(ov + a, -1), (ov + b, -1)] + [(k, c) for k, c in soft]
        add(coefs, zconst - 2, np.inf)
    A = np.array(rows)
    c = np.zeros(N); c[ob:ob + m] = 1
    integ = np.ones(N); integ[ob:] = 0
    res = milp(c, constraints=LinearConstraint(A, lo, hi), integrality=integ, bounds=Bounds(0, 1), options={'time_limit': 120})
    if res.x is None: return None
    x = res.x
    H = [P[i] for i in range(nt) if x[ox + i] > .5]
    return int(round(res.fun)), H

def rand_graph(n, mode, rng):
    if mode == 'gnp':
        p = rng.choice([0.15, 0.2, 0.25, 0.3]); E = [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < p]
    elif mode == 'tree+chords':
        E = [(rng.randrange(i), i) for i in range(1, n)]
        for _ in range(rng.randint(1, 5)):
            a, b = rng.sample(range(n), 2); E.append((a, b))
    elif mode == 'cycle+chords':
        E = [(i, (i + 1) % n) for i in range(n)]
        for _ in range(rng.randint(0, 4)):
            a, b = rng.sample(range(n), 2); E.append((a, b))
    elif mode == 'planted':
        k = rng.randint(2, 4); lab = [rng.randrange(k) for _ in range(n)]; q = rng.choice([0.05, 0.1, 0.15, 0.2])
        E = []
        for a in range(n):
            for b in range(a + 1, n):
                t = lab[a] == lab[b]
                if rng.random() < q: t = not t
                if t: E.append((a, b))
    else:  # cubic-ish
        E = set()
        deg = [0] * n
        for _ in range(3 * n):
            a, b = rng.sample(range(n), 2)
            if deg[a] < 3 and deg[b] < 3 and (min(a, b), max(a, b)) not in E:
                E.add((min(a, b), max(a, b))); deg[a] += 1; deg[b] += 1
        E = list(E)
    return E

if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2]); nlo, nhi = int(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed); stat = collections.Counter(); bad = []
    t0 = time.time()
    for it in range(N):
        n = rng.randint(nlo, nhi); mode = rng.choice(['gnp', 'tree+chords', 'cycle+chords', 'cubic', 'planted', 'planted'])
        E = rand_graph(n, mode, rng); g = G8(n, E)
        if len(g.E) < 2: continue
        gg = nx.Graph(list(g.E))
        if gg.number_of_nodes() < n or not nx.is_connected(gg): continue
        opt, _ = opt_ilp(g); nu, _ = nu_mod(g)
        r = joint(g, opt, nu)
        if r is None: stat['timeout'] += 1; continue
        stat[(opt - nu, r[0])] += 1
        if r[0] > 0: bad.append((n, mode, sorted(g.E), opt, nu, r))
    print('seed', seed, 'n in', (nlo, nhi), 'time %.0fs' % (time.time() - t0))
    for k, v in sorted(stat.items(), key=str): print('  (OPT-nu, min bad) =', k, v)
    for b in bad[:6]: print('  BAD', b)
