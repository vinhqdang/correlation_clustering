"""NF* test (r=2) on structured families at n up to ~21 via joint ILP: paths/cycles with 1-3 chords, ladders, grids, wheels-like, caterpillars."""
import sys, random, time, itertools, collections
import networkx as nx
from w8a_lib import *
from w8a_ilp3 import joint3
def test(name, n, E, stat, bad, r=2):
    g = G8(n, E)
    gg = nx.Graph(list(g.E))
    if gg.number_of_nodes() < n or not nx.is_connected(gg): return
    opt, _ = opt_ilp(g); nu, _ = nu_mod(g)
    if opt == nu: stat[(name, 'tight')] += 1; return
    v = joint3(g, opt, nu, r, tl=60)
    stat[(name, 'gap=%d' % (opt - nu), 'bad=%s' % v)] += 1
    if v: bad.append((name, n, sorted(g.E), opt, nu, v))
if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2]); rng = random.Random(seed); stat = collections.Counter(); bad = []; t0 = time.time()
    for it in range(N):
        fam = rng.choice(['path+chords', 'cycle+chords', 'ladder', 'caterpillar+chords', 'tree+chords', 'planted', 'planted', 'gnp'])
        if fam == 'path+chords':
            n = rng.randint(9, 17); E = [(i, i + 1) for i in range(n - 1)]; k = rng.randint(1, 3)
        elif fam == 'cycle+chords':
            n = rng.randint(9, 17); E = [(i, (i + 1) % n) for i in range(n)]; k = rng.randint(1, 3)
        elif fam == 'ladder':
            L = rng.randint(4, 8); n = 2 * L; E = [(i, i + 1) for i in range(L - 1)] + [(L + i, L + i + 1) for i in range(L - 1)] + [(i, L + i) for i in range(L)]; k = rng.randint(0, 2)
        elif fam == 'grid':
            a, b = rng.choice([(3, 4), (3, 5), (4, 4), (2, 7)]); n = a * b
            E = [(i * b + j, i * b + j + 1) for i in range(a) for j in range(b - 1)] + [(i * b + j, (i + 1) * b + j) for i in range(a - 1) for j in range(b)]; k = rng.randint(0, 2)
        elif fam == 'caterpillar+chords':
            L = rng.randint(4, 8); n = 2 * L; E = [(i, i + 1) for i in range(L - 1)] + [(i, L + i) for i in range(L)]; k = rng.randint(1, 3)
        elif fam == 'planted':
            n = rng.randint(9, 14); kk = rng.randint(2, 4); lab = [rng.randrange(kk) for _ in range(n)]; q = rng.choice([0.08, 0.12, 0.18])
            E = [(a, b) for a in range(n) for b in range(a + 1, n) if (lab[a] == lab[b]) != (rng.random() < q)]; k = 0
        elif fam == 'gnp':
            n = rng.randint(9, 13); pp = rng.choice([0.2, 0.3, 0.4]); E = [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < pp]; k = 0
        else:
            n = rng.randint(10, 16); E = [(rng.randrange(i), i) for i in range(1, n)]; k = rng.randint(1, 4)
        for _ in range(k):
            a, b = rng.sample(range(n), 2); E.append((a, b))
        test(fam, n, E, stat, bad)
    print('seed', seed, 'time %.0fs' % (time.time() - t0))
    for k, v in sorted(stat.items(), key=str): print('  ', k, v)
    for b in bad[:5]: print('  BAD', b)
