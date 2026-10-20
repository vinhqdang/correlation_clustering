"""Property OUT: exists max packing H* and optimum S such that every hard edit of S has an endpoint outside V(H*).
Also property IN0: exists max packing H* and optimum S with NO hard edit at all (i.e. OPT=nu with ... only if g=0 else inherent gap in-member)."""
import random, sys, collections, time
from w8a_lib import *
def run(g, maxpk=4000):
    opt, sols = g.opt(limit=5000)
    parts = [edits_of_partition(g, part_from_mask(s, g.n)) for s in sols]
    nu, _ = nu_mod(g)
    packs = enum_packings(g, minsize=nu, maxcount=maxpk)
    bestboth = None; nooutfail = None; cnt_ok = 0
    for H in packs:
        H = list(H); V = set(v for t in H for v in t)
        mb = None
        for S in parts:
            he = hard_edits(H, S)
            nb = sum(1 for (a, b) in he if a in V and b in V)
            if mb is None or nb < mb: mb = nb
            if mb == 0: break
        if mb == 0: cnt_ok += 1
        if bestboth is None or mb < bestboth: bestboth = mb
    return opt, nu, len(packs), bestboth, cnt_ok
if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2]); random.seed(seed)
    stat = collections.Counter(); bad = []
    t0 = time.time()
    for it in range(N):
        n = random.randint(7, 12)
        mode = random.choice(['gnp', 'tree+chords', 'cycle+chords'])
        if mode == 'gnp':
            p = random.choice([0.2, 0.3, 0.4]); E = [(a, b) for a in range(n) for b in range(a + 1, n) if random.random() < p]
        elif mode == 'tree+chords':
            E = [(random.randrange(i), i) for i in range(1, n)]
            for _ in range(random.randint(1, 4)):
                a, b = random.sample(range(n), 2); E.append((a, b))
        else:
            E = [(i, (i + 1) % n) for i in range(n)]
            for _ in range(random.randint(0, 3)):
                a, b = random.sample(range(n), 2); E.append((a, b))
        g = G8(n, E)
        if not nx.is_connected(nx.Graph(list(g.E) + [(0, 0)])) or len(nx.Graph(list(g.E)).nodes) < n: continue
        opt, nu, npk, mb, ok = run(g)
        stat[(opt - nu, mb)] += 1
        if mb and mb > 0: bad.append((n, sorted(g.E), opt, nu, mb))
    print('seed', seed, 'time %.0fs' % (time.time() - t0))
    for k, v in sorted(stat.items()): print('  (OPT-nu, min #hard edits with both ends in V(H*)) =', k, v)
    for b in bad[:10]: print('  BAD', b)
