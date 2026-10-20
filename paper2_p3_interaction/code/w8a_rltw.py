"""RL test: for random pair-maximal H with forced far hard edits (s(H)>0), min over packings H* (|H*|>=|H|, pair-maximal, s(H*)=0) of tw(I(H*)); compare with tw(I(H)), l."""
import random, sys, time, collections
from w8a_lib import *
from w8a_ilp import rand_graph
from w8a_kappa import far_count, random_pm_packing
if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2]); nlo, nhi = int(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed); st = collections.Counter(); ex = []; t0 = time.time()
    for it in range(N):
        n = rng.randint(nlo, nhi); mode = rng.choice(['gnp', 'tree+chords', 'cycle+chords', 'cubic', 'planted'])
        E = rand_graph(n, mode, rng); g = G8(n, E); gg = nx.Graph(list(g.E))
        if len(g.E) < 3 or gg.number_of_nodes() < n or not nx.is_connected(gg): continue
        H = random_pm_packing(g, rng)
        opt, sols = g.opt(limit=3000)
        optima = [edits_of_partition(g, part_from_mask(s, n)) for s in sols]
        if far_count(H, optima) == 0: continue
        w0 = tw_exact(Igraph(H)); l = opt - len(H)
        packs = enum_packings(g, minsize=len(H), maxcount=100000)
        best = None
        for P in packs:
            P = list(P)
            if not pair_maximal(g, P) or far_count(P, optima) != 0: continue
            w = tw_exact(Igraph(P))
            if best is None or w < best: best = w
        st[('l=%d' % min(l, 4), 'tw(I(H))=%d' % w0, 'min tw(I(H*))=%s' % best)] += 1
        if best is None or best > w0 + l: ex.append((n, sorted(g.E), H, l, w0, best))
    print('seed', seed, 'time %.0fs' % (time.time() - t0))
    for k, v in sorted(st.items(), key=str): print('  ', k, v)
    print('  #instances with min tw(I(H*)) > tw(I(H)) + l :', len(ex))
    for e in ex[:4]: print('  EX', e)
