"""For random (G,H) whose promise (Thm PL class w.r.t. min-fill decomposition) FAILS: is there a packing H* (|H*| >= |H|, pair-maximal) for which the promise holds,
and what is the minimum of tw(I(H*)) (min-fill upper bound) - tw(I(H)) over such H*, compared with l ?"""
import random, sys, collections, time
from w8a_lib import *
from w8a_pl import decomposition, uloc
from networkx.algorithms.approximation import treewidth_min_fill_in
def tw_ub(H):
    I = Igraph(H)
    return 0 if I.number_of_edges() == 0 else treewidth_min_fill_in(I)[0]
rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2]); st = collections.Counter(); ex = []; t0 = time.time()
for it in range(N):
    n = rng.randint(8, 12); E = [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < rng.choice([0.15, 0.2, 0.25])]
    g = G8(n, E)
    gg = nx.Graph(list(g.E))
    if not g.p3s() or len(g.E) < 3 or gg.number_of_nodes() < n or not nx.is_connected(gg): continue
    H = []; used = set(); P = g.p3s()[:]; rng.shuffle(P)
    for t in P:
        if not any(p in used for p in mp(t)): H.append(t); used.update(mp(t))
    opt, sols = g.opt(limit=500); l = opt - len(H)
    optS = [edits_of_partition(g, part_from_mask(s, n)) for s in sols]
    def promise(Hx):
        cand = set(uloc(g, Hx, decomposition(Hx)))
        return any(set(hard_edits(Hx, S)) <= cand for S in optS)
    if promise(H): continue
    w0 = tw_ub(H); best = None
    for Hx in enum_packings(g, minsize=len(H), maxcount=100000):
        Hx = list(Hx)
        if not pair_maximal(g, Hx): continue
        if promise(Hx):
            w = tw_ub(Hx)
            if best is None or w < best: best = w
    st[('l=%d' % min(l, 3), 'tw(I(H))=%d' % w0, 'min tw over promise-packings = %s' % best)] += 1
    if best is None or best > w0 + l: ex.append((n, sorted(g.E), H, l, w0, best))
print('seed', sys.argv[1], 'time %.0fs' % (time.time() - t0))
for k, v in sorted(st.items(), key=str): print('  ', k, v)
print('  #promise-failing instances with no promise-packing of width <= tw(I(H)) + l :', len(ex))
for e in ex[:3]: print('  EX', e)
