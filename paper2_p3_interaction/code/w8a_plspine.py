import random, networkx as nx
from w8a_lib import *
from w8a_pl import class_opt, decomposition, uloc
from w8a_spine import build
rng = random.Random(7)
for name, F in (('triangle', nx.complete_graph(3)), ('K4', nx.complete_graph(4))):
    n, E, H, idx, y, K, arcs = build(F, rng)
    g = G8(n, list(E)); opt, _ = opt_ilp(g)
    bags = decomposition(H); cand = uloc(g, H, bags)
    # restrict candidate set to pairs incident to the trigger for speed? no: full class
    import time; t = time.time()
    co, U, cand = class_opt(g, H, 1, bags)
    print('Sp(%s): n=%d |H|=%d OPT=%d CLASS-OPT=%s best U=%s |U_loc|=%d (%.0fs); y=%d' % (name, n, len(H), opt, co, U, len(cand), time.time() - t, y))
