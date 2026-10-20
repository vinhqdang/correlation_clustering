"""NF* test on symmetric / cyclic families (displacement chains could cycle): circulants, prisms, Mobius ladders, wheels, cube-like, Petersen, subdivided graphs."""
import itertools, collections, time
import networkx as nx
from w8a_lib import *
from w8a_ilp2 import joint2
fams = []
for n in range(7, 17):
    for ks in [(1, 2), (1, 3), (1, 4), (2, 3), (1, 2, 4), (1, 5)]:
        if max(ks) < n / 2: fams.append(('circ%d%s' % (n, ks), nx.circulant_graph(n, list(ks))))
for L in range(3, 9):
    fams.append(('prism%d' % L, nx.circular_ladder_graph(L)))
    fams.append(('mobius%d' % L, nx.Graph([(i, (i + 1) % (2 * L)) for i in range(2 * L)] + [(i, i + L) for i in range(L)])))
    fams.append(('wheel%d' % (L + 3), nx.wheel_graph(L + 3)))
fams.append(('petersen', nx.petersen_graph()))
fams.append(('cube', nx.hypercube_graph(3)))
fams.append(('dodeca', nx.dodecahedral_graph()))
for k in range(3, 8): fams.append(('subdiv-K4-%d' % k, None))
out = collections.Counter(); bad = []
t0 = time.time()
for name, G in fams:
    if G is None: continue
    G = nx.convert_node_labels_to_integers(G)
    n = G.number_of_nodes()
    if n > 20: continue
    g = G8(n, list(G.edges()))
    opt, _ = opt_ilp(g); nu, _ = nu_mod(g)
    if opt == nu: out['tight'] += 1; continue
    v2 = joint2(g, opt, nu, 2, tl=60)
    out[('gap', opt - nu, 'bad(r=2)', v2)] += 1
    if v2: bad.append((name, v2))
    print(name, 'n', n, 'OPT', opt, 'nu', nu, 'bad(r=2)', v2, flush=True)
print(dict(out), 'time %.0fs' % (time.time() - t0)); print('BAD', bad)
