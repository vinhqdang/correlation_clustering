import collections, time
import networkx as nx
from networkx.generators.atlas import graph_atlas_g
from w8a_lib import *
from w8a_ilp import joint
st = collections.Counter(); bad = []
t0 = time.time()
for G in graph_atlas_g():
    n = G.number_of_nodes()
    if n < 3 or not nx.is_connected(G): continue
    g = G8(n, list(G.edges()))
    if not g.p3s(): continue
    opt, _ = opt_ilp(g); nu, _ = nu_mod(g)
    r = joint(g, opt, nu)
    st[(n, r[0])] += 1
    if r[0] > 0: bad.append((n, sorted(g.E), opt, nu, r))
print('time %.0fs' % (time.time() - t0))
for k, v in sorted(st.items()): print(k, v)
for b in bad[:10]: print('BAD', b)
