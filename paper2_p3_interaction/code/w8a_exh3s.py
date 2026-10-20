"""NF** on a uniform random sample of connected graphs on n vertices (graph6 lines from nauty-geng, reservoir-sampled)."""
import sys, subprocess, collections, time, random
import networkx as nx
from w8a_lib import *
from w8a_ilp3 import joint3
n = int(sys.argv[1]); K = int(sys.argv[2]); seed = int(sys.argv[3]); rng = random.Random(seed)
p = subprocess.Popen(['nauty-geng', '-c', '-q', str(n)], stdout=subprocess.PIPE)
res = []; cnt = 0
for line in p.stdout:
    cnt += 1
    if len(res) < K: res.append(line.strip())
    else:
        j = rng.randrange(cnt)
        if j < K: res[j] = line.strip()
st = collections.Counter(); bads = []; t0 = time.time()
for line in res:
    G = nx.from_graph6_bytes(line); g = G8(n, list(G.edges()))
    if not g.p3s(): continue
    nu, _ = nu_mod(g); opt, _ = opt_ilp(g)
    if opt == nu: st['tight'] += 1; continue
    v = joint3(g, opt, nu, 2, tl=60)
    st[('gap>0', 'NF** bad', v)] += 1
    if v: bads.append((sorted(g.E), opt, nu, v))
print('n', n, 'population', cnt, 'sample', len(res), 'seed', seed, 'time %.0fs' % (time.time() - t0))
for k, v in sorted(st.items(), key=str): print('  ', k, v)
for b in bads[:5]: print('  BAD', b)
