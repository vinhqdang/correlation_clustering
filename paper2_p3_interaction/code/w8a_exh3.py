"""NF** (A-local covered edits + hub condition for edits at uncovered vertices) : exhaustive over connected graphs n (geng res/mod), only graphs with OPT>nu."""
import sys, subprocess, collections, time
import networkx as nx
from w8a_lib import *
from w8a_ilp3 import joint3
n = int(sys.argv[1]); res = int(sys.argv[2]); mod = int(sys.argv[3])
p = subprocess.Popen(['nauty-geng', '-c', '-q', str(n), '%d/%d' % (res, mod)], stdout=subprocess.PIPE)
st = collections.Counter(); bads = []; t0 = time.time(); cnt = 0
for line in p.stdout:
    G = nx.from_graph6_bytes(line.strip()); g = G8(n, list(G.edges()))
    if not g.p3s(): continue
    cnt += 1
    nu, _ = nu_mod(g); opt, _ = opt_ilp(g)
    if opt == nu: st['tight'] += 1; continue
    v = joint3(g, opt, nu, 2, tl=60)
    st[('gap>0', 'NF** bad', v)] += 1
    if v: bads.append((sorted(g.E), opt, nu, v))
print('n', n, 'res', res, 'mod', mod, 'graphs', cnt, 'time %.0fs' % (time.time() - t0))
for k, v in sorted(st.items(), key=str): print('  ', k, v)
for b in bads[:5]: print('  BAD', b)
