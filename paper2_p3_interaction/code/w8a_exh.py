"""Exhaustive over connected graphs from nauty-geng (n given, res/mod split).  For graphs with OPT>nu: min bad at r=0 (OUT) [cheap],
and for those with bad0>0 also at r=2 and r=3 (soft distance)."""
import sys, subprocess, collections, time
import networkx as nx
from w8a_lib import *
from w8a_ilp import joint
from w8a_ilp2 import joint2
n = int(sys.argv[1]); res = int(sys.argv[2]); mod = int(sys.argv[3])
p = subprocess.Popen(['nauty-geng', '-c', '-q', str(n), '%d/%d' % (res, mod)], stdout=subprocess.PIPE)
st = collections.Counter(); bads = []; t0 = time.time(); cnt = 0
for line in p.stdout:
    G = nx.from_graph6_bytes(line.strip())
    g = G8(n, list(G.edges()))
    if not g.p3s(): continue
    cnt += 1
    nu, _ = nu_mod(g)
    # cheap upper bound filter: if OPT == nu then nothing to do
    opt, _ = opt_ilp(g)
    if opt == nu: st['tight'] += 1; continue
    r0 = joint(g, opt, nu)[0]
    if r0 == 0: st['gap>0, OUT ok'] += 1; continue
    r2 = joint2(g, opt, nu, 2)
    r3 = joint2(g, opt, nu, 3) if r2 else 0
    st[('OUT fails; r2 bad', r2, 'r3 bad', r3)] += 1
    if r2: bads.append((sorted(g.E), opt, nu, r2, r3))
print('n', n, 'res', res, 'mod', mod, 'graphs', cnt, 'time %.0fs' % (time.time() - t0))
for k, v in sorted(st.items(), key=str): print('  ', k, v)
for b in bads[:5]: print('  BAD', b)
