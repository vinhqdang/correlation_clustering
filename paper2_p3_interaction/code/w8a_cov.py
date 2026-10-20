"""Among maximum packings: relation between coverage |V(H*)| and property OUT (min over optima of #hard edits with both ends in V(H*))."""
import random, sys, collections, time
from w8a_lib import *
from w8a_ilp import rand_graph
def analyze(g, maxpk=3000):
    opt, sols = g.opt(limit=5000)
    parts = [edits_of_partition(g, part_from_mask(s, g.n)) for s in sols]
    nu, _ = nu_mod(g)
    packs = enum_packings(g, minsize=nu, maxcount=maxpk)
    rows = []
    for H in packs:
        H = list(H); V = set(v for t in H for v in t)
        mb = min(sum(1 for (a, b) in hard_edits(H, S) if a in V and b in V) for S in parts)
        rows.append((len(V), mb))
    return opt, nu, rows
if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2]); rng = random.Random(seed)
    st = collections.Counter(); fails = []
    for it in range(N):
        n = rng.randint(7, 11); mode = rng.choice(['gnp', 'tree+chords', 'cycle+chords', 'cubic', 'planted'])
        E = rand_graph(n, mode, rng); g = G8(n, E)
        gg = nx.Graph(list(g.E))
        if len(g.E) < 2 or gg.number_of_nodes() < n or not nx.is_connected(gg): continue
        opt, nu, rows = analyze(g)
        mincov = min(r[0] for r in rows); maxcov = max(r[0] for r in rows)
        okmin = all(mb == 0 for c, mb in rows if c == mincov)
        okany = any(mb == 0 for c, mb in rows)
        allok = all(mb == 0 for c, mb in rows)
        st[('minimal-coverage packings all OUT-good', okmin)] += 1
        st[('some packing OUT-good', okany)] += 1
        st[('ALL max packings OUT-good', allok)] += 1
        if not okmin: fails.append((n, sorted(g.E), rows))
    for k, v in sorted(st.items(), key=str): print(k, v)
    for f in fails[:3]: print('FAIL', f)
