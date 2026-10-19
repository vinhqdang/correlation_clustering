"""Search for a 'relaxable NEQ coupler': members B=(0,1,2), B'=(3,4,5), A on {u,u',x=6}, tunnel stub x'=7 (hard edge 6-7).
Wanted: with x isolated the only zero-excess state pairs (s,s') of (B,B') form a relation without odd closed walks
(e.g. NEQ), a single extra hard edit cannot relax it unless x is non-isolated, and x non-isolated gives some EQ pair at cost +1."""
import itertools, numpy as np, sys, json
from w6c_lib import *
n = 8
lab, same, P = part_table(n)
same = same.astype(np.float32)
pidx = {p: i for i, p in enumerate(P)}
SIG = {(0,1,0):'L', (1,0,0):'R', (1,1,1):'M'}
def sig(cols):
    S = same[:, cols].astype(int)
    out = np.full(len(S), -1)
    for k, (key_, nm) in enumerate(SIG.items()):
        m = (S == np.array(key_)).all(1); out[m] = k
    return out
sB = sig([pidx[(0,1)], pidx[(1,2)], pidx[(0,2)]])
sBp = sig([pidx[(3,4)], pidx[(4,5)], pidx[(3,5)]])
names = ['L','R','M']
def run(u, up, ashape):
    soft = {(0,1),(1,2),(0,2),(3,4),(4,5),(3,5)}
    tri = [u, up, 6]
    # A shape: centre index
    c = tri[ashape]; ends = [t for t in tri if t != c]
    a_edges = {key(c, ends[0]), key(c, ends[1])}
    a_pairs = {key(*p) for p in itertools.combinations(tri, 2)}
    soft |= a_pairs
    fixed_edges = {(0,1),(1,2),(3,4),(4,5),(6,7)} | a_edges
    hard = [p for p in P if p not in soft and p != (6,7) and 7 not in p]
    assert len(hard) == 12, len(hard)
    Hm = np.array([[ (h >> i) & 1 for i in range(12)] for h in range(4096)], np.float32)
    E = np.zeros((4096, 28), np.float32)
    for p in fixed_edges: E[:, pidx[p]] = 1
    for i, p in enumerate(hard): E[:, pidx[p]] = Hm[:, i]
    cost = same @ (1 - 2 * E.T) + E.sum(1)[None, :]        # (4140, 4096)
    iso = np.ones(len(same), bool)
    for v in range(6):
        iso &= ~(same[:, pidx[key(v, 6)]] > 0)
    g = np.full((2, 3, 3, 4096), 99, np.float32)
    for i in (0, 1):
        for s in range(3):
            for t in range(3):
                m = (iso == bool(i)) & (sB == s) & (sBp == t)
                if m.any(): g[i, s, t] = cost[m].min(0)
    return g, hard, E
res = []
tot = 0
for u in (0, 1):
    for up in (3, 4):
        for ash in range(3):
            g, hard, E = run(u, up, ash)
            for h in range(4096):
                G1 = g[1, :, :, h]; G0 = g[0, :, :, h]
                free = (G1 == 3)
                if not free.any(): continue
                F = free.astype(int)
                if np.trace(np.linalg.matrix_power(F, 3)) != 0: continue
                if np.trace(np.linalg.matrix_power(F, 5)) != 0: continue
                if (G1[~free] < 5).any(): continue
                if not ((G0 == 4) & ~free).any(): continue
                res.append((u, up, ash, h, F.tolist(), (G0 == 4).astype(int).tolist()))
print(len(res))
import pickle; pickle.dump(res, open('/tmp/claude-0/scratch/coupler_res.pkl', 'wb'))
for r in res[:20]: print(r)
