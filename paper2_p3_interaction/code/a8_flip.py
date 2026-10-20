"""Audit 8: independent re-verification of (i) the 7-vertex flip gadget G* and its mini-track, (ii) the cost-sharing example of wave6 Sec 4.1."""
import itertools, numpy as np
from a8_lib import *
def mk(n, edges):
    A = np.zeros((n, n), np.int8)
    for u, v in edges: A[u, v] = A[v, u] = 1
    return A
# G*: vertex 2 unused (isolated)
E = [(0,1),(0,4),(0,5),(1,4),(1,5),(4,5),(0,3),(1,3),(3,6),(3,7),(0,6),(5,6)]
Hs = [(6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1)]
I = Inst(mk(8, E), Hs)
print('G*: valid', I.valid(), 'pair-maximal', I.is_pair_maximal(), 'OPT', I.opt(), '|H|', len(Hs))
def blocks(I, k):
    d = {}
    for v in range(I.n): d.setdefault(int(I.lab[k, v]), []).append(v)
    return sorted(b for b in d.values() if len(b) > 1 or True)
for k in I.zero_excess():
    bl = [b for b in blocks(I, k) if b != [2]]
    x_single = all(I.lab[k, v] != I.lab[k, 7] for v in range(8) if v != 7)
    print(bl, 'x singleton' if x_single else 'x not singleton', 'states m0..m4 =', I.state(k))
# mini-track: t=8 hard-adjacent to {0,1,4,5,6}; stub 9 on x=7
E2 = E + [(8,v) for v in (0,1,4,5,6)]
I2 = Inst(mk(9, E2), Hs); print('mini-track no stub: valid', I2.valid(), 'pairmax', I2.is_pair_maximal(), 'OPT', I2.opt(), 'ZE#', len(I2.zero_excess()),
      [ [b for b in blocks(I2,k) if b!=[2]] for k in I2.zero_excess()])
I3 = Inst(mk(10, E2 + [(7, 9)]), Hs); print('mini-track + stub on x: valid', I3.valid(), 'OPT', I3.opt())
# without t but with stub (tunnel respected requires x singleton)
I4 = Inst(mk(9, E + [(7, 8)]), Hs); print('G* + stub on x: OPT', I4.opt())
# --- cost sharing (wave6 4.1): P=(0,1,2) P'=(3,4,5) z=2 has hard neighbour t=6, z'=5 has hard neighbour t'=7
base = [(0,1),(1,2),(3,4),(4,5),(2,6),(5,7)]
Hc = [(0,1,2),(3,4,5)]
def opt_of(n, edges): J = Inst(mk(n, edges), Hc); return J.valid(), J.is_pair_maximal(), J.opt()
print('sharing: no tunnel', opt_of(8, base))
print('sharing: private stubs at x=0 and x\'=3', opt_of(10, base + [(0, 8), (3, 9)]))
print('sharing: one tunnel 0-3', opt_of(8, base + [(0, 3)]))
# same tunnel but states L allowed: remove hard neighbours of z, z'
base2 = [(0,1),(1,2),(3,4),(4,5)]
print('sharing, default state L allowed, one tunnel 0-3', opt_of(6, base2 + [(0, 3)]))
# Sigma(G*) is not 2-decomposable: RRRMM is pairwise consistent but not a solution
vecs = sorted({I.state(k) for k in I.zero_excess()}); print('Sigma(G*) =', vecs)
m = 5; proj = {(i, j): {(v[i], v[j]) for v in vecs} for i in range(m) for j in range(i + 1, m)}
t = 'RRRMM'; print('RRRMM pairwise consistent:', all((t[i], t[j]) in proj[(i, j)] for (i, j) in proj), 'in Sigma:', t in vecs)
print('c(G*) =', max(sum(1 for tt in Hs if v in tt) for v in range(8)))
