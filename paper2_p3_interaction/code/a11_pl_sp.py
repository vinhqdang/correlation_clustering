"""Audit 11: PL DP (real implementation) on Sp(triangle), original packing H, l=1 (no brute force; OPT = |H|+1 by ILP in a11_spine)."""
import time
from a11_lib import *
from a11_pl import *
from a11_spine import build
g, H, d = build([(0, 1), (1, 2), (2, 0)], 3)
N = d['N']
order = [2 * i for i in range(N + 1)] + [2 * i + 1 for i in range(N + 1)]   # members: H list is [P0,X0,P1,X1,...]: leaves P_i = idx 2i, spine X_i = idx 2i+1
Ia = I_graph(H)
bags, edges = elimination_decomposition(Ia, order)
inst = Inst(g, H, 1, bags, edges)
print('n', g.n, '|H|', len(H), 'width', max(len(b) for b in bags) - 1, 'hub(y)=', inst.hub.get(d['y']), 'Elig(b0)=', sorted(inst.elig[d['X'][0][1]]))
t = time.time()
val, cf = dp(inst)
print('DP value m* =', val, ' |H|+1 =', len(H) + 1, ' coverage failures:', cf, ' secs', round(time.time() - t, 1))
