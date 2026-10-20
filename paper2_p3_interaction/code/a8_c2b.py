"""Audit 8: Conjecture C2 -- counterexample verification and search.
(1) re-verify the n=7, I(H)=C4 instance found by a8_l7.py with fresh code (networkx for cluster-graph test);
(2) random search for c<=2 instances whose ZE set is (a) not M-majority closed, (b) not closed under any uniform majority, (c) not 2-decomposable."""
import sys, random, itertools
import numpy as np, networkx as nx
from a8_lib import *

A7 = [[0,1,1,1,0,1,0],[1,0,1,0,0,0,0],[1,1,0,1,0,1,0],[1,0,1,0,0,1,1],[0,0,0,0,0,0,0],[1,0,1,1,0,0,1],[0,0,0,1,0,1,0]]
H7 = [(2,5,6),(1,0,5),(1,2,3),(0,3,6)]
G = nx.from_numpy_array(np.array(A7))
assert all(A7[i][j] == A7[j][i] for i in range(7) for j in range(7))
def is_cluster_graph(Gr):
    return all(Gr.subgraph(c).number_of_edges() == len(c) * (len(c) - 1) // 2 for c in nx.connected_components(Gr))
def apply(Gr, edits):
    H2 = Gr.copy()
    for (u, v) in edits:
        if H2.has_edge(u, v): H2.remove_edge(u, v)
        else: H2.add_edge(u, v)
    return H2
def member_pairs(t, st):
    a, b, c = t
    return {'L': (a, b), 'R': (b, c), 'M': (a, c)}[st]
# validity of packing
for (a, b, c) in H7: assert G.has_edge(a, b) and G.has_edge(b, c) and not G.has_edge(a, c)
ps = [p for t in H7 for p in mpairs(t)]; assert len(ps) == len(set(ps))
sols = []
for st in itertools.product('LRM', repeat=4):
    ed = [member_pairs(t, s) for t, s in zip(H7, st)]
    if is_cluster_graph(apply(G, ed)): sols.append(''.join(st))
print('ZE state vectors (fresh code):', sols)
a, b, c = 'RMMR', 'LRRL', 'MLLM'
mm = ''.join(x if (x == y or x == z) else y if y == z else 'M' for x, y, z in zip(a, b, c))
print('M-majority of', a, b, c, '=', mm, 'is solution?', mm in sols)
# OPT = |H| check by all partitions
I = Inst(A7, H7); print('OPT', I.opt(), '|H|', len(H7), 'ZE partitions', len(I.zero_excess()), 'pair-maximal', I.is_pair_maximal(), 'c =', max(sum(1 for t in H7 if v in t) for v in range(7)))
print('I(H) edges', sorted((i, j) for i in comps_I(H7) for j in comps_I(H7)[i] if i < j))

# search
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1); N = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
def mmaj(a, b, c): return ''.join(x if (x == y or x == z) else y if y == z else 'M' for x, y, z in zip(a, b, c))
def two_dec(vecs):
    m = len(vecs[0]); vset = set(vecs)
    proj = {(i, j): {(v[i], v[j]) for v in vecs} for i in range(m) for j in range(i + 1, m)}
    sol = set()
    def rec(pref):
        i = len(pref)
        if i == m: sol.add(''.join(pref)); return
        for s in 'LRM':
            if all((pref[j], s) in proj[(j, i)] for j in range(i)): rec(pref + [s])
    rec([]); return sol == vset
st = {}; ex = []
done = 0; tried = 0
while done < N and tried < 300 * N:
    tried += 1
    n = rng.choice([int(x) for x in sys.argv[3].split(',')]) if len(sys.argv) > 3 else rng.choice([6, 7, 7, 8, 9])
    A = gen_planted(n, rng, nedit=rng.randint(2, 6)) if rng.random() < 0.85 else gen_random(n, rng)
    H = greedy_pack(A, rng)
    if not H: continue
    I = Inst(A, H)
    if not I.valid() or I.opt() != len(H): continue
    c = max(sum(1 for t in H if v in t) for v in range(n))
    if c > 2: continue
    vecs = sorted({I.state(k) for k in I.zero_excess()})
    if len(vecs) < 3 or len(vecs) > 50: continue
    done += 1; vs = set(vecs); st['inst'] = st.get('inst', 0) + 1
    if not all(mmaj(a, b, cc) in vs for a in vecs for b in vecs for cc in vecs): st['not_Mmajority_n%d' % n] = st.get('not_Mmajority_n%d' % n, 0) + 1; ex.append((A.tolist(), H, vecs))
    if not two_dec(vecs): st['not_2dec'] = st.get('not_2dec', 0) + 1; ex.append(('NOT2DEC', A.tolist(), H, vecs))
print(st)
for e in ex[:4]: print(e)
