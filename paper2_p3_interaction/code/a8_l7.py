"""Audit 8: Lemma 7 (splice) and Lemma 7b (majority closure => cores <= 2; also 2-decomposable => cores <= 2)."""
import sys, random, itertools
import numpy as np, networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in
from a8_lib import *

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
N_INST = int(sys.argv[2]) if len(sys.argv) > 2 else 400
stats = {}; viol = []
def bump(k, v=1): stats[k] = stats.get(k, 0) + v

def anchors_of(I):
    # clean anchors: no hard edge at x, x in >= 1 member
    return [x for x in range(I.n) if not I.hard_nbrs(x) and any(x in t for t in I.H)]

def sigma(I):
    ZE = I.zero_excess()
    vecs = [I.state(k) for k in ZE]
    return ZE, vecs

def singleton(I, k, x): return all(I.lab[k, v] != I.lab[k, x] for v in range(I.n) if v != x)

def faces(I, ZE, anchors):
    """face set: frozenset of anchors -> bool, for subsets up to size 4"""
    sing = {x: np.array([singleton(I, k, x) for k in ZE]) for x in anchors}
    F = {}
    for r in range(0, min(4, len(anchors)) + 1):
        for S in itertools.combinations(anchors, r):
            m = np.ones(len(ZE), bool)
            for x in S: m &= sing[x]
            F[frozenset(S)] = bool(m.any())
    return F

def min_nonfaces(F):
    out = []
    for S, f in F.items():
        if f: continue
        if all(F[S - {x}] for x in S if (S - {x}) in F): out.append(S)
    return out

def decomp(I):
    adj = comps_I(I.H); G = nx.Graph(); G.add_nodes_from(adj)
    for i, s in adj.items():
        for j in s: G.add_edge(i, j)
    if G.number_of_edges() == 0: return None
    w, T = treewidth_min_fill_in(G)
    return T

def lemma7(I, ZE, vecs, anchors, F):
    T = decomp(I)
    if T is None or T.number_of_edges() == 0: return
    sets = [set(v) for v in vecs]; vset = set(vecs); m = len(I.H)
    mf = [S for S in min_nonfaces(F) if len(S) >= 2]
    mem_of = {x: frozenset(i for i, t in enumerate(I.H) if x in t) for x in anchors}
    for (u, v) in list(T.edges()):
        T2 = T.copy(); T2.remove_edge(u, v)
        side1 = nx.node_connected_component(T2, u); side2 = nx.node_connected_component(T2, v)
        mem1 = set().union(*[set(b) for b in side1]); mem2 = set().union(*[set(b) for b in side2])
        B = mem1 & mem2; M1 = mem1 - B; M2 = mem2 - B
        # members in no bag? (isolated members) -> put in M1
        rest = set(range(m)) - mem1 - mem2; M1 |= rest
        Bl = sorted(B)
        # (a) splice closure
        for s1 in vecs:
            for s2 in vecs:
                if all(s1[i] == s2[i] for i in Bl):
                    comb = ''.join(s1[i] if (i in M1 or i in B) else s2[i] for i in range(m))
                    bump('L7a_checks')
                    if comb not in vset: viol.append(('L7a', I.A.tolist(), I.H)); return
        # (b) minimal non-faces with split
        for S in mf:
            for r in range(1, len(S)):
                for S1 in itertools.combinations(sorted(S), r):
                    S2 = tuple(x for x in sorted(S) if x not in S1)
                    if all(mem_of[x] <= M1 for x in S1) and all(mem_of[x] <= M2 for x in S2):
                        # restrictions to B
                        def sig(b, Sminus):
                            ks = [k for k, vec in zip(range(len(ZE)), vecs) if all(singleton(I, ZE[k], x) for x in Sminus)]
                            return {''.join(vecs[k][i] for i in Bl) for k in ks}
                        R1 = set(); R2 = set()
                        for b in S1: R1 |= sig(b, [x for x in S if x != b])
                        for a in S2: R2 |= sig(a, [x for x in S if x != a])
                        bump('L7b_checks')
                        if R1 & R2: viol.append(('L7b', I.A.tolist(), I.H, S, S1)); return
                        bump('L7b_nontriv_Bnonempty' if Bl else 'L7b_B_empty')

# ---- majority closure
PERMS = list(itertools.permutations('LRM'))
SYM = {'L': 0, 'R': 1, 'M': 2}
def all_rules():
    # rules: 729 maps from the 6 orderings of (L,R,M) to a value
    return list(itertools.product('LRM', repeat=6))
RULES = all_rules()
RULEARR = np.array([[SYM[c] for c in r] for r in RULES], dtype=np.int64)   # (729,6)
PIDX = {p: i for i, p in enumerate(PERMS)}

def majority_rules_closed(vecs):
    m = len(vecs[0]); V = np.array([[SYM[c] for c in v] for v in vecs], dtype=np.int64)
    pw = 3 ** np.arange(m, dtype=np.int64)
    codes = set((V * pw).sum(1).tolist()); codeArr = np.array(sorted(codes))
    ok = np.ones(len(RULES), bool)
    N = len(vecs)
    for a, b, c in itertools.product(range(N), repeat=3):
        if not ok.any(): break
        x, y, z = V[a], V[b], V[c]
        base = 0; dist_pos = []
        for i in range(m):
            if x[i] == y[i] or x[i] == z[i]: base += x[i] * pw[i]
            elif y[i] == z[i]: base += y[i] * pw[i]
            else: dist_pos.append(i)
        if not dist_pos: 
            if base not in codes: ok[:] = False
            continue
        img = np.full(len(RULES), base, dtype=np.int64)
        for i in dist_pos:
            perm = ''.join('LRM'[t] for t in (x[i], y[i], z[i]))
            img += RULEARR[:, PIDX[tuple(perm)]] * pw[i]
        ok &= np.isin(img, codeArr)
    return ok

def two_decomposable(vecs, T_idx=None):
    m = len(vecs[0]); vset = set(vecs)
    proj = {(i, j): {(v[i], v[j]) for v in vecs} for i in range(m) for j in range(i + 1, m)}
    # all tuples in {L,R,M}^m with all binary projections consistent
    # backtracking
    cnt = 0; sol = set()
    def rec(pref):
        i = len(pref)
        if i == m: sol.add(''.join(pref)); return
        for s in 'LRM':
            if all((pref[j], s) in proj[(j, i)] for j in range(i)): rec(pref + [s])
    rec([])
    return sol == vset

def lemma7b(I, ZE, vecs, anchors, F):
    if len(vecs) < 3: return
    mf = min_nonfaces(F); mx = max([len(S) for S in mf], default=0)
    bump('nonface_size1', sum(1 for S in mf if len(S)==1)); bump('nonface_size2', sum(1 for S in mf if len(S)==2))
    bump('inst_7b'); 
    if mx >= 3: bump('nonface_size3_found')
    ok = majority_rules_closed(vecs)
    if ok.any():
        bump('inst_majority_closed')
        if len(anchors) >= 3: bump('inst_maj_closed_>=3anchors')
        if mx >= 3: viol.append(('L7b', I.A.tolist(), I.H)); return
        # direct face check for every subset up to size 4: face iff all pairs faces
        for S, f in F.items():
            if len(S) >= 3 and all(F[frozenset(p)] for p in itertools.combinations(S, 2)) != f: viol.append(('L7bsub', I.A.tolist(), I.H, S)); return
    else:
        bump('inst_no_majority')
    if two_decomposable(vecs):
        bump('inst_2dec')
        if mx >= 3: viol.append(('2dec', I.A.tolist(), I.H)); return
    else: bump('inst_not_2dec')
    # c <= 2 and M-majority
    c = max(sum(1 for t in I.H if v in t) for v in range(I.n))
    bump('c%d' % c)
    if not ok.any():
        bump('nomaj_c%d' % c)
        if c <= 2: print('NOMAJ', c, I.A.tolist(), I.H, vecs)

n_ok = 0; tried = 0
while n_ok < N_INST and tried < 100 * N_INST:
    tried += 1
    n = rng.choice([7, 8, 9, 10, 10])
    A = gen_planted(n, rng, nedit=rng.randint(2, 7), bigcluster=rng.random() < 0.3) if rng.random() < 0.8 else gen_random(n, rng)
    H = greedy_pack(A, rng)
    if not H: continue
    I = Inst(A, H)
    if not I.valid() or I.opt() != len(I.H): continue
    ZE, vecs = sigma(I)
    anchors = anchors_of(I)
    if len(anchors) < 2: continue
    if len(vecs) > 40: continue
    F = faces(I, ZE, anchors)
    bump('inst'); n_ok += 1
    lemma7(I, ZE, vecs, anchors, F)
    lemma7b(I, ZE, vecs, anchors, F)
    if viol: break
print(stats); print('VIOLATIONS:', viol[:2])
