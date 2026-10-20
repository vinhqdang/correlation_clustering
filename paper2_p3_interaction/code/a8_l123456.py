"""Audit 8: Lemmas 1-6 and 'inactive anchors cost nothing' (brute force over all partitions)."""
import sys, random, itertools
import numpy as np
from a8_lib import *

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
N_INST = int(sys.argv[2]) if len(sys.argv) > 2 else 600
stats = {}
def bump(k, v=1): stats[k] = stats.get(k, 0) + v
viol = []

def lemma1(I):
    # decomposition: cost = sum c_i + |U|, c_i >= 1 for EVERY partition
    nH = len(I.H); mem_idx = [[I.pidx[p] for p in mpairs(t)] for t in I.H]
    ci = np.stack([I.diff[:, ix].sum(1) for ix in mem_idx], 1)
    U = (I.diff & I.hardmask).sum(1)
    if not (ci.min() >= 1 and (ci.sum(1) + U == I.cost).all()): viol.append(('decomp', I.A.tolist(), I.H)); return
    bump('L1_decomp_partitions', len(I.cost))
    ZE = I.zero_excess()
    comps = hard_components(I)
    for k in ZE:
        bump('L1_ZE')
        # (a) hard edges inside clusters, hard non-edges across
        for (i, j) in I.P:
            if (i, j) not in I.soft:
                if I.together(k, i, j) != bool(I.A[i, j]): viol.append(('L1a', k)); return
        # (b) member states in {L,R,M}
        st = I.state(k)
        if '?' in st: viol.append(('L1b', k)); return
        # (c) every cluster union of whole hard components; cross pairs soft
        for comp in comps:
            if len({I.lab[k, v] for v in comp}) != 1: viol.append(('L1c1', k)); return
        for u, v in itertools.combinations(range(I.n), 2):
            if I.lab[k, u] == I.lab[k, v] and (u, v) not in I.soft:
                if not I.A[u, v]: viol.append(('L1c2', k)); return

def hard_components(I):
    n = I.n; par = list(range(n))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for u, v in itertools.combinations(range(n), 2):
        if I.A[u, v] and (u, v) not in I.soft: par[f(u)] = f(v)
    g = {}
    for v in range(n): g.setdefault(f(v), []).append(v)
    return list(g.values())

def lemma56(I):
    ZE = I.zero_excess()
    for k in ZE:
        for t in I.H:
            for (y, z) in [(t[0], t[1]), (t[1], t[2]), (t[0], t[2])]:
                if I.together(k, y, z):
                    for q in range(I.n):
                        if q in (y, z): continue
                        bump('L5_checks')
                        if I.together(k, q, y) != I.together(k, q, z): viol.append(('L5', k)); return
                        # forced soft pair value when (q,z) hard
                        if key(q, z) not in I.soft and key(q, y) in I.soft:
                            if I.together(k, q, y) != bool(I.A[q, z]): viol.append(('L5forced', k)); return
        # Lemma 6
        for u, w, v in itertools.permutations(range(I.n), 3):
            if u < v and I.together(k, u, w) and I.together(k, w, v):
                bump('L6_checks')
                if not (I.A[u, v] or key(u, v) in I.soft): viol.append(('L6', k)); return

def add_vertices(I, extra_edges, nnew):
    n2 = I.n + nnew; A = np.zeros((n2, n2), np.int8); A[:I.n, :I.n] = I.A
    for (u, v) in extra_edges: A[u, v] = A[v, u] = 1
    return Inst(A, I.H)

def clean_ends(I):
    """vertices x in >=1 member, with no hard edge, as given by the definition of an anchor"""
    out = []
    for x in range(I.n):
        if I.hard_nbrs(x): continue
        if any(x in t for t in I.H): out.append(x)
    return out

def lemma2_and_inactive(I):
    ZE = I.zero_excess()
    anchors = clean_ends(I)
    for x in anchors:
        J = add_vertices(I, [(x, I.n)], 1); xp = I.n
        ZEJ = J.zero_excess()
        # Lemma 2: every ZE solution of the instance with the stub respects the tunnel and forces states
        for k in ZEJ:
            cl = [v for v in range(J.n) if J.lab[k, v] == J.lab[k, x]]
            if sorted(cl) != sorted([x, xp]): viol.append(('L2cluster', I.A.tolist(), I.H, x)); return
            for t in I.H:
                if x in t:
                    if t[1] == x: viol.append(('L2centre',)); return
                    other = [v for v in t if v != x]
                    if not J.together(k, other[0], other[1]): viol.append(('L2yz',)); return
            bump('L2_ZE')
        # inactive anchors cost nothing: OPT(G+stub)=|H|  <=>  some ZE of G has x singleton
        face = any(all(I.lab[k, v] != I.lab[k, x] for v in range(I.n) if v != x) for k in ZE)
        bump('inactive_checks'); bump('inactive_face' if face else 'inactive_nonface')
        if (J.opt() == len(I.H)) != face: viol.append(('inactive', I.A.tolist(), I.H, x)); return
        # more precisely OPT(G+stub) <= OPT(G) if face (adding isolated stub vertex then edge cut costs 0)
    # paired tunnels x1-x2: anchors in different members, Q disjoint, no hard edges, x1x2 non-edge
    for x1, x2 in itertools.combinations(anchors, 2):
        if I.A[x1, x2] or key(x1, x2) in I.soft: continue
        if I.Q(x1) & I.Q(x2): continue
        J = Inst(np.where(np.eye(I.n, dtype=bool), 0, I.A + 0), I.H)
        J.A[x1, x2] = J.A[x2, x1] = 1
        J = Inst(J.A, I.H)
        ZEJ = J.zero_excess()
        for k in ZEJ:
            cl = sorted(v for v in range(J.n) if J.lab[k, v] == J.lab[k, x1])
            if cl != sorted([x1, x2]): viol.append(('L2pair',)); return
        face = any(all(I.lab[k, v] != I.lab[k, x1] for v in range(I.n) if v != x1) and
                   all(I.lab[k, v] != I.lab[k, x2] for v in range(I.n) if v != x2) for k in I.zero_excess())
        bump('pair_checks'); bump('pair_face' if face else 'pair_nonface')
        if (J.opt() == len(I.H)) != face: viol.append(('pair_inactive', I.A.tolist(), I.H, x1, x2)); return

def lemma3_4(I):
    anchors = clean_ends(I)
    if not anchors: return
    xs = rng.sample(anchors, min(len(anchors), rng.randint(1, 2)))
    if I.n + len(xs) > 10: return
    J = add_vertices(I, [(x, I.n + i) for i, x in enumerate(xs)], len(xs))
    nH = len(J.H)
    sel = np.nonzero(J.cost <= nH + 2)[0]
    stubs = [(x, I.n + i) for i, x in enumerate(xs)]
    for k in sel:
        ex = int(J.cost[k]) - nH
        rho = sum(1 for (x, xp) in stubs if not (sorted(v for v in range(J.n) if J.lab[k, v] == J.lab[k, x]) == sorted([x, xp])))
        cuts = sum(1 for (x, xp) in stubs if not J.together(k, x, xp))
        bump('L3_checks')
        if not (2 * ex >= rho and ex >= cuts): viol.append(('L3', k)); return
        # Lemma 4 for first stub: only edited hard pair at x is the cut x-xp (and it is cut)
        x, xp = stubs[0]
        edited_hard_at_x = [p for p in J.P if x in p and p not in J.soft and J.diff[k, J.pidx[p]]]
        if edited_hard_at_x == [key(x, xp)]:
            clx = {v for v in range(J.n) if J.lab[k, v] == J.lab[k, x]}
            bump('L4_checks')
            if not clx <= ({x} | J.Q(x)): viol.append(('L4cluster', k)); return
            mem = [t for t in J.H if x in t]
            if len(mem) == 1:
                t = mem[0]; mi = J.H.index(t)
                ci = int(J.diff[k, [J.pidx[p] for p in mpairs(t)]].sum())
                # paper's sentence: "P is in a cost-1 state and x is together with y"
                if ci != 1: bump('L4_paper_sentence_ci_not_1')
                if ci == 1 and len(clx) == 1: bump('L4_paper_sentence_x_singleton_in_cost1_state')
                if ci == 1 and len(clx) > 1: bump('L4_fixed_sentence_ok')

n_ok = 0; tried = 0
while n_ok < N_INST and tried < 50 * N_INST:
    tried += 1
    n = rng.choice([6, 7, 8, 8, 9])
    I = make_inst(n, rng)
    if I is None or not I.valid(): continue
    bump('inst'); 
    if I.is_pair_maximal(): bump('pairmax')
    lemma1(I)
    if I.opt() == len(I.H):
        bump('tight'); lemma56(I); lemma2_and_inactive(I)
    elif len(I.H) <= 3: pass
    if n <= 8: lemma3_4(I)
    n_ok += 1
    if viol: break
print(stats); print('VIOLATIONS:', viol[:3])
