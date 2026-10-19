"""w6b_nf.py -- test of the CONJECTURE (NF-bad): for every pair-maximal instance there is an optimal partition in which every far
(I(H)-distance >= K, K=2 or 3; infinite counts as far) two-sided edit e=(a,b) is incident to a locally bad member
(a witness member through a, or through b, is bad).  Counter-examples to NF-bad would be instances where EVERY optimum has a far
two-sided edit with both sides at non-bad members (the case that 'all bad members in every bag' cannot cover).
Instance families: random (rand_graph), chain+random edges, odd cycles + chords, 'two chains + tunnels'.
Usage: python3 w6b_nf.py n n_inst seed K family
"""
import sys, random
from collections import Counter
import numpy as np
from w6b_lib import *
from w6b_classify import witnesses
from w6b_bad import bad_members

def far_ff_edits(inst, lab, bad, K):
    U, S = inst.edits(lab)
    out = []
    for (a, b) in U:
        if inst.vertex_dist(a, b) < K: continue
        w = witnesses(inst, lab, (a, b))
        sa = [t for t in w if inst.is_soft(a, t)]; sb = [t for t in w if inst.is_soft(b, t)]
        if not (sa and sb): continue
        Ma = {i for i in inst.member_of[a] if any(t in inst.H[i] for t in sa)}
        Mb = {i for i in inst.member_of[b] if any(t in inst.H[i] for t in sb)}
        if not (Ma & bad) and not (Mb & bad): out.append((a, b))
    return out

def gen_chain(n, rng):
    m = (n - 1) // 2
    A = np.zeros((n, n), dtype=np.int8)
    chain = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
    for (x, y, z) in chain:
        A[x, y] = A[y, x] = A[y, z] = A[z, y] = 1
    pe = rng.choice([0.05, 0.1, 0.15, 0.2])
    for a, b in pair_list(n):
        if A[a, b] or any(a in t and b in t for t in chain): continue
        if rng.random() < pe: A[a, b] = A[b, a] = 1
    return A, chain

def extend(adj, H, rng):
    H = list(H); used = packed_set(H)
    ts = induced_p3s(adj); rng.shuffle(ts)
    for t in ts:
        p = pairs3(t)
        if not p & used: H.append(t); used |= p
    return H

def main():
    n, ninst, seed, K = [int(x) for x in sys.argv[1:5]]
    fam = sys.argv[5] if len(sys.argv) > 5 else 'mixed'
    rng = random.Random(seed)
    stat = Counter(); bad_examples = []
    for _ in range(ninst):
        f = fam if fam != 'mixed' else rng.choice(['rand', 'chain'])
        if f == 'chain':
            adj, ch = gen_chain(n, rng)
            if not is_packing(adj, ch): continue
            H = extend(adj, ch, rng)
        else:
            adj = rand_graph(n, rng)
            H = random_pm_packing(adj, rng, maxc=rng.choice([1, 2, 3, None]), drop=rng.choice([0, 0.3]))
        if not H or not pair_maximal(adj, H) or not is_packing(adj, H): continue
        inst = Inst(adj, H); inst.solve()
        bad = set(bad_members(inst))
        opts = inst.opt_partitions()
        anyfar = False; some_clean = False
        for lab in opts:
            ff = far_ff_edits(inst, lab, bad, K)
            if ff: anyfar = True
            else: some_clean = True
        stat['instances'] += 1
        if anyfar: stat['some optimum has far (F,F) edit'] += 1
        if not some_clean:
            stat['EVERY optimum has far (F,F) edit'] += 1
            bad_examples.append((adj.tolist(), H))
    print('n=%d K=%d family=%s' % (n, K, fam), dict(stat))
    for g, H in bad_examples[:5]:
        edges = [(i, j) for i in range(n) for j in range(i + 1, n) if g[i][j]]
        print('  COUNTEREXAMPLE edges', edges, 'H', H)

if __name__ == '__main__':
    main()
