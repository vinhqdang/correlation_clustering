"""w6b_classify.py -- anatomy of the edited hard pairs of OPTIMAL partitions (pair-maximal H, n <= 9).
For each optimal partition pi and each edited hard pair e=(a,b):
  kind      : 'del' (hard edge, a,b in different clusters) | 'T1' | 'T2' | 'ins?' (inserted, no common neighbour: would
              contradict Lemma half/T1T2)
  witnesses : t such that {a,b,t} is an induced P3 of G xor (S minus e).  (Every edit of an optimum has one: otherwise S - e
              would be a smaller solution.)  A witness is 'soft-anchored' if the triple contains a packed pair.
  sides     : 'a' if some witness t is a soft partner of a, 'b' likewise; two-sided = both.
  dist      : I(H)-distance between the members of a and of b (INF: different components or an end outside V(H)).
Checks (PROVED statements): (AN) every edit has a witness; if no witness is soft-anchored then another edit of U lies in a
witness triple that has no soft pair (weak anchoring lemma).
Usage: python3 w6b_classify.py n n_inst seed
"""
import sys, random, itertools
from collections import Counter
import numpy as np
from w6b_lib import *

INF = 10 ** 6

def witnesses(inst, lab, e):
    a, b = e; n = inst.n
    # final adjacency under S minus e: x'(p,q) = same(p,q), flipped on e
    def xp(p, q):
        s = lab[p] == lab[q]
        if (p, q) == (a, b) or (p, q) == (b, a): s = not s
        return s
    w = []
    for t in range(n):
        if t in (a, b): continue
        k = int(xp(a, b)) + int(xp(a, t)) + int(xp(b, t))
        if k == 2: w.append(t)
    return w

def kind_of(inst, lab, e):
    a, b = e
    if inst.adj[a][b]: return 'del'
    n = inst.n
    zs = [z for z in range(n) if z not in (a, b) and inst.adj[a][z] and inst.adj[z][b] and lab[z] == lab[a]]
    if not zs: return 'ins?'
    for z in zs:
        s1 = inst.is_soft(a, z); s2 = inst.is_soft(z, b)
        if s1 and s2: return 'T1'
    return 'T2'

def analyse(inst, parts):
    rows = []
    viol_AN = 0; viol_w = 0
    for lab in parts:
        U, S = inst.edits(lab)
        Uset = set(U)
        for e in U:
            a, b = e
            w = witnesses(inst, lab, e)
            if not w: viol_w += 1; continue
            anchored = [t for t in w if inst.is_soft(a, t) or inst.is_soft(b, t) or inst.is_soft(a, b)]
            sa = any(inst.is_soft(a, t) for t in w); sb = any(inst.is_soft(b, t) for t in w)
            # (AN): if no witness triple has a soft pair, another edit must lie in some witness triple
            if not anchored:
                other = False
                for t in w:
                    for f in ((a, t), (b, t)):
                        f = tuple(sorted(f))
                        if f in Uset: other = True
                if not other: viol_AN += 1
            d = inst.vertex_dist(a, b)
            rows.append(dict(kind=kind_of(inst, lab, e), anchored=bool(anchored), sa=sa, sb=sb, dist=d, lab=tuple(lab), e=e))
    return rows, viol_AN, viol_w

def main():
    n = int(sys.argv[1]); ninst = int(sys.argv[2]); seed = int(sys.argv[3])
    rng = random.Random(seed)
    cnt = Counter(); VAN = 0; VW = 0; ninst_done = 0
    far2 = []   # instances where EVERY optimum has a two-sided edit at I(H)-distance >= 3 or INF
    for _ in range(ninst):
        adj = rand_graph(n, rng)
        H = random_pm_packing(adj, rng, maxc=rng.choice([1, 2, 3, None]), drop=rng.choice([0, 0.3]))
        if not H: continue
        inst = Inst(adj, H); inst.solve()
        opts = inst.opt_partitions()
        ninst_done += 1
        rows, van, vw = analyse(inst, opts)
        VAN += van; VW += vw
        for r in rows:
            far = 'far' if r['dist'] >= 3 else 'near'
            sides = ('A' if r['sa'] else '') + ('B' if r['sb'] else '')
            cnt[(r['kind'], 'anch' if r['anchored'] else 'UNANCH', 'two' if (r['sa'] and r['sb']) else ('one' if (r['sa'] or r['sb']) else 'zero'), far)] += 1
        # every optimum has a far two-sided edit?
        per = {}
        for r in rows:
            per.setdefault(r['lab'], []).append(r)
        if per and all(any(x['dist'] >= 3 and x['sa'] and x['sb'] for x in v) for v in per.values()) and len(per) == len(opts):
            far2.append((adj.tolist(), H))
    print('n=%d inst=%d witness-missing=%d AN-violations=%d' % (n, ninst_done, VW, VAN))
    for k, v in sorted(cnt.items()): print(k, v)
    print('instances where every optimum has a far two-sided edit:', len(far2))
    for g, H in far2[:3]:
        edges = [(i, j) for i in range(n) for j in range(i + 1, n) if g[i][j]]
        print('  edges', edges, 'H', H)

if __name__ == '__main__':
    main()
