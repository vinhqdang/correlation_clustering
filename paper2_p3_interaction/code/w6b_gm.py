"""w6b_gm.py -- the family G_m = path P_{2m+1} (vertices 0..2m) + the hard edge {1, 2m-1};  H_m = {(2i,2i+1,2i+2): i<m}.
Claims (hand proof in notes/wave6_pointers.md, Prop. FAR):  H_m is a pair-maximal packing, I(H_m) is a path on m nodes (w=1, c=2),
OPT(G_m) = m+1 = |H_m|+1 (l = 1), and in EVERY optimum the hard edge {1,2m-1} is deleted (its ends lie in different clusters);
its ends lie in the members M_1 and M_m, at I(H)-distance m-1.  Checked here exhaustively over all set partitions for m = 3,4,5
(n = 7, 9, 11) together with the key inequality f(C) = C(|C|,2) - 2 e(C) >= -1 for every vertex subset C (equality iff |C| in {2,3}
and G[C] connected), and the 'side-additive' count min_pi [cost(pi) + #two-sided far edits] = OPT + 1.
Usage: python3 w6b_gm.py [m_max]
"""
import sys, itertools
import numpy as np
from w6b_lib import *
from w6b_classify import witnesses

def build(m):
    n = 2 * m + 1
    chain = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
    edges = [(j, j + 1) for j in range(2 * m)] + [(1, 2 * m - 1)]
    return n, mk(n, edges), chain

def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for m in range(3, mmax + 1):
        n, A, H = build(m)
        assert is_packing(A, H) and pair_maximal(A, H)
        inst = Inst(A, H)
        assert inst.c == 2 and len(inst.I[0]) == 1 and all(len(inst.I[i]) == 2 for i in range(1, m - 1))
        opt = inst.solve(); P = inst.opt_partitions()
        assert opt == m + 1, opt
        # every optimum deletes the hard edge
        for lab in P:
            assert lab[1] != lab[2 * m - 1]
            U, S = inst.edits(lab)
            assert U == [(1, 2 * m - 1)], U
        assert inst.vertex_dist(1, 2 * m - 1) == m - 1
        # subset inequality f(C) >= -1
        for r in range(1, n + 1):
            for C in itertools.combinations(range(n), r):
                e = sum(int(A[a, b]) for a, b in itertools.combinations(C, 2))
                f = r * (r - 1) // 2 - 2 * e
                assert f >= -1
                if f == -1:
                    assert r in (2, 3)
        # two-sidedness of the edit in each optimum and the side-additive count
        two = 0
        for lab in P:
            w = witnesses(inst, lab, (1, 2 * m - 1))
            sa = any(inst.is_soft(1, t) for t in w); sb = any(inst.is_soft(2 * m - 1, t) for t in w)
            two += int(sa and sb)
        # side-additive: charge a far two-sided edit once per side
        best = 10 ** 9
        Pall = inst.P; cst = inst.cost
        for idx in np.nonzero(cst <= opt + 1)[0]:
            lab = Pall[idx]; U, S = inst.edits(lab)
            extra = 0
            for e in U:
                if inst.vertex_dist(*e) >= 2:
                    w = witnesses(inst, lab, e)
                    if any(inst.is_soft(e[0], t) for t in w) and any(inst.is_soft(e[1], t) for t in w): extra += 1
            best = min(best, int(cst[idx]) + extra)
        print('m=%d n=%d: OPT=%d=|H|+1, #optima=%d (all delete {1,%d}), I(H)-distance of the edit=%d, two-sided in %d/%d optima, '
              'min side-additive cost=%d (OPT+1=%d)' % (m, n, opt, len(P), 2 * m - 1, m - 1, two, len(P), best, opt + 1))

if __name__ == '__main__':
    main()
