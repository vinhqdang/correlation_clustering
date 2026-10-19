"""w6b_conflict.py -- brute-force check of the wave-6 statements
 (W1) conflict weights: for ANY partition (any cost) and every soft pair (z,b):
        merged   : every y in Delta(z,b) has zy or by edited (both hard pairs);
        unmerged : every y in N_h(z) & N_h(b) has zy or by edited (deleted hard edge).
 (W2) defect bound: sum over soft pairs of c_state(pair) <= 2 D |U|   (U = edited hard pairs)
 (DEL) a deleted hard edge zy forces >= |N_h(z) & N_h(y)| further edited hard pairs: |U| >= 1 + |N_h(z)&N_h(y)|
 (TR) transparency: for any kept hard edge az, rows agree: x(a,t) = x(z,t) for all t (trivial; sanity)
Instances: random pair-maximal (G,H); all partitions with cost <= |H| + slack (slack=3).  Controls: the converse
direction (conflicts must be repaired by an edit of the OTHER type) is not claimed.
Usage: python3 w6b_conflict.py [n] [n_inst] [seed]
"""
import sys, random, itertools
import numpy as np
from w6b_lib import *

def check_instance(inst, slack=3):
    n = inst.n
    inst.solve()
    opt = int(inst.cost.min())
    sel = inst.P[inst.cost <= len(inst.H) + slack]
    pl = pair_list(n)
    soft = sorted(tuple(sorted(p)) for p in inst.soft)
    # precompute conflict candidate y lists
    merged_list = {}; unmerged_list = {}
    for (z, b) in soft:
        merged_list[(z, b)] = sorted(delta(inst, z, b))
        unmerged_list[(z, b)] = sorted(common_hard(inst, z, b))
    stats = dict(parts=0, w1=0, w2=0, dele=0, tr=0, maxratio=0.0, conflicts=0, edits=0, tight2=0)
    for lab in sel:
        stats['parts'] += 1
        U, S = inst.edits(lab)
        Uset = set(U)
        def ed(a, b): return ((a, b) if a < b else (b, a)) in Uset
        tot = 0
        for (z, b) in soft:
            if lab[z] == lab[b]:
                for y in merged_list[(z, b)]:
                    tot += 1
                    if not (ed(z, y) or ed(b, y)): stats['w1'] += 1
            else:
                for y in unmerged_list[(z, b)]:
                    tot += 1
                    if not (ed(z, y) or ed(b, y)): stats['w1'] += 1
        stats['conflicts'] += tot; stats['edits'] += len(U)
        if tot > 2 * inst.D * len(U): stats['w2'] += 1
        if len(U) > 0:
            stats['maxratio'] = max(stats['maxratio'], tot / (2 * inst.D * len(U)))
            if tot == 2 * inst.D * len(U): stats['tight2'] += 1
        for (a, b) in U:
            if inst.adj[a][b] and lab[a] != lab[b]:  # deleted hard edge
                if len(U) < 1 + len(inst.Nh[a] & inst.Nh[b]): stats['dele'] += 1
        # transparency
        for a in range(n):
            for z in inst.Nh[a]:
                if lab[a] == lab[z]:
                    for t in range(n):
                        if t not in (a, z) and (lab[a] == lab[t]) != (lab[z] == lab[t]): stats['tr'] += 1
    return stats

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    ninst = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = random.Random(seed)
    tot = dict(inst=0, parts=0, w1=0, w2=0, dele=0, tr=0, conflicts=0, edits=0, tight2=0); mr = 0.0
    for _ in range(ninst):
        adj = rand_graph(n, rng)
        H = random_pm_packing(adj, rng, maxc=rng.choice([1, 2, 3, None]), drop=rng.choice([0, 0.3]))
        if not H: continue
        inst = Inst(adj, H)
        assert is_packing(adj, H) and pair_maximal(adj, H)
        st = check_instance(inst)
        tot['inst'] += 1
        for k in ('parts', 'w1', 'w2', 'dele', 'tr', 'conflicts', 'edits', 'tight2'): tot[k] += st[k]
        mr = max(mr, st['maxratio'])
    print('n=%d seed=%d' % (n, seed), tot, 'max sum/(2D|U|)=%.3f' % mr)

if __name__ == '__main__':
    main()
