"""w6b_bad.py -- locally bad members.
A member P=(x,y,z) of H has three states of member-cost 1: delete xy ({y,z},{x}), delete yz ({x,y},{z}), insert xz ({x,y,z}); a state
fixes the merged/unmerged status of its three soft pairs.  A pair (p,q) in state 'merged' has conflict weight |Delta(p,q)|, in state
'unmerged' weight |N_h(p) & N_h(q)|.  P is LOCALLY BAD if every cost-1 state has total conflict weight > 0 (static, from (G,H)).
(BAD)  PROVED claim: for every solution S with |S| = |H| + eps, #bad members <= 2 c eps.   Checked on every partition of cost <= |H|+3.
(FAR)  Experiment: in optimal partitions, a far (I(H)-distance >= 3 or infinite) two-sided edit e=(a,b): are the member(s) through
       a (resp. b) that contain a witness t bad?  (EVIDENCE only.)
Usage: python3 w6b_bad.py n n_inst seed
"""
import sys, random
from collections import Counter
import numpy as np
from w6b_lib import *
from w6b_classify import witnesses

def states(t):
    x, y, z = t
    return [dict(m={frozenset((x, y)): 0, frozenset((y, z)): 1, frozenset((x, z)): 0}),   # delete xy : {y,z}
            dict(m={frozenset((x, y)): 1, frozenset((y, z)): 0, frozenset((x, z)): 0}),   # delete yz : {x,y}
            dict(m={frozenset((x, y)): 1, frozenset((y, z)): 1, frozenset((x, z)): 1})]   # insert xz

def conflict_weight(inst, p, q, merged):
    return len(delta(inst, p, q)) if merged else len(common_hard(inst, p, q))

def bad_members(inst):
    bad = []
    for i, t in enumerate(inst.H):
        free = False
        for st in states(t):
            w = sum(conflict_weight(inst, *tuple(pr), mg) for pr, mg in st['m'].items())
            if w == 0: free = True; break
        if not free: bad.append(i)
    return bad

def main():
    n = int(sys.argv[1]); ninst = int(sys.argv[2]); seed = int(sys.argv[3])
    rng = random.Random(seed)
    viol = 0; checked = 0; slack_hist = Counter()
    farcnt = Counter()
    for _ in range(ninst):
        adj = rand_graph(n, rng)
        H = random_pm_packing(adj, rng, maxc=rng.choice([1, 2, 3, None]), drop=rng.choice([0, 0.3]))
        if not H: continue
        inst = Inst(adj, H); inst.solve()
        bad = set(bad_members(inst))
        # (BAD) on every partition with cost <= |H|+3
        for lab, cst in zip(inst.P[inst.cost <= len(H) + 3], inst.cost[inst.cost <= len(H) + 3]):
            eps = int(cst) - len(H)
            checked += 1
            if len(bad) > 2 * inst.c * eps: viol += 1
        if eps >= 0:
            pass
        l = inst.solve() - len(H)
        slack_hist[(len(bad), l)] += 0
        # (FAR)
        for lab in inst.opt_partitions():
            U, S = inst.edits(lab)
            for (a, b) in U:
                if inst.vertex_dist(a, b) < 3: continue
                w = witnesses(inst, lab, (a, b))
                sa = [t for t in w if inst.is_soft(a, t)]; sb = [t for t in w if inst.is_soft(b, t)]
                if not (sa and sb): continue
                # members through a containing a soft partner witness
                Ma = {i for i in inst.member_of[a] if any(t in inst.H[i] for t in sa)}
                Mb = {i for i in inst.member_of[b] if any(t in inst.H[i] for t in sb)}
                farcnt[('a-side has bad member', bool(Ma & bad)), ('b-side has bad member', bool(Mb & bad))] += 1
    print('n=%d: partitions checked for (BAD): %d, violations %d' % (n, checked, viol))
    for k, v in sorted(farcnt.items(), key=str): print('  far two-sided edit:', k, v)

if __name__ == '__main__':
    main()
