"""Verify the closing/joinability lemma (Lemma LOCK) on random tight instances.
Setting: member Q=(x,y,z) (centre y), x clean (all pairs at x except those in Q are hard non-edges, x in exactly one member).
Claim A: if zeta is ZE and x is not a singleton then Q in {R,M}; Q=M -> moving x out gives a ZE solution with x singleton;
         Q=R -> y has no hard edge, y's cluster is {x,y}, y is an end of every other member through it, those members are in state L.
Claim B: in the Q=R case let C = cluster of z in zeta.  Moving x to a singleton and y into C gives a ZE solution iff
         every vertex of C\\{z} is a soft partner of y (pairs (y,c) all soft).
"""
import sys, os, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *

def clusters_of(lab):
    d = {}
    for v, l in enumerate(lab): d.setdefault(l, []).append(v)
    return d

def is_ze(n, E, members, lab):
    soft = set()
    for m in members: soft |= pairs_of(m)
    for p in itertools.combinations(range(n), 2):
        if p not in soft and ((p in E) != (lab[p[0]] == lab[p[1]])): return False
    for (a, b, c) in members:
        cnt = (lab[a] != lab[b]) + (lab[b] != lab[c]) + (lab[a] == lab[c])
        if cnt != 1: return False
    return True

def random_instance(rng, n, M):
    # random graph with random members: pick members greedily, set edges
    for _ in range(200):
        E = set(); members = []; used = set()
        verts = list(range(n))
        for _ in range(M * 10):
            if len(members) >= M: break
            a, b, c = rng.sample(verts, 3)
            ps = pairs_of((a, b, c))
            if ps & used: continue
            members.append((a, b, c)); used |= ps
        if len(members) < M: continue
        soft = used
        # edges: member edges, plus random hard edges
        for (a, b, c) in members: E |= {key(a, b), key(b, c)}
        for p in itertools.combinations(range(n), 2):
            if p not in soft and rng.random() < 0.25: E.add(p)
        ok, _ = check_instance(n, E, members)
        if ok: return E, members
    return None

def main(trials, seed):
    rng = random.Random(seed)
    stats = dict(inst=0, sol=0, caseR=0, caseM=0, condTrue=0, condFalse=0, surgery_ok=0, surgery_fail=0, claimA_fail=0, claimB_mismatch=0)
    for t in range(trials):
        n = rng.randint(7, 10); M = rng.randint(3, 5)
        r = random_instance(rng, n, M)
        if r is None: continue
        E, members = r
        # choose a member and an end x that is clean: x in no other member, no hard edge
        soft = set()
        for m in members: soft |= pairs_of(m)
        for qi, (a, b, c) in enumerate(members):
            for (x, y, z) in ((a, b, c), (c, b, a)):
                cnt = sum(1 for m in members if x in m)
                if cnt != 1: continue
                if any((key(x, v) in E and key(x, v) not in soft) for v in range(n) if v != x): continue
                stats['inst'] += 1
                sols = ze_solve(n, E, members, all_solutions=True, limit=3000)
                for lab in sols:
                    stats['sol'] += 1
                    cl = clusters_of(lab)
                    if len(cl[lab[x]]) == 1: continue   # already singleton
                    tab, tbc, tac = lab[x] == lab[y], lab[y] == lab[z], lab[x] == lab[z]
                    if (tab, tbc, tac) == (True, True, True):
                        stats['caseM'] += 1
                        lab2 = list(lab); lab2[x] = max(lab) + 1
                        if not is_ze(n, E, members, lab2): stats['claimA_fail'] += 1
                        continue
                    if (tab, tbc, tac) != (True, False, False):
                        stats['claimA_fail'] += 1; continue
                    stats['caseR'] += 1
                    if cl[lab[x]] != sorted([x, y]) and sorted(cl[lab[x]]) != sorted([x, y]): stats['claimA_fail'] += 1; continue
                    C = [v for v in range(n) if lab[v] == lab[z]]
                    cond = all(key(y, cc) in soft for cc in C if cc != z)
                    lab2 = list(lab); lab2[x] = max(lab) + 1; lab2[y] = lab[z]
                    ok = is_ze(n, E, members, lab2)
                    stats['condTrue' if cond else 'condFalse'] += 1
                    if cond and ok: stats['surgery_ok'] += 1
                    elif cond and not ok: stats['surgery_fail'] += 1
                    elif (not cond) and ok: stats['claimB_mismatch'] += 1
    print(stats)

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
