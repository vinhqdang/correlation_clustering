"""Exp 1/2/4: l = OPT-|H|; Conjecture C1 search; optimal pointer endpoints. Seeds fixed. usage: w2d_exp1.py NMIN NMAX PER_CELL"""
import sys, json, random, itertools, collections
from multiprocessing import Pool
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from w2d_common import *

def orbits_of(n, E, H, verts):
    G = nx.Graph()
    for v in range(n): G.add_node(('v', v), t='v')
    for q in E: G.add_edge(('v', q[0]), ('v', q[1]), t='g')
    for i, P in enumerate(H):
        G.add_node(('p', i), t='p')
        for r, v in zip('abc', P): G.add_edge(('p', i), ('v', v), t=r)
    def marked(v):
        X = G.copy(); X.nodes[('v', v)]['t'] = 'm'; return X
    reps = []
    for v in verts:
        Xv = marked(v); hit = False
        for r in reps:
            if nx.is_isomorphic(Xv, marked(r), node_match=lambda a, b: a['t'] == b['t'], edge_match=lambda a, b: a['t'] == b['t']): hit = True; break
        if not hit: reps.append(v)
    return len(reps)

def analyse(args):
    fam, n, seed = args
    rng = random.Random(f"{fam}-{n}-{seed}")
    E = FAMILIES[fam](rng, n)
    O = opt(n, E)
    ps = list(itertools.combinations(range(n), 2)); idx = {p: i for i, p in enumerate(ps)}
    mp, Hmax = max_packing_size(n, E)
    out = []
    packs = [("rand1", rand_maximal_packing(rng, n, E)), ("rand2", rand_maximal_packing(rng, n, E)), ("max", Hmax)]
    for tag, H in packs:
        l = O - len(H); V = vh(H); R = [v for v in range(n) if v not in V]; PP = packed_pairs(H)
        c = max([sum(v in P for P in H) for v in range(n)] or [0])
        rec = dict(fam=fam, n=n, seed=seed, tag=tag, opt=O, H=len(H), l=l, R=len(R), c=c, maxpack=mp)
        if l > 0 and len(R) >= 1:
            RR = np.array([1.0 if (p not in PP and p[0] not in V and p[1] not in V) else 0 for p in ps])
            RT = np.array([1.0 if (p not in PP and (p[0] not in V or p[1] not in V)) else 0 for p in ps])
            rec["minRR"] = int(sum(RR[idx[p]] for p in solve(n, E, w=RR, cap=O)[1]))
            rec["minRT"] = int(round(sum(RT[idx[p]] for p in solve(n, E, w=RT, cap=O)[1])))
            # endpoints of unpacked edits in some optimum
            ends = []
            for v in range(n):
                row = np.array([1.0 if (p not in PP and v in p) else 0 for p in ps])
                if row.sum() == 0: continue
                if solve(n, E, cap=O, extra=[(row, 1, 100)]) is not None: ends.append(v)
            rec["ends"] = ends
            rec["end_orbits"] = orbits_of(n, E, H, ends) if ends else 0
            rec["ends_in_R"] = sum(v in R for v in ends)
            if rec["minRR"] > 0:
                json.dump(dict(n=n, E=sorted(E), H=H, opt=O, l=l, note="C1 counterexample: every optimum has an unpacked R-R edit"),
                          open(f"instances_w2d/c1_counter_{fam}_{n}_{seed}_{tag}.json", "w"))
            if rec["minRT"] > 0 and n <= 10:
                json.dump(dict(n=n, E=sorted(E), H=H, opt=O, l=l, note="every optimum has an unpacked edit touching R"),
                          open(f"instances_w2d/forcedRT_{fam}_{n}_{seed}_{tag}.json", "w"))
        out.append(rec)
    return out

if __name__ == "__main__":
    nmin, nmax, per = map(int, sys.argv[1:4])
    jobs = [(f, n, s) for n in range(nmin, nmax + 1) for f in FAMILIES for s in range(per)]
    with Pool(4) as p: res = [r for rs in p.imap_unordered(analyse, jobs, chunksize=2) for r in rs]
    json.dump(res, open(f"instances_w2d/exp1_results_{nmin}_{nmax}_{per}.json", "w"))
    print("done", len(res))
