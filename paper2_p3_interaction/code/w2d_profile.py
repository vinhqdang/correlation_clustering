"""Exp 3: is membership of an unpacked pair in SOME optimum determined by its radius-r rooted profile?
Profile(p, r) = ball of radius r around the two endpoints in the incidence graph (G-edges + one node per packed P3,
joined to its a,b,c with labelled edges); endpoints marked. Search within one instance for two unpacked pairs p, q,
both Lemma-2 viable, with isomorphic profiles but p in some optimum and q in none.  usage: w2d_profile.py NMIN NMAX PER"""
import sys, json, random, itertools
from multiprocessing import Pool
import networkx as nx
from w2d_common import *

def inc_graph(n, E, H):
    G = nx.Graph()
    for v in range(n): G.add_node(('v', v), t='v')
    for q in E: G.add_edge(('v', q[0]), ('v', q[1]), t='g')
    for i, P in enumerate(H):
        G.add_node(('p', i), t='p')
        for r, v in zip('abc', P): G.add_edge(('p', i), ('v', v), t=r)
    return G
def profile(G, p, r):
    src = [('v', p[0]), ('v', p[1])]
    dist = {}
    for s in src:
        for k, d in nx.single_source_shortest_path_length(G, s, cutoff=r).items(): dist[k] = min(d, dist.get(k, 99))
    B = G.subgraph(dist).copy()
    for s in src: B.nodes[s]['t'] = 'm'
    return B
def same(A, B):
    return nx.is_isomorphic(A, B, node_match=lambda a, b: a['t'] == b['t'], edge_match=lambda a, b: a['t'] == b['t'])
def viable(n, E, p, PP, V, Hl):
    """Lemma 2 necessary condition for unpacked p=(u,v) with l excess: |D(p)| <= l-1 where D counts z outside Q(u)uQ(v)."""
    return True

def run(args):
    fam, n, seed = args
    rng = random.Random(f"prof-{fam}-{n}-{seed}")
    if fam == "path": E = frozenset((i, i + 1) for i in range(n - 1))
    elif fam == "cycle": E = frozenset(key(i, (i + 1) % n) for i in range(n))
    elif fam == "caterpillar": E = frozenset([(i, i + 1) for i in range(n // 2 - 1)] + [(i, n // 2 + i) for i in range(n - n // 2) if i < n // 2])
    else: E = FAMILIES[fam](rng, n)
    O = opt(n, E); res = []
    H = rand_maximal_packing(rng, n, E) if seed % 3 else max_packing_size(n, E)[1]
    l = O - len(H)
    if l <= 0: return res
    V = vh(H); PP = packed_pairs(H); ps = list(itertools.combinations(range(n), 2)); m = len(ps)
    cand = [p for p in ps if p not in PP]
    inopt = {}
    for p in cand:
        row = np.zeros(m); row[ps.index(p)] = 1
        inopt[p] = solve(n, E, cap=O, extra=[(row, 1, 1)]) is not None
    if all(inopt.values()) or not any(inopt.values()): return res
    G = inc_graph(n, E, H)
    A = [p for p in cand if inopt[p]]; B = [p for p in cand if not inopt[p]]
    best = {}
    for r in (1, 2, 3, 4):
        PA = {p: profile(G, p, r) for p in A}; PB = {q: profile(G, q, r) for q in B}
        for p in A:
            for q in B:
                if PA[p].number_of_nodes() == PB[q].number_of_nodes() and same(PA[p], PB[q]):
                    best[r] = (p, q); break
            if r in best: break
    if best:
        res.append(dict(fam=fam, n=n, seed=seed, l=l, opt=O, H=H, E=sorted(E), r_pairs={r: v for r, v in best.items()}, ndiam=nx.diameter(G) if nx.is_connected(G) else None))
    return res

if __name__ == "__main__":
    nmin, nmax, per = map(int, sys.argv[1:4])
    fams = list(FAMILIES) + ["path", "cycle", "caterpillar"]
    jobs = [(f, n, s) for n in range(nmin, nmax + 1) for f in fams for s in range(per)]
    with Pool(3) as p: out = [r for rs in p.imap_unordered(run, jobs, chunksize=2) for r in rs]
    json.dump(out, open(f"instances_w2d/profile_hits_{nmin}_{nmax}_{per}.json", "w"))
    summ = {}
    for o in out:
        for r in o["r_pairs"]:
            r = int(r); b = summ.get(r)
            if b is None or o["n"] < b["n"]: summ[r] = o
    print("hits", len(out))
    for r in sorted(summ): print("radius", r, "min n", summ[r]["n"], summ[r]["fam"], "diam", summ[r]["ndiam"], "l", summ[r]["l"])
