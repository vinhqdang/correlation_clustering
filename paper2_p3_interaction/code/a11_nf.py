"""Audit 11: naive normal form (all hard edits have an endpoint outside V(H*)) and NF* (A-local or uncovered endpoint), exhaustive over
connected graphs, maximum packings x all optima, ground truth by full partition enumeration.
usage: python3 a11_nf.py n [res mod]   (res/mod: shard of the geng list)   or  python3 a11_nf.py n sample SEED COUNT"""
import sys, random
from a11_lib import *
import networkx as nx, subprocess

def analyse(g):
    nu, MP = max_packings(g)
    o, opts = opt_and_optima(g)
    naive = False; nfs = False
    for M in MP:
        own = soft_owner(M); cov = {v for t in M for v in t}
        W = {v: set() for v in range(g.n)}
        for (a, b) in own: W[a].add(b); W[b].add(a)
        for S in opts:
            U = [p for p in S if p not in own]
            nv = all(a not in cov or b not in cov for a, b in U)
            nf = all(a not in cov or b not in cov or (W[a] & W[b]) for a, b in U)
            naive |= nv; nfs |= nf
            if naive and nfs: return nu, o, True, True
    return nu, o, naive, nfs

def graphs(n, res=None, mod=None):
    args = ['nauty-geng', '-q', '-c', str(n)]
    if mod: args.append(f'{res}/{mod}')
    out = subprocess.run(args, capture_output=True, text=True).stdout.split()
    for s in out:
        yield nx.from_graph6_bytes(s.encode())

if __name__ == '__main__':
    n = int(sys.argv[1])
    if len(sys.argv) > 2 and sys.argv[2] == 'sample':
        rng = random.Random(int(sys.argv[3])); cnt = int(sys.argv[4])
        gs = []
        while len(gs) < cnt:
            p = rng.choice([0.25, 0.35, 0.5, 0.65])
            G = nx.gnp_random_graph(n, p, seed=rng.randint(0, 10**9))
            if nx.is_connected(G): gs.append(G)
    else:
        res, mod = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (None, None)
        gs = graphs(n, res, mod)
    tot = withp3 = gt = naivefail = nffail = 0
    ex = None; exnf = None
    for G in gs:
        g = Gr(n, G.edges()); tot += 1
        if not g.p3(): continue
        withp3 += 1
        nu, o, naive, nfs = analyse(g)
        if o > nu: gt += 1
        if not naive:
            naivefail += 1
            if ex is None or len(g.E) < len(ex.E): ex = g
        if not nfs:
            nffail += 1; exnf = g
    print(f'n={n} {sys.argv[2:]}: connected={tot} withP3={withp3} OPT>nu={gt} naive-fail={naivefail} NF*-fail={nffail}', 'minimal naive example edges:', sorted(ex.E) if ex else None, 'NF* counterexample:', sorted(exnf.E) if exnf else None)
