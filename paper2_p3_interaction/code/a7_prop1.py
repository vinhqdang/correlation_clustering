"""Audit 7: Prop 1 of wave6_ml: for optimal S (l=|S|-|H|): hard component K not in BIG: if some edge of K is deleted, every piece lies in an exceptional cluster
and #pieces <= l+1 ; if none deleted, K in one cluster that is exceptional (|C|<=1+2D+u) or non-exceptional (union of whole comps).
Also |K| <= Z* when S is optimal with |S|<=|H|+l."""
import sys, random, itertools
import numpy as np, networkx as nx
import a7_lib as L
rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2]); n0 = [int(x) for x in sys.argv[3].split(',')]
st = dict(inst=0, opt_parts=0, split=0, fail_piece=0, fail_count=0, fail_size=0, fail_exc=0, maxK=0, unclean=0)
for it in range(N):
    n = rng.choice(n0); A, H = L.gen_instance(n, rng)
    if not H: continue
    S = L.soft_of(n, H); P = L.partitions(n); cs, X, pl = L.costs(A, P); idx = {p: k for k, p in enumerate(pl)}
    o = cs.min(); c = max(sum(1 for t in H if v in t) for v in range(n)); D = 2 * c
    l = int(o) - len(H); B = 2 * D + 3 * l; Z = B + D + 1; Zs = max(Z, (l + 1) * (2 * D + 2 * l + 1))
    hard = A.astype(bool) & ~S; dh = hard.sum(1)
    Gh = nx.Graph(); Gh.add_nodes_from(range(n)); Gh.add_edges_from((i, j) for i, j in L.prs(n) if hard[i, j])
    comps = [set(k) for k in nx.connected_components(Gh)]
    st['inst'] += 1
    for r in np.where(cs == o)[0]:
        lab = P[r]; st['opt_parts'] += 1
        U = [p for k, p in enumerate(pl) if not S[p] and (X[r][k] != bool(A[p]))]
        u = len(U); ends = set(v for p in U for v in p)
        exc = {int(lab[v]) for v in ends}
        for K in comps:
            if any(dh[v] > B for v in K): continue
            st['maxK'] = max(st['maxK'], len(K))
            if len(K) > Zs: st['fail_size'] += 1
            kept = nx.Graph(); kept.add_nodes_from(K); kept.add_edges_from((i, j) for i, j in itertools.combinations(sorted(K), 2) if hard[i, j] and lab[i] == lab[j])
            pieces = list(nx.connected_components(kept))
            dele = [(i, j) for i, j in itertools.combinations(sorted(K), 2) if hard[i, j] and lab[i] != lab[j]]
            if dele:
                st['split'] += 1
                if len(pieces) > l + 1: st['fail_count'] += 1
                for pc in pieces:
                    if int(lab[next(iter(pc))]) not in exc: st['fail_piece'] += 1
                    if sum(1 for v in range(n) if lab[v] == lab[next(iter(pc))]) > 2 * D + l + 1: st['fail_size'] += 1
            else:
                cl = int(lab[next(iter(K))])
                members = [v for v in range(n) if lab[v] == cl]
                if cl in exc and len(members) > 1 + 2 * D + u: st['fail_exc'] += 1
print(st)
