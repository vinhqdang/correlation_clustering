"""Audit 7: Proposition FAR (path + chord) and the NF-bad refutation examples, from scratch."""
import itertools, sys
import numpy as np, networkx as nx
import a7_lib as L

def inter_graph(H):
    I = nx.Graph(); I.add_nodes_from(range(len(H)))
    for i, j in itertools.combinations(range(len(H)), 2):
        if set(H[i]) & set(H[j]): I.add_edge(i, j)
    return I

def locally_bad(A, H):
    n = len(A); S = L.soft_of(n, H); Dl, Gm = L.delta_gamma(A, S); bad = []
    for t in H:
        x, y, z = t
        sp = lambda a, b: tuple(sorted((a, b)))
        allsoft = [sp(x, y), sp(y, z), sp(x, z)]
        states = [{sp(y, z)}, {sp(x, y)}, set(allsoft)]     # delete xy / delete yz / insert xz  (merged soft pairs)
        b = True
        for s in states:
            if sum(len(Dl[p]) if p in s else len(Gm[p]) for p in allsoft) == 0: b = False
        bad.append(b)
    return bad

def optima(A):
    n = len(A); P = L.partitions(n); cs, X, pl = L.costs(A, P)
    o = cs.min(); return int(o), P[cs == o], X[cs == o], pl

def edits_info(A, H, row, pl, K):
    """for one optimal partition (X row): list of (pair, deleted?, two_sided, far_dist, ends_bad)"""
    n = len(A); S = L.soft_of(n, H); idx = {p: k for k, p in enumerate(pl)}
    I = inter_graph(H); dist = dict(nx.all_pairs_shortest_path_length(I))
    mem = [[i for i, t in enumerate(H) if v in t] for v in range(n)]
    bad = locally_bad(A, H)
    Q = [set(y for y in range(n) if S[v, y]) for v in range(n)]
    out = []
    def stat(p, drop=None):
        if p == drop: return bool(A[p[0]][p[1]])
        return bool(row[idx[p]])
    for k, p in enumerate(pl):
        if S[p]: continue
        if bool(row[k]) == bool(A[p]): continue
        a, b = p
        def witnesses(side):
            res = []
            for t in range(n):
                if t in (a, b): continue
                sp = lambda x, y: tuple(sorted((x, y)))
                e = [stat(sp(a, b), p), stat(sp(a, t), p), stat(sp(b, t), p)]
                if sum(e) == 2 and True:
                    # induced P3 on {a,b,t}: exactly two edges
                    res.append(t)
            return res
        W = witnesses(0)
        ta = [t for t in W if t in Q[a]]; tb = [t for t in W if t in Q[b]]
        two = bool(ta) and bool(tb)
        d = min((dist[i].get(j, 10 ** 6)) for i in mem[a] for j in mem[b])
        endsbad = any(bad[i] for i in mem[a] + mem[b])
        out.append((p, bool(A[p]), two, d, endsbad))
    return out

if __name__ == "__main__":
    # ---------------- FAR
    for m in (3, 4, 5):
        n = 2 * m + 1; A = np.zeros((n, n), int)
        for i in range(n - 1): A[i][i + 1] = A[i + 1][i] = 1
        A[1][2 * m - 1] = A[2 * m - 1][1] = 1
        H = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
        assert L.is_packing(A, H) and L.is_pair_maximal(A, H)
        I = inter_graph(H); assert nx.is_isomorphic(I, nx.path_graph(m))
        o, P, X, pl = optima(A)
        sep = all(P[r][1] != P[r][2 * m - 1] for r in range(len(P)))
        bad = locally_bad(A, H)
        info = [edits_info(A, H, X[r], pl, 3) for r in range(len(P))]
        far = [[e for e in inf if e[3] >= m - 1] for inf in info]
        print("m", m, "n", n, "OPT", o, "|H|+1", len(H) + 1, "#opt", len(P), "all separate chord ends:", sep, "bad members", [i for i, b in enumerate(bad) if b],
              "every optimum has a far (>=m-1) hard edit:", all(len(f) > 0 for f in far), "two-sided:", [any(e[2] for e in f) for f in far])
    # larger m by constrained subset DP: OPT and best cost with 1 and 2m-1 in the same cluster
    def dp_with_constraint(A, a, b):
        n = len(A); N = 1 << n
        adj = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
        inner = [0] * N
        for mm in range(1, N):
            low = (mm & -mm).bit_length() - 1; r = mm & (mm - 1)
            inner[mm] = inner[r] + bin(adj[low] & r).count("1")
        f = [bin(mm).count("1") * (bin(mm).count("1") - 1) // 2 - 2 * inner[mm] for mm in range(N)]
        INF = 10 ** 9; dp = [INF] * N; dp[0] = 0
        for mm in range(1, N):
            low = mm & -mm; rest = mm ^ low; sub = rest; best = INF
            while True:
                c = sub | low; v = dp[mm ^ c] + f[c]
                if v < best: best = v
                if sub == 0: break
                sub = (sub - 1) & rest
            dp[mm] = best
        tot_e = sum(bin(x).count("1") for x in adj) // 2
        # constrained: cluster containing a also contains b
        best = INF; full = N - 1; ab = (1 << a) | (1 << b)
        for mm in range(N):
            if mm & ab == ab: best = min(best, f[mm] + dp[full ^ mm])
        return tot_e + dp[full], tot_e + best
    for m in (6, 7):
        n = 2 * m + 1; A = np.zeros((n, n), int)
        for i in range(n - 1): A[i][i + 1] = A[i + 1][i] = 1
        A[1][2 * m - 1] = A[2 * m - 1][1] = 1
        o, oc = dp_with_constraint(A, 1, 2 * m - 1)
        print("m", m, "n", n, "OPT", o, "(m+1 =", m + 1, ") best with chord kept (1,2m-1 in same cluster):", oc)
    # ---------------- NF-bad examples
    ex = [
        ("dist2 n=10", [(0,1),(0,6),(1,2),(2,3),(2,9),(3,4),(4,5),(5,6),(6,7),(7,8),(7,9)], [(0,1,2),(2,3,4),(4,5,6),(6,7,8),(2,9,7)], 2),
        ("dist3 n=10", [(0,1),(0,8),(1,2),(1,4),(2,3),(3,4),(4,5),(4,9),(5,6),(6,7),(7,8)], [(0,1,2),(2,3,4),(4,5,6),(6,7,8),(1,4,9)], 3),
    ]
    for name, E, H, K in ex:
        n = 10; A = np.zeros((n, n), int)
        for i, j in E: A[i][j] = A[j][i] = 1
        # centre check: H triples are given as (end, centre?, end?) in the note: find the induced-P3 ordering
        Hn = []
        for t in H:
            for perm in itertools.permutations(t):
                if A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and not A[perm[0]][perm[2]]: Hn.append(perm); break
            else: raise SystemExit("not P3 " + str(t))
        assert len(Hn) == len(H)
        print(name, "packing", L.is_packing(A, Hn), "pair-maximal", L.is_pair_maximal(A, Hn), "I(H) edges", sorted(inter_graph(Hn).edges()), "c =", max(sum(1 for t in Hn if v in t) for v in range(n)))
        o, P, X, pl = optima(A)
        bad = locally_bad(A, Hn)
        print("  OPT", o, "|H|", len(Hn), "#optima", len(P), "bad members", [i for i, b in enumerate(bad) if b])
        for r in range(len(P)):
            inf = edits_info(A, Hn, X[r], pl, K)
            print("  optimum", P[r].tolist(), [(e[0], "del" if e[1] else "ins", "two-sided" if e[2] else "one-sided", "dist", e[3], "touches-bad", e[4]) for e in inf])
