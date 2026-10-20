"""Wave 6 hardness attempt: shared library (written fresh).
Graph = dict/array; vertices 0..n-1.  Cost of a partition = #pairs with [edge] != [same block].
Members are triples (a,b,c) with centre b (ab, bc edges, ac non-edge)."""
import itertools, os, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

def prs(n): return list(itertools.combinations(range(n), 2))
def key(a, b): return (a, b) if a < b else (b, a)

_PT = {}
def part_table(n):
    if n in _PT: return _PT[n]
    res = []
    def rec(i, lab, m):
        if i == n: res.append(tuple(lab)); return
        for c in range(m + 1):
            lab.append(c); rec(i + 1, lab, max(m, c + 1)); lab.pop()
    rec(0, [], 0)
    lab = np.array(res, dtype=np.int8)
    P = prs(n)
    same = np.stack([lab[:, i] == lab[:, j] for i, j in P], 1)
    _PT[n] = (lab, same, P); return _PT[n]

def adj_from_edges(n, edges):
    A = np.zeros((n, n), np.int8)
    for a, b in edges: A[a, b] = A[b, a] = 1
    return A

def pairs_of(t):
    a, b, c = t; return {key(a, b), key(b, c), key(a, c)}

def check_packing(A, H):
    seen = set()
    for (a, b, c) in H:
        if not (A[a][b] and A[b][c] and not A[a][c]): return False, ('notP3', (a, b, c))
        ps = pairs_of((a, b, c))
        if ps & seen: return False, ('overlap', (a, b, c))
        seen |= ps
    return True, None

def p3_list(A):
    n = len(A); out = []
    for b in range(n):
        nb = [x for x in range(n) if A[b][x]]
        for a, c in itertools.combinations(nb, 2):
            if not A[a][c]: out.append((a, b, c))
    return out

def is_pair_maximal(A, H):
    used = set()
    for t in H: used |= pairs_of(t)
    return all(pairs_of(t) & used for t in p3_list(A))

def brute_opt(A):
    n = len(A); lab, same, P = part_table(n)
    e = np.array([A[i][j] for i, j in P], bool)
    cost = (same ^ e).sum(1)
    return int(cost.min()), lab, same, cost, P

def ilp_opt(A, extra_obj=None, cost_le=None, fixed=None, time_limit=600):
    """Cluster editing ILP with triangle inequalities. x_ij=1 iff together.
    returns (cost, labels) of an optimum (cost counted w.r.t. A); extra_obj = tie-break vector over pairs."""
    n = len(A); P = prs(n); idx = {p: i for i, p in enumerate(P)}
    m = len(P)
    e = np.array([A[i][j] for i, j in P], float)
    c = 1 - 2 * e                    # cost = sum e + sum x*(1-2e)
    rows = []; lo = []; hi = []
    import scipy.sparse as sp
    r = []; cc = []; vv = []; k = 0
    for i, j, l in itertools.combinations(range(n), 3):
        a, b, d = idx[(i, j)], idx[(j, l)], idx[(i, l)]
        for (p, q, s) in ((a, b, d), (a, d, b), (b, d, a)):
            # x_p + x_q - x_s <= 1
            r += [k, k, k]; cc += [p, q, s]; vv += [1, 1, -1]; k += 1
    Acon = sp.csr_matrix((vv, (r, cc)), shape=(k, m))
    cons = [LinearConstraint(Acon, -np.inf, 1)]
    lb = np.zeros(m); ub = np.ones(m)
    if fixed:
        for (a, b), v in fixed.items():
            lb[idx[key(a, b)]] = v; ub[idx[key(a, b)]] = v
    obj = c.copy()
    base = e.sum()
    if cost_le is not None:
        cons.append(LinearConstraint(sp.csr_matrix(c.reshape(1, -1)), -np.inf, cost_le - base))
        obj = extra_obj if extra_obj is not None else np.zeros(m)
    elif extra_obj is not None:
        obj = c * 1000 + extra_obj
    res = milp(obj, constraints=cons, integrality=np.ones(m), bounds=Bounds(lb, ub), options={'time_limit': time_limit})
    if res.x is None: return None
    x = np.round(res.x).astype(int)
    lab = [-1] * n; nxt = 0
    for v in range(n):
        if lab[v] == -1:
            lab[v] = nxt
            for w in range(v + 1, n):
                if x[idx[(v, w)]] == 1: lab[w] = nxt
            nxt += 1
    cost = int(base + c @ x)
    return cost, lab

def cost_of(A, lab):
    n = len(A); return sum(1 for i, j in prs(n) if bool(A[i][j]) != (lab[i] == lab[j]))

def tw_exact(adjm):
    """exact treewidth by subset DP (<= ~16 nodes); adjm = list of bitmasks"""
    n = len(adjm)
    if n == 0: return -1
    arr = [0] * (1 << n)
    def Q(S, v):
        reach = 1 << v; frontier = reach; out = 0
        while frontier:
            nxt = 0; f = frontier
            while f:
                b = f & -f; f ^= b; nxt |= adjm[b.bit_length() - 1]
            nxt &= ~reach; out |= nxt & ~S; nxt &= S; reach |= nxt; frontier = nxt
        return bin(out).count('1')
    for S in range(1, 1 << n):
        best = 10 ** 9; f = S
        while f:
            b = f & -f; f ^= b; v = b.bit_length() - 1; T = S ^ b
            val = max(arr[T], Q(T, v))
            if val < best: best = val
        arr[S] = best
    return arr[(1 << n) - 1]

def tw_heur(adj_sets):
    """min-fill upper bound on treewidth for a graph given as dict node->set"""
    import networkx as nx
    from networkx.algorithms.approximation import treewidth_min_fill_in
    G = nx.Graph(); G.add_nodes_from(adj_sets)
    for u, s in adj_sets.items():
        for v in s: G.add_edge(u, v)
    if G.number_of_edges() == 0: return 0
    return treewidth_min_fill_in(G)[0]

def ih_graph(H):
    m = len(H); adj = {i: set() for i in range(m)}
    for i in range(m):
        for j in range(i + 1, m):
            if set(H[i]) & set(H[j]): adj[i].add(j); adj[j].add(i)
    return adj

def outpath(name):
    d = os.environ.get('W6C_OUT', '/tmp/w6c_out'); os.makedirs(d, exist_ok=True); return os.path.join(d, name)
