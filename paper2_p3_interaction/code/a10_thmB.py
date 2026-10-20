"""Theorem B (LPPH W[1]-hard via Grid Tiling with <=): independent construction + two independent solvers.
GT<=: s_ij=(a,b) in S_ij; a_ij <= a_{i,j+1}; b_ij <= b_{i+1,j}.   LPPH: paths over nodes, ops = single node or link (pair), budget L."""
import sys, random, itertools
def gt_solve(k, n, S):
    cells = [(i, j) for i in range(k) for j in range(k)]
    def rec(idx, ch):
        if idx == len(cells): return True
        i, j = cells[idx]
        for (a, b) in S[(i, j)]:
            if j > 0 and not (ch[(i, j - 1)][0] <= a): continue
            if i > 0 and not (ch[(i - 1, j)][1] <= b): continue
            ch[(i, j)] = (a, b)
            if rec(idx + 1, ch): return True
        return False
    return rec(0, {})
def build(k, n, S):
    node = {}; 
    def nid(key):
        if key not in node: node[key] = len(node)
        return node[key]
    arms = {}
    for i in range(k):
        for j in range(k):
            arms[('A', i, j)] = [nid(('A', i, j, e)) for e in sorted(S[(i, j)], key=lambda e: (e[0], e[1]))]
            arms[('B', i, j)] = [nid(('B', i, j, e)) for e in sorted(S[(i, j)], key=lambda e: (e[1], e[0]))]
    paths = [frozenset(v) for v in arms.values()]
    # rows: A arms; constraint a_ij <= a_{i,j+1}: path = {A_{i,j+1}: a>=alpha+1} U {A_ij: a<=alpha}
    for i in range(k):
        for j in range(k - 1):
            for al in range(1, n):
                p = [nid(('A', i, j + 1, e)) for e in S[(i, j + 1)] if e[0] >= al + 1] + [nid(('A', i, j, e)) for e in S[(i, j)] if e[0] <= al]
                paths.append(frozenset(p))
    for j in range(k):
        for i in range(k - 1):
            for be in range(1, n):
                p = [nid(('B', i + 1, j, e)) for e in S[(i + 1, j)] if e[1] >= be + 1] + [nid(('B', i, j, e)) for e in S[(i, j)] if e[1] <= be]
                paths.append(frozenset(p))
    links = [(nid(('A', i, j, e)), nid(('B', i, j, e))) for i in range(k) for j in range(k) for e in S[(i, j)]]
    return len(node), paths, links
def lpph_solve(N, paths, links, L):
    """hitting set with ops = single nodes or links; generic branching on an unhit path (smallest first). Independent of construction."""
    lk = {}
    for (x, y) in links: lk.setdefault(x, []).append((x, y)); lk.setdefault(y, []).append((x, y))
    def rec(unhit, budget):
        if not unhit: return True
        if budget == 0: return False
        p = min(unhit, key=len)
        if not p: return False
        for v in p:
            for op in [(v,)] + lk.get(v, []):
                nu = [q for q in unhit if not (set(op) & q)]
                if rec(nu, budget - 1): return True
        return False
    return rec(list(paths), L)
def lpph_solve_bf(N, paths, links, L):
    """pure brute force over all op subsets of size L (second, even more naive solver; tiny instances)"""
    ops = [(v,) for v in range(N)] + links
    for comb in itertools.combinations(ops, L):
        hit = set()
        for op in comb: hit.update(op)
        if all(p & hit for p in paths): return True
    return False
if __name__ == "__main__":
    seed, cnt, k, n, dens = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
    bf = len(sys.argv) > 6
    rng = random.Random(seed); y = nn = bad = 0
    for _ in range(cnt):
        S = {(i, j): [(a, b) for a in range(1, n + 1) for b in range(1, n + 1) if rng.random() < dens] for i in range(k) for j in range(k)}
        gt = gt_solve(k, n, S)
        N, paths, links = build(k, n, S)
        r = lpph_solve(N, paths, links, k * k)
        if bf: r2 = lpph_solve_bf(N, paths, links, k * k); assert r == r2
        if gt != r: bad += 1; print("MISMATCH", S)
        y += gt; nn += (not gt)
    print(f"k={k} n={n} dens={dens}: {cnt} instances, YES {y}, NO {nn}, mismatches {bad}")
