"""DFS version of tight_choices (all one-edit-per-member solutions S_F of G[W] that make G[W] xor S a cluster graph), with early triple pruning.
Used to test whole-component Rule 6.2 reducibility on 14-20 member components (brute force 3^m is too slow)."""
import itertools
from w6a_lib import *
from w6a_rule import members_vertices
def tight_choices_dfs(adj, H, F, W, limit=10**6):
    vs = [v for v in range(len(adj)) if W >> v & 1]
    mem = [H[i] for i in F]
    # BFS order over F by sharing
    order = [0]; seen = {0}
    while len(order) < len(mem):
        for i in range(len(mem)):
            if i not in seen and any(set(mem[i]) & set(mem[j]) for j in order): order.append(i); seen.add(i); break
        else:
            for i in range(len(mem)):
                if i not in seen: order.append(i); seen.add(i); break
    pair_owner = {}
    for i, t in enumerate(mem):
        for p in pairs3(t): pair_owner[p] = i
    # status of pair under assignment: edited? (pair edited iff chosen[owner]==pair)
    chosen = {}
    # triples within W with at least one pair owned: check P3 in G xor S when all owners of owned pairs assigned
    triples_by_owner = {}
    Wl = vs
    import itertools as it
    trip = []
    for x, y, z in it.combinations(Wl, 3):
        owners = {pair_owner[p] for p in ((x, y), (y, z), (x, z)) if p in pair_owner}
        trip.append(((x, y, z), owners))
    # index triples by the last-assigned owner (in order)
    pos = {m: k for k, m in enumerate(order)}
    by_last = {}
    for t, owners in trip:
        if owners: by_last.setdefault(max(pos[o] for o in owners), []).append(t)
        else: by_last.setdefault(-1, []).append(t)
    def edge(a, b):
        e = adj[a] >> b & 1
        p = key(a, b)
        if p in pair_owner and chosen.get(pair_owner[p]) == p: e ^= 1
        return e
    def bad(t):
        x, y, z = t
        e1, e2, e3 = edge(x, y), edge(y, z), edge(x, z)
        return e1 + e2 + e3 == 2
    for t in by_last.get(-1, []):
        if bad(t): return []
    res = []
    def rec(k):
        if len(res) >= limit: return
        if k == len(order): res.append(tuple(chosen[i] for i in range(len(mem)))); return
        m = order[k]
        for p in pairs3(mem[m]):
            chosen[m] = p
            if not any(bad(t) for t in by_last.get(k, [])): rec(k + 1)
            del chosen[m]
    rec(0)
    return res
if __name__ == "__main__":
    import random
    from w6a_rule import tight_choices, conn_subfamilies
    rng = random.Random(5); cnt = 0
    for _ in range(300):
        n = rng.randint(6, 9); adj = random_instance(n, rng); H = maximal_pair_packing(adj, rng)
        if not H: continue
        for F in conn_subfamilies(H, 3):
            W = members_vertices(H, F)
            a = sorted(tuple(sorted(s)) for s in tight_choices(adj, H, F, W))
            b = sorted(tuple(sorted(s)) for s in tight_choices_dfs(adj, H, F, W))
            assert a == b, (adj, H, F); cnt += 1
    print("dfs == brute on", cnt, "subfamilies")
