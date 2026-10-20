"""RM2 (abstract model of wave 7): forest, conflicts = pairs of nodes (path), operations: break (node, hits every conflict whose path contains it)
and tunnel (pair of nodes from a link graph, hits every conflict with an endpoint in the pair).  Check the simple bounded search tree
(deepest-LCA conflict: break at its LCA, or one of the <= 2d tunnels at its two endpoints) against exhaustive search."""
import random, itertools, sys

def gen(rng, N, nconf, ntun, deg):
    par = [-1] * N
    for v in range(1, N):
        par[v] = rng.randrange(v) if rng.random() < 0.9 else -1
    depth = [0] * N
    for v in range(N):
        d = 0; u = v
        while par[u] != -1:
            u = par[u]; d += 1
        depth[v] = d
    def root(v):
        while par[v] != -1: v = par[v]
        return v
    def path(x, y):
        if root(x) != root(y): return None
        a, b = x, y; pa = [a]; pb = [b]
        while a != b:
            if depth[a] >= depth[b]:
                a = par[a]; pa.append(a)
            else:
                b = par[b]; pb.append(b)
        # a == b == lca
        i = pa.index(a)
        j = pb.index(b)
        return pa[:i + 1] + pb[:j][::-1], a
    conf = []
    tries = 0
    while len(conf) < nconf and tries < 1000:
        tries += 1
        x, y = rng.randrange(N), rng.randrange(N)
        p = path(x, y)
        if p is None: continue
        conf.append((x, y, tuple(p[0]), p[1]))
    tun = []
    dg = [0] * N
    tries = 0
    while len(tun) < ntun and tries < 1000:
        tries += 1
        x, y = rng.randrange(N), rng.randrange(N)
        if x == y or dg[x] >= deg or dg[y] >= deg: continue
        if (x, y) in tun or (y, x) in tun: continue
        tun.append((x, y)); dg[x] += 1; dg[y] += 1
    return par, depth, conf, tun

def exhaustive(N, conf, tun, L):
    best = 10**9
    for k in range(0, len(tun) + 1):
        for F in itertools.combinations(tun, k):
            ends = set(v for t in F for v in t)
            rem = [c for c in conf if c[0] not in ends and c[1] not in ends]
            # min hitting set of node paths by exhaustive search over breaks
            nodes = sorted(set(v for c in rem for v in c[2]))
            if not rem:
                b = 0
            else:
                b = None
                for kk in range(0, len(rem) + 1):
                    for B in itertools.combinations(nodes, kk):
                        Bs = set(B)
                        if all(Bs & set(c[2]) for c in rem):
                            b = kk; break
                    if b is not None: break
            best = min(best, k + b)
    return best

def branching(depth, conf, tun, L):
    """return min cost <= L via search tree, or None if > L"""
    tunadj = {}
    for (x, y) in tun:
        tunadj.setdefault(x, []).append((x, y)); tunadj.setdefault(y, []).append((x, y))
    best = [None]
    def rec(unhit, cost):
        if best[0] is not None and cost >= best[0]: return
        if cost > L: return
        if not unhit:
            best[0] = cost; return
        # deepest lca
        q = max(unhit, key=lambda c: depth[c[3]])
        lca = q[3]
        # option 1: break at lca
        rec([c for c in unhit if lca not in c[2]], cost + 1)
        # option 2: tunnels at endpoints
        for e in (q[0], q[1]):
            for t in tunadj.get(e, []):
                ends = set(t)
                rec([c for c in unhit if c[0] not in ends and c[1] not in ends], cost + 1)
    rec(list(conf), 0)
    return best[0]

if __name__ == '__main__':
    seed = int(sys.argv[1]); T = int(sys.argv[2])
    rng = random.Random(seed)
    bad = 0; withtun = 0
    for it in range(T):
        N = rng.randint(4, 9)
        par, depth, conf, tun = gen(rng, N, rng.randint(1, 6), rng.randint(0, 5), 2)
        ex = exhaustive(N, conf, tun, 10)
        br = branching(depth, conf, tun, 10)
        if ex != br:
            bad += 1; print('MISMATCH', ex, br, par, conf, tun)
        # does the tunnel ever help?
        exn = exhaustive(N, conf, [], 10)
        if exn > ex: withtun += 1
    print('instances', T, 'mismatches', bad, 'tunnel strictly helps in', withtun)
