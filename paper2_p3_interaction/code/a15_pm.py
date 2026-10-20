import sys
from a15_lib import *
def check(I):
    soft = set()
    for t in I.members: soft |= pairs_of(t)
    n = I.n
    hn = [set() for _ in range(n)]
    for (u, v) in I.E:
        if (u, v) not in soft:
            hn[u].add(v); hn[v].add(u)
    bad = 0; ex = None
    for v in range(n):
        nb = sorted(hn[v])
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                u, w = nb[i], nb[j]
                pq = key(u, w)
                if pq not in soft and pq not in I.E:
                    bad += 1
                    if ex is None: ex = [I.W.B.names[x] for x in (u, v, w)]
    return bad, ex
if __name__ == '__main__':
    for (k, n, p) in [(1, 2, 8), (2, 2, 10)]:
        S = [[{(1, 1), (2, 1)} for _ in range(k)] for _ in range(k)]
        W, meta = build_gt(k, n, S, p, thin_mode=False)
        I = Inst(W)
        print(k, n, p, check(I), flush=True)
