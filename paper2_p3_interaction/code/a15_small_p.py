import sys
from a15_lib import *
from w10_exact import exact_decide
p = int(sys.argv[1]); mode = sys.argv[2]; lmax = int(sys.argv[3])
n = 1
U, m, sg, sep = layout2(n, 1)
tiles = {sg[(1, 1)]} if mode == 'yes' else set()
tun = [((('A', 0, 0), sg[(1, 1)]), (('B', 0, 0), sg[(1, 1)]))] if mode == 'yes' else []
W, meta = build_general(n, p, {('A', 0, 0): set(tiles), ('B', 0, 0): set(tiles)}, [], tun, extra=1)
I = Inst(W)
print('p', p, 'n', I.n, '|H|', len(I.members), 'valid', I.valid, 'u_max', max(W.used.values()), flush=True)
cls, cid = classes_of(I)
inv = {v: k for k, v in cid.items()}
for l in range(int(sys.argv[4]) if len(sys.argv) > 4 else 0, lmax + 1):
    r, sol = exact_decide(I.n, I.E, I.members, l, verbose=False, timeout=3000)
    print('l', l, 'OPT<=|H|+l:', r, flush=True)
    if r is True:
        lab, cost = sol
        print('  cost', cost, 'excess', cost - len(I.members))
        split = {}
        for v in range(I.n):
            split.setdefault(cls[v], set()).add(lab[v])
        print('  split classes:', [(inv[c], len(s)) for c, s in split.items() if len(s) > 1 and isinstance(inv[c], tuple) and len(inv[c]) == 2 and inv[c][1][0] in ('s', 'e')])
        break
