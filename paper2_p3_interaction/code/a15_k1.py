import sys
from a15_lib import *
p = int(sys.argv[1]); mode = sys.argv[2]
n = 1
U, m, sg, sep = layout2(n, 1)
tiles = {sg[(1, 1)]} if mode == 'yes' else set()
tun = [((('A', 0, 0), sg[(1, 1)]), (('B', 0, 0), sg[(1, 1)]))] if mode == 'yes' else []
W, meta = build_general(n, p, {('A', 0, 0): set(tiles), ('B', 0, 0): set(tiles)}, [], tun, extra=1)
I = Inst(W)
print('n', I.n, '|H|', len(I.members), 'valid', I.valid, 'c', I.c, 'u_max', max(W.used.values()))
N, size, e, ne, cid = contract(I)
print('classes', N)
cost, lab = wopt(N, size, e, ne)
print('contracted OPT', cost, 'excess over |H|', cost - len(I.members))
