import sys, time
from w11_lib import *
from w11_cell import *
from w11_gen import make_inst
cell={'k': 3, 'groups': [[(0, 0), (1, 1)], [(2, 1), (1, 0)], [(2, 0), (0, 1)], [(2, 2)]], 'out': [(1, 2)], 'inn': [(0, 2)], 'centres': [0, 0, 2], 'perm': [[2, 1, 0], [2, 1, 0], [2, 1, 0]], 'pend': 2, 'nlab': 3, 'lab': None}
la={('g', 1): (0, 0), ('g', 2): (0, 0), ('g', 3): (2, 1), ('w', 1): (1, 0), ('ws', 0): (1, 0), ('ws', 1): (0, 1), ('we', 1): (0, 0)}
l=int(sys.argv[1]); Ns=[int(x) for x in sys.argv[2].split(',')]
for N in Ns:
    nv,H,npd,cl=build_chain(cell,N,la,None,[0,0],[0])
    # start pendants: sn=2, end pendants en=1 as in the cell search run? infer from labels: ws 0,1 and we 1 -> en=2? try both
    for en in (1,2):
        nv,H,npd,cl=build_chain(cell,N,la,None,[0,0],[0]*en)
        I=make_inst(nv,H,npd,cl)
        if not (I.valid() and I.pair_maximal() and I.cmax()<=2): print('N',N,'en',en,'invalid'); continue
        t0=time.time()
        sz,T=min_certificate(I,l,max_size=60,workers=4)
        print('N',N,'en',en,'n',I.n,'members',len(H),'l',l,'min cert',sz,round(time.time()-t0,1),'s',flush=True)
