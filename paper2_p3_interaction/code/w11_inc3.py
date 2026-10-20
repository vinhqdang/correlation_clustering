import sys
from w11_lib import *
from w11_gen import *
nv,H,npend,cl=(10, [(2, 1, 5), (7, 4, 8), (1, 6, 3), (6, 0, 4), (5, 7, 9), (0, 2, 3)], 5, [3, 1, 5, 3, 4, 3, -1, -1, 1, 0, 4, 0, 5, 6, 2])
I=make_inst(nv,H,npend,cl)
print('n',I.n,'valid',I.valid(),'pm',I.pair_maximal(),'c',I.cmax(),'IH edges',sorted(I.IH().edges()))
print('|H|',len(H),'remote pairs',len(I.remote()),'OPT?', opt(I)-len(H) if I.n<=16 else '')
Phi=set(); 
for step in range(6):
    sz,C,t=min_increment(I,Phi,workers=4)
    print('step',step,'F\'(Phi)-|H| =',t-len(H),'min increment size',sz,C)
    if sz is None: break
    Phi|=set(C)
print('F(N)-|H| =',fprime_value(I,set(I.remote()))-len(H))
for l in range(0,4):
    s,T=min_certificate(I,l); print('min certificate for l=%d:'%l,s,T)
for l in range(4,6):
    s,T=min_certificate(I,l); print('min certificate for l=%d:'%l,s,T)
