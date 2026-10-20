import sys, time
from w10_lib import *
m=10
W=World(2,thin_mode=True)
A=W.add_chain('A',m); Bc=W.add_chain('B',m)
host=('s',10)
print(A.perm_guest('zA',host,{4,5}), Bc.perm_guest('zB',host,{4,5}))
print(cmp_guest(A,host,{a for a in range(m) if a>3},Bc,host,{a for a in range(m) if a<6},'y','zA','zB'))
I=ThinInst(W); print('n',I.n,I.valid)
t0=time.time(); bad=0
for a in range(m):
    row=''
    for b in range(m):
        lab=I.solve(together=I.force_pos(A,a)+I.force_pos(Bc,b))
        exp=(a not in (4,5)) and (b not in (4,5)) and (a<=3 or b>=6)
        row+='#' if lab is not None else '.'
        bad+= (lab is not None)!=exp
    print(a,row)
print('mismatch',bad,round(time.time()-t0,1))
