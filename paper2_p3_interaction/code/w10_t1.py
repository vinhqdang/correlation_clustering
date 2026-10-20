"""Test 1: comparator guest between two chains (m=10).  separator group u=2 (positions 4,5; host s_10)."""
import sys, time
from w10_lib import *
m=10; p=int(sys.argv[1]) if len(sys.argv)>1 else 8
W=World(p)
A=W.add_chain('A',m); Bc=W.add_chain('B',m)
host=('s',10)
effz=[A.perm_guest('zA',host,{4,5}), Bc.perm_guest('zB',host,{4,5})]
print('perm eff forbid',effz)
effs=cmp_guest(A,host,{a for a in range(m) if a>3},Bc,host,{a for a in range(m) if a<6},'y','zA','zB')
print('cmp eff forbid A,B',effs)
t0=time.time()
I=Inst(W); print('n',I.n,'|H|',len(I.members),I.valid,'c',I.c,round(time.time()-t0,1),flush=True)
bad=0
for a in range(m):
    row=''
    for b in range(m):
        lab=I.solve(together=I.force_pos(A,a)+I.force_pos(Bc,b))
        feas=lab is not None
        exp=(a not in (4,5)) and (b not in (4,5)) and (a<=3 or b>=6)
        row+= ('#' if feas else '.')
        if feas!=exp: bad+=1
    print(a,row,flush=True)
print('mismatches vs expectation',bad,round(time.time()-t0,1))
