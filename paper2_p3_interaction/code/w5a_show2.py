import numpy as np, itertools
from w5a_lib import *
from w5a_loc import w2local
n=8
E=[(0, 2), (0, 6), (1, 4), (1, 6), (1, 7), (2, 3), (2, 5), (3, 7), (4, 5), (4, 6), (4, 7), (5, 7)]
H=[(0, 6, 4), (3, 7, 4), (2, 5, 4), (0, 2, 3), (6, 1, 7)]
A=np.zeros((n,n),np.int8)
for i,j in E: A[i,j]=A[j,i]=1
pk=packed_matrix(n,H); L=w2local(pk)
lab_t,same,P=part_table(n); lab,_,cost=all_costs(A); opt=cost.min()
print('opt',opt,'maximal',is_pair_maximal(A,H))
for i in np.nonzero(cost==opt)[0]:
    l=lab_t[i]; cl={}
    for v in range(n): cl.setdefault(l[v],[]).append(v)
    pid,U,F,I=pieces_of(A,pk,lambda a,b:l[a]==l[b])
    r=anatomy(A,pk,l)
    print(list(cl.values()),'U',U,'F',F,'I',I,'t1,t2',r['t1'],r['t2'],'nonW2local inserted',[p for p in I if not L[p]])
