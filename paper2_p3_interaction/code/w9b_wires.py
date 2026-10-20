import itertools
from w9b_poly import *
# role of q in A: position a/b/c ; pair roles for (q,p1),(q,p2): q=a: (ab,ac) ; q=b: (ab,bc); q=c: (bc,ac)
rp={'a':('ab','ac'),'b':('ab','bc'),'c':('bc','ac')}
names='LRM'
res={}
for qa in 'abc':
    for qb in 'abc':
        for eps in itertools.product((0,1),repeat=4):
            e=[[eps[0],eps[1]],[eps[2],eps[3]]]
            R=cross_relation(rp[qa],rp[qb],e,[(0,0),(0,1),(1,0),(1,1)])
            # bijection between 2-subsets?
            if len(R)==2:
                (a1,b1),(a2,b2)=sorted(R)
                if a1!=a2 and b1!=b2:
                    res.setdefault((qa,qb),set()).add(tuple(sorted(R)))
for k,v in sorted(res.items()):
    print(k,[ ''.join(names[x]+names[y]+' ' for x,y in r) for r in sorted(v)])
