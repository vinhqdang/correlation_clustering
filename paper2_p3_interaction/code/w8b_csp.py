"""CSP-style relaxation: triples containing a non-local pair outside T are dropped (not the partition relaxation).
F_csp(T) = min cost over x in {0,1}^pairs with triangle constraints only on triples all of whose pairs are counted or in T;
T pairs fixed to their G-status (respected); counted pairs cost 1 per disagreement; free pairs unconstrained."""
import numpy as np, itertools
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from w8b_ilp import ILP
class ILPcsp(ILP):
    def __init__(self,inst):
        super().__init__(inst)
        n=inst.n
        self.trip=[]
        for a,b,c in itertools.combinations(range(n),3):
            self.trip.append((self.idx[(a,b)],self.idx[(b,c)],self.idx[(a,c)]))
    def F(self,flip,T,budget=None):
        m=len(self.pr); Tset=set(T)
        free=set(i for i in self.N if i not in Tset)
        c=np.zeros(m); lb=np.zeros(m); ub=np.ones(m); const=0
        for i,p in enumerate(self.pr):
            e=(p in self.I.E)^bool(flip.get(i,False))
            if i in Tset: lb[i]=ub[i]=1.0 if e else 0.0
            elif self.counted[i]:
                if e: c[i]=-1; const+=1
                else: c[i]=1
        rows=[]
        for (ab,bc,ac) in self.trip:
            if ab in free or bc in free or ac in free: continue
            rows+=[(ab,bc,ac),(ab,ac,bc),(bc,ac,ab)]
        A=lil_matrix((len(rows),m))
        for r,(x,y,z) in enumerate(rows): A[r,x]=1;A[r,y]=1;A[r,z]=-1
        res=milp(c,constraints=LinearConstraint(A.tocsr(),-np.inf,np.ones(len(rows))),integrality=np.ones(m),bounds=Bounds(lb,ub))
        if res.status!=0: return 10**9
        return int(round(res.fun+const))
if __name__=="__main__":
    from w8b_ilp import far
    for m in (3,4,5,6):
        I=far(m); o=ILPcsp(I); r,st=o.run(1); print("FARcsp",m,r,st['cores'],st['nodes'])
    from w8b_nstar import nstar
    I=nstar([(0,'*'),(0,1),(0,2),(1,3),(1,2),(3,0),(2,3)],0); o=ILPcsp(I); r,st=o.run(1); print("N*K4 csp",r,st)
