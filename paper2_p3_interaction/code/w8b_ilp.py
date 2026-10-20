"""ILP oracle version of the far-pair branching (scipy milp). Used for larger n (FAR_m, N*(F))."""
import numpy as np, itertools, sys
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from w8b_lib import Inst

class ILP:
    def __init__(self,inst):
        self.I=inst; n=inst.n; self.n=n
        self.pr=inst.pairs; self.idx=inst.pidx
        self.loc=inst.local_matrix(); self.width=inst.width
        self.counted=[(p in inst.soft) or self.loc[p] for p in self.pr]
        self.N=[i for i,c in enumerate(self.counted) if not c]
        rows=[]
        for a,b,c in itertools.combinations(range(n),3):
            ab,bc,ac=self.idx[(a,b)],self.idx[(b,c)],self.idx[(a,c)]
            rows.append((ab,bc,ac)); rows.append((ab,ac,bc)); rows.append((bc,ac,ab))
        A=lil_matrix((len(rows),len(self.pr)))
        for r,(x,y,z) in enumerate(rows):
            A[r,x]=1;A[r,y]=1;A[r,z]=-1   # x_xy + x_yz - x_xz <= 1
        self.A=A.tocsr(); self.ub=np.ones(len(rows))
    def F(self,flip,T,budget=None):
        """min cost over partitions; x_p=1 iff together. cost of counted pair: [edge xor flip] ? (1-x):x ; T pairs forced."""
        m=len(self.pr); c=np.zeros(m); lb=np.zeros(m); ubv=np.ones(m); const=0
        Tset=set(T)
        for i,p in enumerate(self.pr):
            e=(p in self.I.E)^bool(flip.get(i,False))
            if i in Tset:
                lb[i]=ubv[i]=1.0 if e else 0.0
            elif self.counted[i]:
                if e: c[i]=-1; const+=1
                else: c[i]=1
        res=milp(c,constraints=LinearConstraint(self.A,-np.inf,self.ub),integrality=np.ones(m),bounds=Bounds(lb,ubv))
        if res.status!=0: return 10**9
        return int(round(res.fun+const))
    def solve(self,flip,l,stats):
        B=len(self.I.H)+l; stats['nodes']+=1
        if self.F(flip,self.N)<=B: return True
        if l==0: return False
        T=list(self.N)
        for q in list(T):
            T2=[x for x in T if x!=q]
            if self.F(flip,T2)>B: T=T2
        stats['cores'].append(len(T))
        for e in T:
            f2=dict(flip); f2[e]=not f2.get(e,False)
            if self.solve(f2,l-1,stats): return True
        return False
    def run(self,l):
        st={'nodes':0,'cores':[]}; r=self.solve({},l,st); return r,st

def far(m,extra=((1,None),)):
    n=2*m+1; E=set((i,i+1) for i in range(2*m)); E.add((1,2*m-1))
    H=[(2*i,2*i+1,2*i+2) for i in range(m)]
    return Inst(n,E,H)
if __name__=="__main__":
    for m in range(3,9):
        I=far(m); assert I.valid()
        o=ILP(I)
        r,st=o.run(1)
        print("FAR",m,"n",I.n,"width",o.width,"nonlocal",len(o.N),"l=1:",r,"cores",st['cores'],"nodes",st['nodes'])
        r0,st0=o.run(0); print("   l=0:",r0)
