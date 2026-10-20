"""NU / majority search over the scope language (numpy-checked CEGAR)."""
import pickle, itertools, sys, time
import numpy as np
from pysat.solvers import Cadical153
D=(0,1,2)
def load(which):
    P,T=pickle.load(open('/tmp/w9b_scopes.pkl','rb'))
    rels=[]   # list of (arity, np.array rows (m x k))
    for L in [(0,),(1,),(2,),(0,1),(0,2),(1,2)]:
        rels.append((1,np.array([[x] for x in L])))
    if which in ('pairs','all','pairs_only'):
        for R in P: 
            if len(R): rels.append((2,np.array(sorted(R))))
    if which in ('all','tri'):
        for R in T:
            if len(R): rels.append((3,np.array(sorted(R))))
    return rels
class Search:
    def __init__(self,r,rels,conservative=True,nu=True,extra=None):
        self.r=r; self.rels=rels
        self.tuples=np.array(list(itertools.product(D,repeat=r)))
        self.n=len(self.tuples)
        self.pw=3**np.arange(r-1,-1,-1)
        self.S=Cadical153()
        v=lambda i,a:3*i+a+1
        self.v=v
        for i in range(self.n):
            self.S.add_clause([v(i,a) for a in D])
            for a,b in itertools.combinations(D,2): self.S.add_clause([-v(i,a),-v(i,b)])
        for i,t in enumerate(self.tuples):
            t=t.astype(int)
            s=set(t.tolist())
            if len(s)==1: self.S.add_clause([v(i,int(t[0]))])
            if conservative:
                for a in D:
                    if a not in s: self.S.add_clause([-v(i,int(a))])
            if nu:
                from collections import Counter
                a,m=Counter(t.tolist()).most_common(1)[0]
                if m>=r-1: self.S.add_clause([v(i,int(a))])
        if extra: extra(self)
        self.added=0
    def check(self,f,maxv=30):
        r=self.r; viol=[]
        for ri,(k,R) in enumerate(self.rels):
            m=len(R)
            if m**r>6_000_000: 
                # chunk over first index
                pass
            Rset=set(map(tuple,R.tolist()))
            # all r-tuples of rows
            grids=np.indices((m,)*r).reshape(r,-1)    # r x m^r
            rowsv=R[grids]                             # r x N x k
            codes=(rowsv*self.pw[:,None,None]).sum(0)  # N x k
            img=f[codes]                               # N x k
            code_img=(img*(3**np.arange(k-1,-1,-1))).sum(1)
            Rcodes=np.zeros(3**k,bool)
            for t in Rset: Rcodes[sum(a*3**(k-1-j) for j,a in enumerate(t))]=True
            bad=np.nonzero(~Rcodes[code_img])[0]
            if len(bad):
                for b in bad[:maxv]:
                    viol.append((ri,codes[b].copy()))
                if len(viol)>=maxv*3: break
        return viol
    def add_viol(self,ri,cols):
        k,R=self.rels[ri]
        Rset=set(map(tuple,R.tolist()))
        for bad in itertools.product(D,repeat=k):
            if bad not in Rset:
                self.S.add_clause([-self.v(int(cols[i]),bad[i]) for i in range(k)])
        self.added+=1
    def solve(self,maxit=10**6,verbose=True):
        it=0
        while it<maxit:
            it+=1
            if not self.S.solve(): return None
            mod=self.S.get_model()
            f=np.zeros(self.n,int)
            for i in range(self.n):
                for a in D:
                    if mod[3*i+a]>0: f[i]=a
            viol=self.check(f)
            if not viol: return f
            for ri,cols in viol: self.add_viol(ri,cols)
            if verbose and it%20==0: print(' it',it,'added',self.added,flush=True)
        return 'timeout'
if __name__=='__main__':
    r=int(sys.argv[1]); which=sys.argv[2]
    rels=load(which); print('relations',len(rels),flush=True)
    t=time.time(); s=Search(r,rels); res=s.solve()
    print('NU arity',r,which,'->','NONE' if res is None else 'FOUND','clauses',s.added,'t',round(time.time()-t,1))
    if res is not None and res is not None: print(list(res))
