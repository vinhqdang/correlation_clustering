import itertools, sys
from pysat.solvers import Cadical153
def bool_nu(rels, r, extra=None):
    """rels: list of sets of tuples over {0,1}; find f:{0,1}^r->{0,1}, NU, preserving all rels. CEGAR."""
    tuples=list(itertools.product((0,1),repeat=r)); tid={t:i+1 for i,t in enumerate(tuples)}
    S=Cadical153()
    for t in tuples:
        c=sum(t)
        if c<=1: S.add_clause([-tid[t]])
        if c>=r-1: S.add_clause([tid[t]])
    while True:
        if not S.solve(): return None
        m=S.get_model(); f={t:m[tid[t]-1]>0 for t in tuples}
        found=False
        for R in rels:
            k=len(next(iter(R))); Rl=list(R)
            for rows in itertools.product(Rl,repeat=r):
                cols=[tuple(rows[j][i] for j in range(r)) for i in range(k)]
                img=tuple(int(f[c]) for c in cols)
                if img not in R:
                    for bad in itertools.product((0,1),repeat=k):
                        if bad not in R:
                            S.add_clause([(-tid[cols[i]] if bad[i] else tid[cols[i]]) for i in range(k)])
                    found=True; break
            if found: break
        if not found: return f
if __name__=="__main__":
    R5={(0,0,0),(0,1,1),(1,0,1),(1,1,0),(1,1,1)}
    OR3={t for t in itertools.product((0,1),repeat=3) if sum(t)>0}
    for r in range(3,9):
        print(r,'OR3',bool_nu([OR3],r) is not None,'R5',bool_nu([R5],r) is not None,flush=True)
