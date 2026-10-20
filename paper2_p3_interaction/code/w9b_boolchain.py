import itertools
R5={(0,0,0),(0,1,1),(1,0,1),(1,1,0),(1,1,1)}
def chain_sol(N,t0,tN):
    """solutions x in {0,1}^N for which exist t_1..t_{N-1} with R5(t_{i-1},x_i,t_i), t_0=t0,t_N=tN"""
    sols=set()
    for x in itertools.product((0,1),repeat=N):
        # DP
        S={t0}
        for i in range(N):
            S2=set()
            for t in S:
                for t2 in (0,1):
                    if (t,x[i],t2) in R5: S2.add(t2)
            S=S2
        if tN in S: sols.add(x)
    return sols
def max_core(sols,N):
    # lists: singleton-value lists (x_i=v). find minimal inconsistent families
    lits=[(i,v) for i in range(N) for v in (0,1)]
    best=0;bestf=None
    for k in range(1,N+1):
        for fam in itertools.combinations(lits,k):
            if len({i for i,v in fam})<k: continue
            def cons(f): return any(all(s[i]==v for i,v in f) for s in sols)
            if cons(fam): continue
            if all(cons(tuple(x for x in fam if x!=y)) for y in fam):
                if k>best: best=k;bestf=fam
    return best,bestf
for t0 in (0,1):
    for tN in (0,1):
        out=[]
        for N in range(2,8):
            s=chain_sol(N,t0,tN); out.append((N,len(s),max_core(s,N)[0]))
        print('t0',t0,'tN',tN,out)
