# audit9: Lemma tmaj (cores have <=2 anchors if Sigma closed under some majority op), Example ex:maj
import sys, random, itertools
from r2_lib import *
perms=list(itertools.permutations((1,2,3)))
def states(H,E,lab,n):
    st=[]
    for (a,b,c) in H:
        tab=lab[a]==lab[b]; tbc=lab[b]==lab[c]; tac=lab[a]==lab[c]
        # edits: ab edited iff not together; bc likewise; ac edited iff together
        ed=[not tab, not tbc, tac]
        if sum(ed)!=1: return None
        st.append(ed.index(True)+1)
    return tuple(st)
def sigma(n,E,H):
    cs=costs(n,E); L,prs,T=parts(n)
    S=set(); parts_=[]
    for k in np.nonzero(cs==len(H))[0]:
        s=states(H,E,L[k],n)
        if s is not None: S.add(s); parts_.append(L[k])
    return S,parts_
def closed(S,op):
    # op: dict perm->value
    S=list(S)
    if not S: return True
    m=len(S[0]); Ss=set(S)
    for x,y,z in itertools.product(S,repeat=3):
        r=[]
        for i in range(m):
            a,b,c=x[i],y[i],z[i]
            if a==b or a==c: r.append(a)
            elif b==c: r.append(b)
            else: r.append(op[(a,b,c)])
        if tuple(r) not in Ss: return False
    return True
def some_majority(S):
    if len(S)<=1: return True
    # restrict search: ops 3^6; prune using necessary condition via small search
    for vals in itertools.product((1,2,3),repeat=6):
        op=dict(zip(perms,vals))
        if closed(S,op): return True
    return False
if __name__=="__main__":
    # Example ex:maj
    n=6; E=set(frozenset(map(int,e)) for e in "03 05 12 14 23 24 25 34 35 45".split())
    H=[(0,5,4),(1,2,5),(0,3,2),(1,4,3)]
    S,_=sigma(n,E,H); print("ex:maj Sigma",sorted(S))
    op3={p:3 for p in perms}; op1={p:1 for p in perms}
    print("closed m3",closed(S,op3),"closed m1",closed(S,op1))
    # random test of lemma tmaj
    seed=int(sys.argv[1]); N=int(sys.argv[2]); rng=random.Random(seed)
    cnt=dict(tight=0,maj=0,nomaj=0,core_ge3_maj=0,core_ge3_nomaj=0,maxcore=0)
    for it in range(N):
        n=rng.randint(6,9); p=rng.choice([0.3,0.45,0.6])
        E=rand_graph(n,p,rng); H,used=rand_packing(n,E,rng,maximal=True)
        if not H: continue
        if int(costs(n,E).min())!=len(H): continue
        cnt['tight']+=1
        S,pl=sigma(n,E,H)
        hasH=lambda v: any(frozenset((v,u)) in E and frozenset((v,u)) not in used for u in range(n) if u!=v)
        V=sorted({v for t in H for v in t}); anchors=[v for v in V if not hasH(v)]
        if len(anchors)>8 or len(anchors)<3: continue
        single=[ {v for v in range(n) if all(lab[v]!=lab[u] for u in range(n) if u!=v)} for lab in pl]
        def face(Sa): return any(set(Sa)<=s for s in single)
        cores=[]
        for k in range(0,len(anchors)+1):
            for Sa in itertools.combinations(anchors,k):
                if not face(Sa) and all(face(Sa[:i]+Sa[i+1:]) for i in range(len(Sa))): cores.append(Sa)
        mx=max([len(c) for c in cores],default=0)
        cnt['maxcore']=max(cnt['maxcore'],mx)
        if mx>=3:
            m=some_majority(S)
            cnt['core_ge3_maj' if m else 'core_ge3_nomaj']+=1
            if m: print("VIOL tmaj",seed,n,sorted(map(sorted,E)),H,cores)
    print(seed,cnt)
