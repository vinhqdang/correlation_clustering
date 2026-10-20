import sys, itertools, random
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from r1_lib import *
rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 1)
NI=int(sys.argv[2]) if len(sys.argv)>2 else 150
fail={}
cnt={}
def chk(name,cond,info=None):
    cnt[name]=cnt.get(name,0)+1
    if not cond:
        fail.setdefault(name,[]).append(info)
def lab_of(p):
    lab={}
    for i,b in enumerate(p):
        for v in b: lab[v]=i
    return lab
def edits(n,E,lab):
    return {frozenset((a,b)) for a,b in itertools.combinations(range(n),2) if (frozenset((a,b)) in E)!=(lab[a]==lab[b])}
for it in range(NI):
    n=rng.choice([5,6,7,8]); pr=rng.choice([0.2,0.3,0.4,0.5,0.6])
    E=rand_graph(n,pr,rng)
    if not p3s(n,E): continue
    H,used=rand_packing(n,E,rng,maximal=True,keep=rng.choice([0.3,0.7,1.0]))
    if not H: continue
    nH=len(H); VH=set(v for t in H for v in t)
    R=[v for v in range(n) if v not in VH]
    Q={v:set() for v in range(n)}
    for P in used:
        a,b=tuple(P); Q[a].add(b); Q[b].add(a)
    D=max(len(Q[v]) for v in range(n)); q=D
    adj={v:{u for u in range(n) if frozenset((u,v)) in E} for v in range(n)}
    hardadj={v:{u for u in adj[v] if frozenset((u,v)) not in used} for v in range(n)}
    # classes
    seen=set(); classes=[]
    for r in R:
        if r in seen: continue
        comp={r}; st=[r]
        while st:
            x=st.pop()
            for y in adj[x]:
                if y in R and y not in comp: comp.add(y); st.append(y)
        seen|=comp; classes.append(comp)
    # Lemma R
    for r,s in itertools.combinations(R,2):
        if s in adj[r]: chk('LemR',(adj[r]|{r})==(adj[s]|{s}),(n,E))
        else: chk('LemR',not(adj[r]&adj[s]),(n,E))
    A={tuple(sorted(X)):set(adj[next(iter(X))])-set(X) for X in classes}
    ap=allparts_cost(n,E); m=min(c for c,_ in ap)
    opts=[p for c,p in ap if c==m]
    for p in opts:
        lab=lab_of(p); S=edits(n,E,lab); U={e for e in S if e not in used}; u=len(U)
        ell=m-nH
        # C1
        for b in p:
            inter=set(b)&set(R)
            if inter:
                ok=any(inter==X for X in classes)
                chk('C1',ok,(n,E,H,p))
        # N2, N3, P-props
        for X in classes:
            AX=A[tuple(sorted(X))]; x0=next(iter(X)); Bx=set(p[lab[x0]])-set(R)
            if len(X)>ell: chk('N2',set(p[lab[x0]])==set(X)|AX,(n,E,H,p))
            if len(X)>=1:
                for z,z2 in itertools.combinations(AX,2):
                    if z2 not in adj[z]: chk('N3a',frozenset((z,z2)) in used,(n,E,H))
                for a in set(range(n))-set(X)-AX:
                    if a in VH:
                        for z in AX:
                            if z in adj[a]: chk('N3b',frozenset((a,z)) in used,(n,E,H))
            if Bx and not (Bx&AX): chk('P2',False,(n,E,H,p))
            if (Bx-AX): chk('P3',bool(Bx&AX),(n,E,H,p))
            if len(AX)>=ell+q and q>0: chk('P4',Bx<=AX,(n,E,H,p,ell))
        # Lemma half
        for b in p:
            for v in b:
                d=len(adj[v]&set(b)); chk('half_a',2*d>=len(b)-1,(n,E,p))
            for a,c in itertools.combinations(b,2):
                if c not in adj[a]: chk('half_b',any((z in adj[a] and z in adj[c]) for z in b if z not in (a,c)),(n,E,p))
        # T1T2
        for e in U:
            a,b=tuple(e)
            if e not in E:
                C=p[lab[a]]
                ok=False
                for z in C:
                    if z in (a,b): continue
                    if z in adj[a] and z in adj[b]:
                        h1=frozenset((a,z)) not in used; h2=frozenset((z,b)) not in used
                        if not(h1 and h2): ok=True
                chk('T1T2',ok,(n,E,H,p))
    # all partitions with cost<=nH+2 : near, EX, W1,W2,DEL,F1, BAD
    for c,p in ap:
        if c>nH+2: continue
        lab=lab_of(p); S=edits(n,E,lab); U={e for e in S if e not in used}; u=len(U)
        chk('decomp',c==sum(1 for e in S if e in used)+u,None)
        for z,b in itertools.combinations(range(n),2):
            if lab[z]==lab[b]:
                Aset={a for a in range(n) if a not in(z,b) and frozenset((a,z)) not in used and frozenset((a,b)) not in used and ((a in adj[z])!=(a in adj[b]))}
                chk('near',len(Aset)<=u,(n,E,H,p))
        # W1/W2
        tot=0
        for P in used:
            zz,bb=tuple(P)
            same=lab[zz]==lab[bb]
            if same: conf=[y for y in range(n) if y not in (zz,bb) and frozenset((zz,y)) not in used and frozenset((bb,y)) not in used and ((y in adj[zz])!=(y in adj[bb]))]
            else: conf=[y for y in hardadj[zz]&hardadj[bb]]
            for y in conf:
                chk('W1',(frozenset((zz,y)) in U) or (frozenset((bb,y)) in U),(n,E,H,p))
            tot+=len(conf)
        chk('W2',tot<=2*D*u,(n,E,H,p,tot,D,u))
        # EX
        for e in U:
            x,y=tuple(e)
            for (x_,y_) in ((x,y),(y,x)):
                Cx=set(p[lab[x_]])
                dele=lambda v:sum(1 for f in U if v in f and f in E)
                insr=lambda v:sum(1 for f in U if v in f and f not in E)
                if e in E:
                    chk('EXa',len(hardadj[x_])<=D+dele(x_)+dele(y_)-1,(n,E,H,p))
                else:
                    chk('EXa',len(hardadj[x_])<=D+insr(y_)-1+dele(x_),(n,E,H,p))
                chk('EXb',len(Cx)<=1+D+len(hardadj[x_])+insr(x_),(n,E,H,p))
                chk('EXb2',len(Cx)<=1+2*D+u,(n,E,H,p))
                for zc in Cx: chk('EXc',len(hardadj[zc])<=2*D+2*u,(n,E,H,p))
        # DEL
        for e in U:
            if e in E:
                z,y=tuple(e); chk('DEL',u>=1+len(hardadj[z]&hardadj[y]),(n,E,H,p))
        # SM: pieces
        K=[(a,b) for a in range(n) for b in adj[a] if a<b and frozenset((a,b)) in E and frozenset((a,b)) not in used and frozenset((a,b)) not in U]
        par=list(range(n))
        def f(x):
            while par[x]!=x: par[x]=par[par[x]]; x=par[x]
            return x
        for a,b in K: par[f(a)]=f(b)
        for b in p:
            pcs={}
            for v in b: pcs.setdefault(f(v),[]).append(v)
            if len(pcs)>=2:
                chk('SM_iii',len(b)<=max(D+u+1,2*D,2*u),(n,E,H,p))
                for A1,A2 in itertools.combinations(list(pcs.values()),2):
                    for a in A1:
                        for bb in A2:
                            e=frozenset((a,bb)); chk('SM_i',(e in used) or (e in U and e not in E),(n,E,H,p))
        # BAD
        ncl=0
        for t in H:
            x,y,z=t
            states=[({frozenset((x,y))},),({frozenset((y,z))},),({frozenset((x,z))},)]
            allbad=True
            for st in ([frozenset((x,y))],[frozenset((y,z))],[frozenset((x,z))]):
                # merged status of each pair of t under this state
                bad=False
                for pr_ in pairs3(t):
                    a_,b_=tuple(pr_)
                    merged=((frozenset((a_,b_)) in E)!=(pr_ in st))
                    if merged: w_=len([y_ for y_ in range(n) if y_ not in (a_,b_) and frozenset((a_,y_)) not in used and frozenset((b_,y_)) not in used and ((y_ in adj[a_])!=(y_ in adj[b_]))])
                    else: w_=len(hardadj[a_]&hardadj[b_])
                    if w_>0: bad=True
                if not bad: allbad=False
            if allbad: ncl+=1
        cmax=max(sum(1 for t in H if v in t) for v in range(n))
        chk('BAD',ncl<=2*cmax*(c-nH),(n,E,H,p,ncl,c-nH))
    # HUB
    for v in range(n):
        Nh=sorted(hardadj[v])
        edgesG=[(a,b) for a,b in itertools.combinations(Nh,2) if frozenset((a,b)) not in E and frozenset((a,b)) not in used]
        # max matching brute
        def mm(es):
            best=0
            for k in range(1,min(len(es),4)+1):
                for comb in itertools.combinations(es,k):
                    vs=[x for e in comb for x in e]
                    if len(set(vs))==2*k: best=k;break
                else: break
            return best
        nu=mm(edgesG)
        ell=m-nH
        chk('HUBa',nu<=ell,(n,E,H))
        # alpha
        al=0
        for k in range(1,len(Nh)+1):
            found=False
            for comb in itertools.combinations(Nh,k):
                if all(frozenset(x) not in E for x in itertools.combinations(comb,2)): found=True;break
            if found: al=k
            else: break
        cH=max(sum(1 for t in H if v_ in t) for v_ in range(n))
        chk('HUBc',al<=2*ell+cH+1,(n,E,H,v,al,ell,cH))
    # N4 shrink
    for X in classes:
        if len(X)>=2:
            x=max(X)
            # G-x
            keep=[v for v in range(n) if v!=x]; mp={v:i for i,v in enumerate(keep)}
            E2={frozenset((mp[a],mp[b])) for e in E for a,b in [tuple(e)] if a!=x and b!=x}
            o2=opt(n-1,E2)
            for ell in range(0,4):
                if len(X)>=ell+2:
                    chk('N4',(m<=nH+ell)==(o2<=nH+ell),(n,E,H,ell,X))
print('checks',sum(cnt.values()))
for k in sorted(cnt): print(k,cnt[k],'FAIL' if k in fail else 'ok', (fail[k][:1] if k in fail else ''))
