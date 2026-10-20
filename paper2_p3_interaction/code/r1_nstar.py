import sys, itertools; sys.path.insert(0,'.')
from r1_lib import pairs3
def build(Fedges, t):
    # Fedges: list of oriented (tail,head) on vertex names; add (t,'*')
    edges=[(t,'*')]+list(Fedges)
    N=len(edges)-1
    V={}
    def vid(name):
        if name not in V: V[name]=len(V)
        return V[name]
    a=[vid(('a',i)) for i in range(N+1)]; b=[vid(('b',i)) for i in range(N+1)]
    c=[vid(('c',i)) for i in range(N+1)]; u=[vid(('u',i)) for i in range(N+1)]
    z=vid('z'); y=vid('y'); u_ext=u+[z]
    E=set()
    verts=set(x for e in edges for x in e)
    for v in verts:
        K=[a[i] for i,(tl,hd) in enumerate(edges) if tl==v]+[x for i,(tl,hd) in enumerate(edges) if hd==v for x in (b[i],c[i],u[i])]
        for p in itertools.combinations(K,2): E.add(frozenset(p))
    for i in range(N+1): E.add(frozenset((a[i],b[i])))
    for i in range(N+1): E.add(frozenset((b[i],u_ext[i+1])))
    E.add(frozenset((y,u[0]))); E.add(frozenset((y,c[0])))
    H=[]
    for i in range(N+1): H.append((a[i],b[i],c[i]))
    for i in range(N+1): H.append((u[i],b[i],u_ext[i+1]))
    return len(V),E,H,y
def is_cluster(n,E,S,vs):
    # induced on vs after flipping S
    A=lambda p: (frozenset(p) in E)!=(frozenset(p) in S)
    for x,yv,zv in itertools.permutations(vs,3):
        if x<zv and A((x,yv)) and A((yv,zv)) and not A((x,zv)): return False
    return True
def tight_sols(H,members):
    for st in itertools.product(range(3),repeat=len(members)):
        S=set()
        for m,s in zip(members,st):
            x,yv,zv=H[m]
            S.add(frozenset([(x,yv),(yv,zv),(x,zv)][s]))
        yield S
def run(Fe,t,name):
    n,E,H,y=build(Fe,t)
    members=list(range(len(H)))
    # pair set validity
    used=set()
    for tt in H:
        x,yv,zv=tt
        assert frozenset((x,yv)) in E and frozenset((yv,zv)) in E and frozenset((x,zv)) not in E
        P=pairs3(tt); assert not(P&used); used|=P
    allv=range(n)
    # pair maximality
    pm=[tt for tt in itertools.permutations(allv,3) if tt[0]<tt[2] and frozenset((tt[0],tt[1])) in E and frozenset((tt[1],tt[2])) in E and frozenset((tt[0],tt[2])) not in E and not(pairs3(tt)&used)]
    print(name,'n',n,'|H|',len(H),'pair-maximal',not pm)
    # tight solution of whole graph?
    anyok=False
    for S in tight_sols(H,members):
        if is_cluster(n,E,S,list(allv)): anyok=True;break
    print(' tight solution exists on whole G (would make OPT=|H|):',anyok)
    # tight on W=V-y
    W=[v for v in allv if v!=y]
    anyW=any(is_cluster(n,E,S,W) for S in tight_sols(H,members))
    print(' tight on V\\{y}:',anyW)
    # pair-closed subfamilies, irreducibility
    Wof=lambda F:set(v for m in F for v in H[m])
    nred=0;tot=0;bad=[]
    for r in range(1,len(H)+1):
        for F in itertools.combinations(members,r):
            WF=Wof(F)
            outside=[m for m in members if m not in F]
            if any(len(set(H[m])&WF)>=2 for m in outside): continue
            tot+=1
            red=False
            for S in tight_sols(H,F):
                vs=sorted(WF)
                if not is_cluster(n,E,S,vs): continue
                # clusters of G[WF]+S
                A=lambda p: (frozenset(p) in E)!=(frozenset(p) in S)
                clusters=[];seen=set()
                for v in vs:
                    if v in seen: continue
                    cl={x for x in vs if x==v or A((x,v))}; seen|=cl; clusters.append(frozenset(cl))
                ok=True
                for yy in allv:
                    if yy in WF: continue
                    Ny={x for x in WF if frozenset((x,yy)) in E}
                    if Ny and frozenset(Ny) not in clusters: ok=False;break
                if ok: red=True;break
            if red: nred+=1; bad.append(F)
    print(' pair-closed families',tot,'reducible',nred,bad[:3])
# triangle F, cyclic
run([('v0','v1'),('v1','v2'),('v2','v0')],'v0','triangle')
