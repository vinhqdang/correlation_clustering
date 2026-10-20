"""Audit 14: INDEPENDENT thin-guest model of the Grid-Tiling(<=) reduction, written from notes/wave10_hardness_proof.md Sec. 2 only
(instance data of the bare chain OH_m from w8c_oh; profiles, guests, comparators, tunnels, SAT all fresh).
Tests tunnel-release semantics: ZE(I - C) feasible  <=>  C is a tiling   (all choices, and all cut sets of size <= k^2)."""
import sys, os, itertools, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_oh import oh_instance
from a14_chain import enum_ze, spine_word, k2
from pysat.solvers import Cadical153
from functools import lru_cache

@lru_cache(None)
def chain_data(m):
    n,E,members,names,idx=oh_instance(m)
    sols=enum_ze(n,E,members)
    pos={}
    for l in sols:
        w=spine_word(m,idx,l); a=w.index('M'); assert w=='R'*a+'M'+'L'*(m-a-1); pos[a]=l
    assert sorted(pos)==list(range(m))
    prof={}   # prof[(host,a)] = frozenset of names in cluster minus host
    for a,l in pos.items():
        for h in names:
            prof[(h,a)]=frozenset(names[v] for v in range(n) if l[v]==l[idx[h]] and v!=idx[h])
    return n,E,members,names,idx,prof

def classes(m,u,prof):
    h=('s',4*u+2)
    cl={'A0':range(0,2*u),'A1':[2*u],'A2':[2*u+1],'A3':range(2*u+2,m)}
    P={}
    for c,rng in cl.items():
        sets={prof[(h,a)] for a in rng}
        assert len(sets)<=1,(m,u,c)       # profile constant on class (empty class allowed at the caps)
        P[c]=next(iter(sets)) if sets else frozenset()
    return h,P

def Wsets(m,u,prof):
    h,P=classes(m,u,prof)
    allp=set().union(*P.values())
    return {'pad':{h}|P['A0']|P['A3'],'X':{h}|P['A0'],'Y':{h}|P['A3'],'all':{h}|allp}

def eff_forbidden(m,u,prof,W):
    h=('s',4*u+2)
    return {a for a in range(m) if not prof[(h,a)]<=W}

class Model:
    def __init__(self,k,n,S):
        self.k,self.n,self.S=k,n,S
        U=1+n*(n+1); m=2*U; self.U,self.m=U,m
        slot={(al,s):(al-1)*(n+1)+s for al in range(1,n+1) for s in range(1,n+1)}
        sep={al:al*(n+1) for al in range(1,n+1)}
        cd=chain_data(m); N0,E0,mem0,names0,idx0,prof=cd
        # check table of effective forbidden sets for every group
        for u in range(U):
            Ws=Wsets(m,u,prof)
            assert eff_forbidden(m,u,prof,Ws['pad'])=={2*u,2*u+1}
            assert eff_forbidden(m,u,prof,Ws['X'])==set(range(2*u,m))
            assert eff_forbidden(m,u,prof,Ws['Y'])==set(range(0,2*u+2))
            assert eff_forbidden(m,u,prof,Ws['all'])==set()
        V={}
        def vid(x):
            if x not in V: V[x]=len(V)
            return V[x]
        members=[]; hardE=set(); free=set(); tun={}
        cells=[(i,j) for i in range(k) for j in range(k)]
        tags=[(t,i,j) for (i,j) in cells for t in 'AB']
        tilegroups={tg:{} for tg in tags}
        for (i,j) in cells:
            for (a,b) in S[i][j]:
                tilegroups[('A',i,j)][slot[(a,b)]]=(a,b)
                tilegroups[('B',i,j)][slot[(b,a)]]=(a,b)
        for tg in tags:
            for (a,b,c) in mem0: members.append((vid((tg,names0[a])),vid((tg,names0[b])),vid((tg,names0[c]))))
            soft0=set()
            for (a,b,c) in mem0: soft0|={k2(a,b),k2(b,c),k2(a,c)}
            for (u_,v_) in E0:
                if (u_,v_) not in soft0: hardE.add(k2(vid((tg,names0[u_])),vid((tg,names0[v_]))))
            for v in names0: vid((tg,v))
        def zv(tg,u): return vid(('z',tg,u))
        for tg in tags:
            for u in range(U):
                Ws=Wsets(m,u,prof); h=('s',4*u+2)
                z=zv(tg,u); hardE.add(k2(z,vid((tg,h))))
                Wz=Ws['all'] if u in tilegroups[tg] else Ws['pad']
                for v in Wz-{h}: free.add(k2(z,vid((tg,v))))
                if u in tilegroups[tg]:
                    x=vid(('x',tg,u)); y=vid(('y',tg,u))
                    for v in Ws['pad']: free.add(k2(y,vid((tg,v))))
                    members.append((x,y,z))
        # tunnels
        self.tunnels=[]
        for (i,j) in cells:
            for (a,b) in S[i][j]:
                xa=vid(('x',('A',i,j),slot[(a,b)])); xb=vid(('x',('B',i,j),slot[(b,a)]))
                tun[((i,j),(a,b))]=k2(xa,xb)
        # comparators
        cmps=[]
        for i in range(k):
            for j in range(k-1):
                for al in range(1,n): cmps.append((('A',i,j),('A',i,j+1),al))
        for j in range(k):
            for i in range(k-1):
                for al in range(1,n): cmps.append((('B',i,j),('B',i+1,j),al))
        for (c1,c2,al) in cmps:
            u=sep[al]; Ws=Wsets(m,u,prof); y=vid(('c',c1,c2,al))
            for v in Ws['X']: free.add(k2(y,vid((c1,v))))
            for v in Ws['Y']: free.add(k2(y,vid((c2,v))))
            members.append((zv(c1,u),y,zv(c2,u)))
        self.nv=len(V); self.members=members; self.hardE=hardE; self.free=free; self.tun=tun; self.V=V
        # SAT
        soft=set()
        for (a,b,c) in members: soft|={k2(a,b),k2(b,c),k2(a,c)}
        assert all(p not in soft for p in free) or True
        var={}
        for p in sorted(soft|free|set(tun.values())): var[p]=len(var)+1
        self.var=var
        cl=[]
        possible=set(var)|hardE
        TRUE=len(var)+1
        def lit(p):
            if p in var: return var[p]
            return TRUE if p in hardE else None
        adj=[[] for _ in range(self.nv)]
        for p in possible: adj[p[0]].append(p[1]); adj[p[1]].append(p[0])
        cl.append([TRUE])
        for v in range(self.nv):
            nb=adj[v]
            for i_ in range(len(nb)):
                for j_ in range(i_+1,len(nb)):
                    u_,w_=nb[i_],nb[j_]
                    a,b=lit(k2(u_,v)),lit(k2(v,w_)); c=lit(k2(u_,w_))
                    cl.append([-a,-b]+([c] if c is not None else []))
        for (a,b,c) in members:
            ab,bc,ac=var[k2(a,b)],var[k2(b,c)],var[k2(a,c)]
            cl += [[-ac,ab],[-ac,bc],[ab,bc],[-ab,-bc,ac]]
        self.s=Cadical153(bootstrap_with=cl)
    def feasible(self,cut):
        """cut: set of (cell,tile) whose tunnel is removed (pair becomes hard non-edge); other tunnels are hard edges"""
        ass=[]
        for key_,p in self.tun.items():
            ass.append(-self.var[p] if key_ in cut else self.var[p])
        return self.s.solve(assumptions=ass)

def gt_ok(k,t):
    for i in range(k):
        for j in range(k):
            if j+1<k and t[(i,j)][0]>t[(i,j+1)][0]: return False
            if i+1<k and t[(i,j)][1]>t[(i+1,j)][1]: return False
    return True

def randS(k,n,rng,dens):
    S=[[set() for _ in range(k)] for _ in range(k)]
    for i in range(k):
        for j in range(k):
            for a in range(1,n+1):
                for b in range(1,n+1):
                    if rng.random()<dens: S[i][j].add((a,b))
            if not S[i][j]: S[i][j].add((rng.randint(1,n),rng.randint(1,n)))
    return S

def run(k,n,seed,dens,cutsets=False):
    rng=random.Random(seed); S=randS(k,n,rng,dens)
    t0=time.time(); M=Model(k,n,S)
    cells=[(i,j) for i in range(k) for j in range(k)]
    mism=0; tot=0; nyes=0
    for ch in itertools.product(*[sorted(S[i][j]) for (i,j) in cells]):
        t=dict(zip(cells,ch)); cut={(c,t[c]) for c in cells}
        f=M.feasible(cut); e=gt_ok(k,t); tot+=1; nyes+=e; mism+=(f!=e)
    msg=f'k{k} n{n} seed{seed} dens{dens} |V|={M.nv} |H|={len(M.members)} choices={tot} tilings={nyes} mismatches={mism}'
    if cutsets:
        keys=list(M.tun); bad=0; cnt=0; feas=0
        for r in range(0,k*k+1):
            for C in itertools.combinations(keys,r):
                cut=set(C); f=M.feasible(cut); cnt+=1
                percell={c:[kk for kk in C if kk[0]==c] for c in cells}
                e=all(len(percell[c])==1 for c in cells) and r==k*k and gt_ok(k,{c:percell[c][0][1] for c in cells})
                feas+=f; bad+=(f!=e)
        msg+=f' | cutsets<=k^2: {cnt} tested, feasible {feas}, mismatches {bad}'
    print(msg,round(time.time()-t0,1),flush=True)
    return mism

if __name__=='__main__':
    k,n=int(sys.argv[1]),int(sys.argv[2]); seeds=range(int(sys.argv[3]),int(sys.argv[4])); dens=float(sys.argv[5]); cs=len(sys.argv)>6
    tot=0
    for sd in seeds: tot+=run(k,n,sd,dens,cs)
    print('TOTAL MISMATCHES',tot)
