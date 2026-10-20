"""Audit 14: independent check of Lemma 1 (chain) and Lemma 2 (profile table) for OH_m.
Only the *instance data* (members, hard edges) is taken from w8c_oh.oh_instance; everything else (validity, ZE enumeration, profiles) is written fresh."""
import sys, os, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_oh import oh_instance
from pysat.solvers import Cadical153

def k2(a,b): return (a,b) if a<b else (b,a)

def check_valid(n,E,members):
    soft=set(); ok=True
    for (a,b,c) in members:
        ok &= (k2(a,b) in E and k2(b,c) in E and k2(a,c) not in E)
        ps=[k2(a,b),k2(b,c),k2(a,c)]
        for p in ps:
            if p in soft: ok=False
            soft.add(p)
    return ok, soft

def enum_ze(n,E,members,limit=10**6):
    ok,soft=check_valid(n,E,members); assert ok
    possible=set(soft)|set(E)
    var={}; 
    for p in sorted(possible): var[p]=len(var)+1
    cl=[]
    for p in E:
        if p not in soft: cl.append([var[p]])      # hard edge: together
    adj=[[] for _ in range(n)]
    for (u,v) in possible: adj[u].append(v); adj[v].append(u)
    for v in range(n):
        nb=adj[v]
        for i in range(len(nb)):
            for j in range(i+1,len(nb)):
                u,w=nb[i],nb[j]
                a,b=var[k2(u,v)],var[k2(v,w)]
                if k2(u,w) in var: cl.append([-a,-b,var[k2(u,w)]])
                else: cl.append([-a,-b])
    for (a,b,c) in members:
        ab,bc,ac=var[k2(a,b)],var[k2(b,c)],var[k2(a,c)]
        # allowed togetherness patterns (ab,bc,ac): L=(0,1,0) R=(1,0,0) M=(1,1,1)
        cl += [[ab,bc,ac,],]  # placeholder replaced below
        cl.pop()
        cl += [[-ac,ab],[-ac,bc],[ab,bc],[-ab,-bc,ac]]
    s=Cadical153(bootstrap_with=cl); sols=[]
    while s.solve():
        mod=set(l for l in s.get_model() if l>0)
        parent=list(range(n))
        def f(x):
            while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
            return x
        for p,vv in var.items():
            if vv in mod: parent[f(p[0])]=f(p[1])
        sols.append([f(v) for v in range(n)])
        s.add_clause([(-v if v in mod else v) for v in var.values()])
        if len(sols)>=limit: break
    s.delete(); return sols

def spine_word(m,idx,lab):
    out=[]
    for j in range(m):
        a,b,c=idx[('s',2*j)],idx[('s',2*j+1)],idx[('s',2*j+2)]
        t=(lab[a]==lab[b],lab[b]==lab[c],lab[a]==lab[c])
        out.append({(True,False,False):'R',(False,True,False):'L',(True,True,True):'M'}.get(t,'?'))
    return ''.join(out)

if __name__=='__main__':
    ms=[int(x) for x in sys.argv[1:]] or list(range(7,41))
    bad=0
    for m in ms:
        t0=time.time()
        n,E,members,names,idx=oh_instance(m)
        sols=enum_ze(n,E,members)
        words=sorted(spine_word(m,idx,l) for l in sols)
        exp=sorted('R'*a+'M'+'L'*(m-a-1) for a in range(m))
        distinct=len({tuple(sorted(set(l))) and tuple(l) for l in sols})
        # canonical partition check
        parts={frozenset(frozenset(v for v in range(n) if l[v]==c) for c in set(l)) for l in sols}
        good=(words==exp and len(parts)==len(sols)==m)
        bad+= (not good)
        print('m',m,'n',n,'|H|',len(members),'#ZE',len(sols),'one-hot words',words==exp,'distinct partitions',len(parts),'OK' if good else 'FAIL',round(time.time()-t0,1),flush=True)
    print('failures',bad)
