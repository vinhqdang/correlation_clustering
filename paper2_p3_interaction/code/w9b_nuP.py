"""existence of a conservative NU of arity r with property P: f(t) occurs at least twice in t (whenever some value does)."""
import pickle, itertools, sys
from collections import Counter
from w9b_poly import *
from w9b_scope import *
def propP(ps):
    r=ps.r
    for i,t in enumerate(ps.tuples):
        c=Counter(t.tolist() if hasattr(t,'tolist') else t)
        rep=[v for v in D if c[v]>=2]
        if rep:
            for v in D:
                if v not in rep: ps.S.add_clause([-ps.var(i,v)])
def lang_c2():
    hardpairs=[(2,4),(2,5),(4,5),(2,3),(4,0),(5,1)]
    K3=set()
    for A in role_orders((0,1,2)):
        for B in role_orders((1,3,4)):
            for C in role_orders((3,0,5)):
                for st in itertools.product((0,1),repeat=6):
                    hard={frozenset(p):s for p,s in zip(hardpairs,st)}
                    R=scope_relation([A,B,C],hard,6)
                    if len(R): K3.add(R)
    dia=pickle.load(open('/tmp/w9b_diamond_rels.pkl','rb'))
    rels=dict(lists())
    for k,v in language('sf').items():
        if v[0]==2 and len(v[1]): rels[k]=v
    for i,R in enumerate(K3): rels[('K3',i)]=(3,frozenset(R))
    for i,R in enumerate(dia): rels[('dia',i)]=(4,frozenset(R))
    return rels
if __name__=='__main__':
    which=sys.argv[1]
    if which=='c2': rels=lang_c2()
    elif which=='k4':
        st=pickle.load(open('/tmp/w9b_k4.pkl','rb')); rels=dict(lists())
        for i,fs in enumerate(st): rels[('K4',i)]=(4,frozenset(fs))
    elif which=='c2k4':
        rels=lang_c2(); st=pickle.load(open('/tmp/w9b_k4.pkl','rb'))
        for i,fs in enumerate(st): rels[('K4',i)]=(4,frozenset(fs))
    print(which,'language size',len(rels),flush=True)
    for r in (4,5,6):
        ps=PolySearch(r,rels,nu=True,conservative=True,extra=propP)
        res=ps.solve()
        print('conservative NU arity',r,'with property P ->',res is not None,flush=True)
