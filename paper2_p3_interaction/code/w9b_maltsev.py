"""existence of a (conservative) Maltsev polymorphism p(x,y,y)=x=p(y,y,x) for the languages."""
import itertools, pickle, sys
from w9b_poly import *
from w9b_nuP import lang_c2
def maltsev(ps):
    for i,t in enumerate(ps.tuples):
        x,y,z=int(t[0]),int(t[1]),int(t[2])
        if y==z: ps.S.add_clause([ps.var(i,x)])
        if x==y: ps.S.add_clause([ps.var(i,z)])
def test(name,rels,cons):
    ps=PolySearch(3,rels,nu=False,idem=True,conservative=cons,extra=maltsev)
    res=ps.solve()
    print(name,'(conservative)' if cons else '(idempotent)','Maltsev:',res is not None,flush=True)
test('sf',language('sf'),True)
test('c2 (K3,diamond,links,lists)',lang_c2(),True)
st=pickle.load(open('/tmp/w9b_k4.pkl','rb')); r=dict(lists())
for i,fs in enumerate(st): r[('K4',i)]=(4,frozenset(fs))
test('K4',r,True)
test('scope (8188 rel.)',None,True) if False else None
