"""Per-coordinate majority rules for the real edge relations R_conj (conjunction of the 4 triples at a shared vertex).
rule = 6 values (order of PERMS). Enumerate pairs (rA,rB) = 3^12 vectors, vectorised."""
import itertools, io, contextlib, pickle
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from w7_hub3 import cfgs, S
perms=list(itertools.permutations(S)); pidx={p:i for i,p in enumerate(perms)}
sidx={s:i for i,s in enumerate(S)}
allv=np.array(list(itertools.product(range(3),repeat=12)),dtype=np.int8)  # cols 0-5 rA, 6-11 rB
def compat_mask(R):
    R=sorted(R); Rs=set(R)
    good=np.ones(len(allv),dtype=bool)
    relmask=np.zeros((3,3),bool)
    for s,t in R: relmask[sidx[s],sidx[t]]=True
    for u,v,w in itertools.product(R,repeat=3):
        def res(col, off):
            a,b,c=u[col],v[col],w[col]
            if a==b or a==c: return np.full(len(allv),sidx[a],dtype=np.int8)
            if b==c: return np.full(len(allv),sidx[b],dtype=np.int8)
            return allv[:,off+pidx[(a,b,c)]]
        x=res(0,0); y=res(1,6)
        good&=relmask[x,y]
    return good
# collect unique relations by role pair
uniq={}
for (qm,qn,eps),(R,LP,LC,E) in cfgs.items():
    uniq.setdefault((qm,qn),set()).add(frozenset(R))
allowed={}  # (qm,qn,R) -> (set of rA, set of rB) as tuples
res={}
for (qm,qn),rels in uniq.items():
    for R in rels:
        if len(R)==9: continue
        g=compat_mask(R)
        sub=allv[g]
        pa={tuple(x[:6]) for x in sub}; pb={tuple(x[6:]) for x in sub}
        prod=len(sub)==len(pa)*len(pb)
        res[(qm,qn,R)]=(pa,pb,prod,len(sub))
pickle.dump(res,open('/tmp/w7_rel5.pkl','wb'))
print('relations:',len(res),'all product:',all(v[2] for v in res.values()),'any empty:',[k[:2] for k,v in res.items() if v[3]==0][:5])
# universal rule for role q (as parent-side qm=q or child-side qn=q)
U={}
for q in 'abc':
    sets=[]
    for (qm,qn,R),(pa,pb,prod,cnt) in res.items():
        if qm==q: sets.append(pa)
        if qn==q: sets.append(pb)
    U[q]=set.intersection(*sets) if sets else set()
    print('role',q,'universal rules:',len(U[q]))
print('abc common:',len(U['a']&U['b']&U['c']),' ab:',len(U['a']&U['b']),' ac:',len(U['a']&U['c']),' bc:',len(U['b']&U['c']))
