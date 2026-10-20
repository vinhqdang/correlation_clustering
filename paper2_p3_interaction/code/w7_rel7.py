"""Independent re-verification of Theorem A on REAL graphs: every 5-vertex configuration of two members sharing exactly one vertex,
with all 4 cross pairs hard (edge / non-edge) -> ZE set Sigma (via cluster-graph test, w7_lib.Inst.sigma, no use of the ok3 formula).
Check: Sigma is closed under EVERY majority operation (all 729 uniform rules at once)."""
import itertools, numpy as np
from w7_lib import Inst, ST
def closed_all_majorities(vecs):
    vs=set(vecs)
    for u in vecs:
        for v in vecs:
            for w in vecs:
                opts=[]
                for x,y,z in zip(u,v,w):
                    if x==y or x==z: opts.append(x)
                    elif y==z: opts.append(y)
                    else: opts.append(ST)       # all-distinct: every output is realised by some rule
                for t in itertools.product(*opts):
                    if ''.join(t) not in vs: return False
    return True
n=0; bad=0; sizes={}
verts=[0,1,2,3,4]   # q=0, A-others 1,2 ; B-others 3,4
for oA in itertools.permutations([0,1,2]):
    for oB in itertools.permutations([0,3,4]):
        for eps in itertools.product((0,1),repeat=4):
            A=np.zeros((5,5),bool)
            for (a,b,c) in (oA,oB):
                A[a,b]=A[b,a]=A[b,c]=A[c,b]=True
            for (p,r),e in zip([(1,3),(1,4),(2,3),(2,4)],eps):
                A[p,r]=A[r,p]=bool(e)
            I=Inst(A,[oA,oB])
            if not I.valid(): continue
            sg=I.sigma(); n+=1
            sizes[len(sg)]=sizes.get(len(sg),0)+1
            if not closed_all_majorities(sg): bad+=1
print('configs (valid):',n,'|Sigma| histogram:',dict(sorted(sizes.items())),'not closed under all majorities:',bad)
# sanity: the same test on the soft triangle (K_4 example) must fail for M-majority
