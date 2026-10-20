"""Linked-pair path hitting (LPPH) and the reduction Grid Tiling with <= -> LPPH (own proof in notes/wave7_special.md).
LPPH: a forest (here: disjoint paths), a family of 'obstruction' paths, a set of links (pairs of nodes); an operation is a single node or a link;
a node set is hit by an operation if the node belongs to it; choose <= L operations so that every obstruction path contains a hit node.
This script builds the instance for k x k Grid Tiling with <= and checks equivalence with brute force on random tiny instances."""
import itertools, random, sys
def gt_solve(k,n,S):
    cells=[(i,j) for i in range(k) for j in range(k)]
    lists=[sorted(S[c]) for c in cells]
    for choice in itertools.product(*lists):
        t=dict(zip(cells,choice)); ok=True
        for i in range(k):
            for j in range(k):
                if j+1<k and not t[(i,j)][0]<=t[(i,j+1)][0]: ok=False;break
                if i+1<k and not t[(i,j)][1]<=t[(i+1,j)][1]: ok=False;break
            if not ok: break
        if ok: return True
    return False
def build(k,n,S):
    nid=0; arms={}     # ('A',i,j)-> list of (tile,node) in position order
    node_of={}
    for i in range(k):
        for j in range(k):
            for typ,key in (('A',lambda e:(e[0],e[1])),('B',lambda e:(e[1],e[0]))):
                tiles=sorted(S[(i,j)],key=key)
                arms[(typ,i,j)]=[]
                for e in tiles:
                    node_of[(typ,i,j,e)]=nid; arms[(typ,i,j)].append((e,nid)); nid+=1
    paths=[]
    for key,lst in arms.items(): paths.append(frozenset(v for _,v in lst))     # obstructions: arms
    # row constraints on A arms: a_ij <= a_{i,j+1}
    for i in range(k):
        for j in range(k-1):
            for alpha in range(1,n):
                suf=[v for e,v in arms[('A',i,j+1)] if e[0]>=alpha+1]
                pre=[v for e,v in arms[('A',i,j)] if e[0]<=alpha]
                paths.append(frozenset(suf+pre))
    # column constraints on B arms: b_ij <= b_{i+1,j}
    for j in range(k):
        for i in range(k-1):
            for beta in range(1,n):
                suf=[v for e,v in arms[('B',i+1,j)] if e[1]>=beta+1]
                pre=[v for e,v in arms[('B',i,j)] if e[1]<=beta]
                paths.append(frozenset(suf+pre))
    links=[(node_of[('A',i,j,e)],node_of[('B',i,j,e)]) for i in range(k) for j in range(k) for e in S[(i,j)]]
    return nid,paths,links
def lpph_solve(nid,paths,links,L):
    ops=[frozenset([v]) for v in range(nid)]+[frozenset(l) for l in links]
    byn={}
    for o in ops:
        for v in o: byn.setdefault(v,[]).append(o)
    def rec(hit,left):
        for P in paths:
            if not (P&hit): break
        else: return True
        if left==0: return False
        if not P: return False
        cand=set()
        for v in P:
            for o in byn[v]: cand.add(o)
        for o in cand:
            if rec(hit|o,left-1): return True
        return False
    return rec(frozenset(),L)
if __name__=='__main__':
    rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); k=2; n=3
    stats={True:0,False:0}; mism=0
    for _ in range(N):
        S={}
        for i in range(k):
            for j in range(k):
                allt=[(a,b) for a in range(1,n+1) for b in range(1,n+1)]
                S[(i,j)]=set(rng.sample(allt,rng.randint(2,4)))
        gt=gt_solve(k,n,S)
        nid,paths,links=build(k,n,S)
        ph=lpph_solve(nid,paths,links,k*k)
        stats[gt]+=1
        if gt!=ph: mism+=1; print('MISMATCH',S)
    print('instances: YES %d NO %d; mismatches %d'%(stats[True],stats[False],mism))
