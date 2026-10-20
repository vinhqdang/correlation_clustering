"""critical-structure grammar with explicit witness trees (K3 / diamond blocks, over-approximate links)."""
import itertools, pickle, sys, time
from w9b_nuP import lang_c2
D=(0,1,2)
SUBSETS=[frozenset(s) for k in (1,2,3) for s in itertools.combinations(D,k)]
FULL=frozenset(D); ALLQ=SUBSETS+[FULL]
def Pre(R,Z): return frozenset(a for (a,b) in R if b in Z)
def prepare(blocks):
    prep=[]
    for bi,(rel,ports) in enumerate(blocks):
        rel=[tuple(t) for t in rel]; k=len(rel[0])
        for root in ports:
            tab={}
            for Q in itertools.product(ALLQ,repeat=k):
                sol=[t for t in rel if all(t[i] in Q[i] for i in range(k))]
                tab[Q]=frozenset(t[root] for t in sol) if sol else frozenset()
            prep.append((bi,root,k,rel,tab))
    return prep
def run(blocks,links,rounds,verbose=True):
    links=list({frozenset(R) for R in links}|{frozenset((b,a) for a,b in R) for R in links})
    prep=prepare(blocks)
    allY=[y for y in [frozenset()]+SUBSETS if y!=FULL]
    M={}   # (Y,Z) -> (size,tree)
    tops=[]
    best_top=(0,None)
    for rnd in range(rounds):
        cache={}
        Mitems=list(M.items())
        def childbest(Ti,Qi):
            key=(Ti,Qi)
            if key in cache: return cache[key]
            best=None
            for R in links:
                Yi=frozenset(z for z in D if Pre(R,{z})<=Ti)
                for (Yc,Zc),(sz,tr) in Mitems:
                    if Pre(R,Zc)==Qi and Zc<=Yi and Yi<=Yc:
                        if best is None or sz>best[0]: best=(sz,('child',R,tr))
            cache[key]=best; return best
        newM=dict(M); top=(0,None)
        for bi,root,k,rel,tab in prep:
            ports=blocks[bi][1]
            for Q,Z in tab.items():
                for Y in allY:
                    if not Z<=Y: continue
                    if (not Z) and Y: continue
                    total=0; ok=True; provs=[]
                    for i in range(k):
                        if Q[i]==FULL: continue
                        Ti=frozenset(s for s in D if frozenset(t[root] for t in rel if t[i]==s and all(t[j] in Q[j] for j in range(k) if j!=i))<=Y)
                        if Ti==FULL: ok=False;break
                        best=(1,('list',Q[i])) if Q[i]<=Ti else None
                        if i in ports:
                            cb=childbest(Ti,Q[i])
                            if cb is not None and (best is None or cb[0]>best[0]): best=cb
                        if best is None: ok=False;break
                        total+=best[0]; provs.append((i,best[1]))
                    if not ok or total==0: continue
                    tree=('blk',bi,root,Q,tuple(provs))
                    if not Z:
                        if total>top[0]: top=(total,tree)
                    else:
                        key=(Y,Z)
                        if total>newM.get(key,(0,None))[0]: newM[key]=(total,tree)
        M=newM
        if top[0]>best_top[0]: best_top=top
        tops.append(top[0])
        if verbose: print(rnd,'top',top[0],'types',len(M),flush=True)
    return M,best_top
def count_lists(tree):
    if tree[0]=='blk': return sum(count_lists(p) for i,p in tree[4])
    if tree[0]=='list': return 1
    return count_lists(tree[2])
def show(tree,blocks,ind=0):
    pad='  '*ind
    if tree[0]=='blk':
        print(pad+'block',tree[1],'root',tree[2],'Q',[sorted(q) for q in tree[3]])
        for i,p in tree[4]: print(pad+' member',i,':'); show(p,blocks,ind+2)
    elif tree[0]=='list': print(pad+'list',sorted(tree[1]))
    else:
        print(pad+'link',sorted(tree[1])); show(tree[2],blocks,ind+1)
if __name__=='__main__':
    which=sys.argv[1]; rounds=int(sys.argv[2])
    rels=lang_c2(); blocks=[]; seen=set()
    for k,(a,R) in rels.items():
        if which in ('K3','both') and k[0]=='K3' and (frozenset(R),'K3') not in seen: seen.add((frozenset(R),'K3')); blocks.append((set(R),[0,1,2]))
        if which in ('dia','both') and k[0]=='dia' and (frozenset(R),'dia') not in seen: seen.add((frozenset(R),'dia')); blocks.append((set(R),[0,3]))
    links=[R for k,(a,R) in rels.items() if k[0]=='B' and len(R)]
    M,bt=run(blocks,links,rounds)
    print('best top',bt[0]); show(bt[1],blocks)
    pickle.dump((blocks,bt),open('/tmp/w9b_gram3_%s.pkl'%which,'wb'))
