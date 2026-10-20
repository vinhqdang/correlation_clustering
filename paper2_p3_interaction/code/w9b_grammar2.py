"""faster version of the critical-structure grammar on trees of blocks (see w9b_grammar.py)."""
import itertools, pickle, sys, time
from collections import defaultdict
from w9b_nuP import lang_c2
D=(0,1,2)
SUBSETS=[frozenset(s) for k in (1,2,3) for s in itertools.combinations(D,k)]
FULL=frozenset(D)
ALLQ=SUBSETS+[FULL]
def Pre(R,Z): return frozenset(a for (a,b) in R if b in Z)
def prepare(blocks):
    """for every (block,root): dict  Qtuple -> (Z, {Y: [T_i per member or None]})  computed once"""
    prep=[]
    for bi,(rel,ports) in enumerate(blocks):
        rel=[tuple(t) for t in rel]; k=len(rel[0])
        for root in ports:
            tab={}
            for Q in itertools.product(ALLQ,repeat=k):
                sol=[t for t in rel if all(t[i] in Q[i] for i in range(k))]
                if not sol: Z=frozenset()
                else: Z=frozenset(t[root] for t in sol)
                tab[Q]=Z
            prep.append((bi,root,k,rel,tab))
    return prep
def run(blocks,links,rounds,verbose=True):
    links=list({frozenset(R) for R in links}|{frozenset((b,a) for a,b in R) for R in links})
    prep=prepare(blocks)
    allY=SUBSETS[:-1]+[frozenset()]    # all subsets except FULL (7 incl. empty minus full: singles,pairs,empty)
    allY=[y for y in [frozenset()]+SUBSETS if y!=FULL]
    M={}
    top_hist=[]
    global PROV,TOPPROV
    PROV={}; TOPPROV=None
    # precompute Pre tables for links: for each link and each Z
    for rnd in range(rounds):
        # child-best cache for this round: key (Ti,Qi) -> size
        cache={}
        Mitems=list(M.items())
        def childbest(Ti,Qi):
            key=(Ti,Qi)
            if key in cache: return cache[key]
            best=None; bestp=None
            for R in links:
                Yi=frozenset(z for z in D if Pre(R,{z})<=Ti)
                for (Yc,Zc),sz in Mitems:
                    if Pre(R,Zc)==Qi and Zc<=Yi and Yi<=Yc:
                        if best is None or sz>best: best=sz; bestp=(R,(Yc,Zc))
            cache[key]=(best,bestp); return cache[key]
        newM=dict(M); top=0
        for bi,root,k,rel,tab in prep:
            ports=blocks[bi][1]
            for Q,Z in tab.items():
                for Y in allY:
                    if not Z<=Y: continue
                    if (not Z) and Y: continue          # empty output only for target empty
                    total=0; ok=True; provs=[]
                    for i in range(k):
                        if Q[i]==FULL: continue
                        # T_i
                        Ti=frozenset(s for s in D if frozenset(t[root] for t in rel if t[i]==s and all(t[j] in Q[j] for j in range(k) if j!=i))<=Y)
                        if Ti==FULL: ok=False;break
                        best=1 if Q[i]<=Ti else None; prov=('list',Q[i]) if best else None
                        if i in ports:
                            cb,cp=childbest(Ti,Q[i])
                            if cb is not None and (best is None or cb>best): best=cb; prov=('child',)+cp
                        if best is None: ok=False;break
                        total+=best; provs.append((i,prov))
                    if not ok or total==0: continue
                    if not Z:
                        if total>top: top=total; topprov=(bi,root,provs)
                    else:
                        key=(Y,Z)
                        if total>newM.get(key,0): newM[key]=total; PROV[key]=(bi,root,Q,provs)
        M=newM
        if top and (TOPPROV is None or top>TOPPROV[0]): TOPPROV=(top,topprov)
        mx=max(M.values(),default=0)
        top_hist.append((rnd,top,len(M),mx))
        if verbose: print(top_hist[-1],flush=True)
    return M,top_hist
if __name__=='__main__':
    which=sys.argv[1]; rounds=int(sys.argv[2])
    rels=lang_c2()
    blocks=[]; seen=set()
    for k,(a,R) in rels.items():
        if which in ('K3','both') and k[0]=='K3' and (frozenset(R),'K3') not in seen: seen.add((frozenset(R),'K3')); blocks.append((set(R),[0,1,2]))
        if which in ('dia','both') and k[0]=='dia' and (frozenset(R),'dia') not in seen: seen.add((frozenset(R),'dia')); blocks.append((set(R),[0,3]))
    links=[R for k,(a,R) in rels.items() if k[0]=='B' and len(R)]
    print(which,'blocks',len(blocks),'links',len(links),flush=True)
    t0=time.time()
    run(blocks,links,rounds)
    import pickle
    pickle.dump((blocks,PROV,TOPPROV),open('/tmp/w9b_gram_%s.pkl'%which,'wb'))
    print('time',round(time.time()-t0))
