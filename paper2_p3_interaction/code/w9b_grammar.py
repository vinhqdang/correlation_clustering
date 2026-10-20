"""Critical (irredundant) list structures on trees of blocks: value iteration on types (Y target, Z actual root set).
Blocks (c=2): K3 relations (3 members, all are ports), diamond relations (4 members, ports = members 0 and 3).  Links: any of the 40 binary relations (over-approximation: roles ignored).
Returns the maximum size of a critical structure after t rounds; unbounded growth <=> unbounded list cores on such trees."""
import itertools, pickle, sys
from w9b_poly import *
from w9b_nuP import lang_c2
D=(0,1,2)
SUBSETS=[frozenset(s) for k in (1,2,3) for s in itertools.combinations(D,k)]
FULL=frozenset(D)
def sset(x): return frozenset(x)

def build(blocks,links):
    return blocks,links
def run(blocks,links,rounds=12,verbose=True):
    # blocks: list of (rel (set of tuples), ports (indices))
    # links: list of relations (set of pairs) -- orientation: R(port_state, child_root_state); include inverses
    links=list({frozenset(R) for R in links}|{frozenset((b,a) for a,b in R) for R in links})
    pre={}
    def Pre(R,Z):
        return frozenset(a for (a,b) in R if b in Z)
    # M[(Y,Z)] = max size of Y-critical structure with root set exactly Z ; lists are structures of size 1: list L as a 'child' provider
    M={}
    # base: none. A list on a member is handled inside productions.
    # precompute production tables per block/root
    prods=[]
    for rel,ports in blocks:
        k=len(next(iter(rel)))
        rel=[tuple(t) for t in rel]
        for root in ports:
            others=[i for i in range(k) if i!=root]
            prods.append((rel,k,root,others,[i for i in others if i in ports]))
    hist=[]
    for rnd in range(rounds):
        newM=dict(M)
        for rel,k,root,others,oports in prods:
            for Qs in itertools.product(SUBSETS+[FULL],repeat=len(others)):
                # also root provider: list on root itself (Qroot) or none
                for Qr in SUBSETS+[FULL]:
                    Q={root:Qr}
                    for o,q in zip(others,Qs): Q[o]=q
                    sol=[t for t in rel if all(t[i] in Q[i] for i in range(k))]
                    Z=frozenset(t[root] for t in sol)
                    # targets Y (including empty meaning refutation) : Y ⊇ Z
                    for Y in SUBSETS+[frozenset()]:
                        if not Z<=Y or Y==FULL: continue
                        # necessity + providers
                        total=0; ok=True
                        for i in range(k):
                            if Q[i]==FULL: continue
                            # T_i: states s of member i so that output stays in Y when other Q fixed and x_i=s
                            Ti=frozenset(s for s in D if frozenset(t[root] for t in rel if t[i]==s and all(t[j] in Q[j] for j in range(k) if j!=i))<=Y)
                            if Ti==FULL: ok=False;break    # provider i unnecessary
                            # provider options for member i: list (size 1) realising Q[i] if Q[i] within Ti ; or child
                            best=None
                            if i!=root or True:
                                # a list gives Q exactly L; any L with L ⊆ Ti is allowed in principle, we need Q[i] itself to be list-realisable: Q[i] is a SUBSET (any nonempty subset is a list)
                                best=1 if Q[i]<=Ti else None
                            if i in oports or i==root:
                                for R in links:
                                    for (Yc,Zc),sz in list(M.items()):
                                        if Pre(R,Zc)!=Q[i]: continue
                                        Yi=frozenset(z for z in D if Pre(R,{z})<=Ti)
                                        if not (Zc<=Yi and Yi<=Yc): continue
                                        if best is None or sz>best: best=sz
                            if best is None: ok=False;break
                            total+=best
                        if not ok: continue
                        key=(Y,Z)
                        if Z==frozenset() and Y!=frozenset(): continue
                        if total>newM.get(key,0) and Z: newM[key]=total
                        if not Z:
                            newM[('top',)]=max(newM.get(('top',),0),total)
        M=newM
        top=M.get(('top',),0)
        hist.append((rnd,top,len(M),max([v for k,v in M.items() if k!=('top',)],default=0)))
        if verbose: print(hist[-1],flush=True)
    return M,hist
if __name__=='__main__':
    rels=lang_c2()
    blocks=[]
    for k,(a,R) in rels.items():
        if k[0]=='K3': blocks.append((R,[0,1,2]))
        elif k[0]=='dia': blocks.append((R,[0,3]))   # order of members in diamond Sigma: A,B,C,D -> ports A(0),D(3)
    # reduce: dedupe blocks
    seen={}
    for R,p in blocks: seen[(frozenset(R),tuple(p))]=1
    blocks=[(set(R),list(p)) for (R,p) in seen]
    links=[R for k,(a,R) in rels.items() if k[0]=='B' and len(R)]
    print('blocks',len(blocks),'links',len(links),flush=True)
    run(blocks,links,rounds=int(sys.argv[1]) if len(sys.argv)>1 else 8)
