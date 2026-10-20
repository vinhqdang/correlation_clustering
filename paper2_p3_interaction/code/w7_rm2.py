"""Abstract model RM2 (tunnels hit conflict endpoints, breaks hit path nodes): test the packing-endpoint hypothesis.
cost(F)=|F|+nu(Omega minus conflicts with an endpoint released by F), nu = max number of node-disjoint conflict paths in a forest (= min breaks)."""
import random, itertools, sys
def tree_paths(parent,s,t):
    # path nodes between s and t in rooted forest given parent array (root parent=-1)
    anc=[];x=s
    while x!=-1: anc.append(x); x=parent[x]
    seen={a:i for i,a in enumerate(anc)}
    y=t;down=[]
    while y not in seen:
        down.append(y); y=parent[y]
    lca=y
    return set(anc[:seen[lca]+1])|set(down), lca
def depth(parent,x):
    d=0
    while parent[x]!=-1: x=parent[x]; d+=1
    return d
def nu(parent,confs):
    # greedy by deepest lca: returns max packing (list of conflicts) of node-disjoint paths
    items=[]
    for (s,t) in confs:
        P,l=tree_paths(parent,s,t); items.append((depth(parent,l),P,(s,t)))
    items.sort(key=lambda x:-x[0])
    used=set(); pack=[]
    for d,P,c in items:
        if not (P&used): pack.append(c); used|=P
    return pack
def solve(parent,confs,tunnels,restrict=None):
    best=10**9
    T=tunnels if restrict is None else [t for t in tunnels if t in restrict]
    for r in range(len(T)+1):
        for F in itertools.combinations(T,r):
            rel=set(v for t in F for v in t)
            res=[c for c in confs if c[0] not in rel and c[1] not in rel]
            best=min(best,len(F)+len(nu(parent,res)))
    return best
def main():
    rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); fails=0; tot=0; gain=0
    for _ in range(N):
        n=rng.randint(8,16)
        parent=[-1]+[rng.randrange(max(0,i-rng.choice([1,2,3])),i) for i in range(1,n)]
        k=rng.randint(4,9)
        confs=[]
        while len(confs)<k:
            s,t=rng.sample(range(n),2); confs.append((s,t))
        ends=sorted({v for c in confs for v in c})
        # tunnels: random matching-ish graph with fan-out<=2 on endpoint nodes, plus some on non-endpoints
        tun=set()
        for _ in range(rng.randint(3,9)):
            a,b=rng.sample(range(n),2)
            if sum(1 for t in tun if a in t)<2 and sum(1 for t in tun if b in t)<2: tun.add(tuple(sorted((a,b))))
        tun=sorted(tun)
        full=solve(parent,confs,tun)
        pk=nu(parent,confs); pe={v for c in pk for v in c}
        restr={t for t in tun if t[0] in pe or t[1] in pe}
        r1=solve(parent,confs,tun,restr)
        tot+=1; gain+= full<len(pk)
        if r1!=full: fails+=1; ex=(parent,confs,tun,full,r1)
    print('instances',tot,'with tunnel gain',gain,'restricted-to-packing-endpoints != optimum:',fails)
    if fails: print(ex)
if __name__=='__main__':
    main()
