import itertools, random
def partitions(n):
    def rec(i, blocks):
        if i == n:
            yield [list(b) for b in blocks]; return
        for b in blocks:
            b.append(i); yield from rec(i+1, blocks); b.pop()
        blocks.append([i]); yield from rec(i+1, blocks); blocks.pop()
    yield from rec(0, [])
def mk(n, edges):
    E = set(frozenset(e) for e in edges)
    return E
def cost(n, E, blocks):
    lab = {}
    for i,b in enumerate(blocks):
        for v in b: lab[v]=i
    c = 0
    for a,b in itertools.combinations(range(n),2):
        if (frozenset((a,b)) in E) != (lab[a]==lab[b]): c+=1
    return c
def allparts_cost(n,E):
    res=[]
    for p in partitions(n):
        res.append((cost(n,E,p),p))
    return res
def opt(n,E):
    return min(c for c,_ in allparts_cost(n,E))
def p3s(n,E):
    out=[]
    for y in range(n):
        for x,z in itertools.combinations([v for v in range(n) if v!=y],2):
            if frozenset((x,y)) in E and frozenset((y,z)) in E and frozenset((x,z)) not in E:
                out.append((x,y,z))
    return out
def pairs3(t):
    x,y,z=t
    return {frozenset((x,y)),frozenset((y,z)),frozenset((x,z))}
def rand_graph(n,p,rng):
    return set(frozenset(e) for e in itertools.combinations(range(n),2) if rng.random()<p)
def rand_packing(n,E,rng,maximal=True,keep=1.0):
    ts=p3s(n,E); rng.shuffle(ts); used=set(); H=[]
    for t in ts:
        if rng.random()>keep: continue
        P=pairs3(t)
        if not (P&used): H.append(t); used|=P
    if maximal:
        for t in ts:
            P=pairs3(t)
            if not (P&used): H.append(t); used|=P
    return H,used
