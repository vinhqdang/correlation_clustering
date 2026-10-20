"""Table L(S): for a member A=(a,b,c) (centre b) and an outside vertex t with N(t) cap {a,b,c} = S (all three pairs {x,t} hard and unedited by the member),
the partitions (states) of {a,b,c} compatible with a cluster graph on {a,b,c,t}.  States: L={a}{bc}, R={ab}{c}, M={abc}, O={a}{b}{c}, X={ac}{b}.
Also checks that no list contains O without also containing L or R... precisely: every list not containing O is a subset of one of {L},{R},{X},{M}."""
import itertools
def cluster_graph(adj,verts):
    for x,y,z in itertools.combinations(verts,3):
        e=[adj[frozenset((x,y))],adj[frozenset((y,z))],adj[frozenset((x,z))]]
        if sum(e)==2: return False
    return True
parts={'O':[{'a'},{'b'},{'c'}],'L':[{'a'},{'b','c'}],'R':[{'a','b'},{'c'}],'M':[{'a','b','c'}],'X':[{'a','c'},{'b'}]}
def inpart(P,x,y): return any(x in B and y in B for B in P)
table={}
for S in [frozenset(c) for k in range(4) for c in itertools.combinations('abc',k)]:
    ok=[]
    for name,P in parts.items():
        adj={}
        for x,y in itertools.combinations('abc',2): adj[frozenset((x,y))]=inpart(P,x,y)
        for x in 'abc': adj[frozenset((x,'t'))]=(x in S)
        if cluster_graph(adj,list('abct')): ok.append(name)
    table[''.join(sorted(S))]=''.join(sorted(ok))
print('L(S):',table)
ze={k:''.join(c for c in v if c in 'LRM') for k,v in table.items()}
print('L(S) on zero-excess states:',ze)
assert all(len(v)<=1 for v in ze.values()) or all(len(v)<=1 for k,v in ze.items() if k!='')
lists=set(table.values())
noO=[l for l in lists if 'O' not in l]
print('lists not containing O:',noO)
assert all(not ('L' in l and 'R' in l) for l in noO)
print('OK: a hard vertex cannot exclude O unless it also excludes L or R')
# views of the five states at the three vertices (y_ab,y_bc,y_ac)
ST={'L':(0,1,0),'R':(1,0,0),'M':(1,1,1),'O':(0,0,0),'X':(0,0,1)}
for s,(ab,bc,ac) in ST.items(): print(s,'view at a (ab,ac)=',(ab,ac),' at b (ab,bc)=',(ab,bc),' at c (bc,ac)=',(bc,ac))
