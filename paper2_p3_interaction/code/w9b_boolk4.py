import pickle, itertools
from w9b_bool import bool_nu
stored=pickle.load(open('/tmp/w9b_k4.pkl','rb'))
for dom in [(0,1),(0,2),(1,2)]:
    rels=set()
    for fs in stored:
        sub=frozenset(tuple(0 if x==dom[0] else 1 for x in t) for t in fs if all(x in dom for x in t))
        if len(sub)>0: rels.add(sub)
    # remove full? keep
    print('domain',dom,'distinct boolean 4-ary relations',len(rels))
    for r in (3,4,5,6):
        res=bool_nu(list(rels),r)
        print('  NU arity',r,'->',res is not None,flush=True)
