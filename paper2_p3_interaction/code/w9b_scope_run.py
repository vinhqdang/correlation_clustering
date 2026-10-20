import pickle, time, itertools, sys
from multiprocessing import Pool
from w9b_scope import *
hardpairs = [(2, 4), (2, 5), (4, 5), (2, 3), (4, 0), (5, 1)]
def work(args):
    A,B=args
    out={}
    for C in role_orders((3,0,5)):
        for st in itertools.product((0,1,None),repeat=6):
            hard={frozenset(p):s for p,s in zip(hardpairs,st)}
            R=scope_relation([A,B,C],hard,6)
            if R not in out: out[R]=(A,B,C,st)
    return out
if __name__=='__main__':
    t=time.time()
    jobs=[(A,B) for A in role_orders((0,1,2)) for B in role_orders((1,3,4))]
    allr={}
    with Pool(4) as p:
        for o in p.imap_unordered(work,jobs):
            for R,v in o.items(): allr.setdefault(R,v)
    print('triangle scope relations',len(allr),time.time()-t)
    P={}
    cross=[(1,3),(1,4),(2,3),(2,4)]
    for A in role_orders((0,1,2)):
        for B in role_orders((0,3,4)):
            for st in itertools.product((0,1,None),repeat=4):
                hard={frozenset(p):s for p,s in zip(cross,st)}
                R=scope_relation([A,B],hard,5)
                P.setdefault(R,(A,B,st))
    print('pair scope relations',len(P))
    pickle.dump((P,allr),open('/tmp/w9b_scopes.pkl','wb'))
    from collections import Counter
    print('tri sizes',sorted(Counter(len(r) for r in allr).items()))
