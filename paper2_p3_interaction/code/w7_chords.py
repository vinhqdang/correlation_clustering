import itertools, sys
from w7_dpy import opt
m=int(sys.argv[1]); n=2*m+1
base=[(i,i+1) for i in range(2*m)]
cent=[2*i+1 for i in range(m)]
chs=[(a,b) for a,b in itertools.combinations(cent,2) if b-a>=4]
print('m',m,'centres',cent,'chords',len(chs))
from collections import Counter
cnt=Counter()
for k in (1,2,3):
    for C in itertools.combinations(chs,k):
        r,sols=opt(n,base+list(C),20000)
        ex=r-m
        # which chords are deleted in ALL optima
        dele=[c for c in C if all(s[c[0]]!=s[c[1]] for s in sols)]
        some=[c for c in C if any(s[c[0]]!=s[c[1]] for s in sols)]
        cnt[(k,ex,len(dele))]+=1
for key in sorted(cnt): print('chords',key[0],'excess',key[1],'chords deleted in every optimum',key[2],':',cnt[key],'instances')
