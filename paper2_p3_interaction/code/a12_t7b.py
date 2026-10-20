import itertools
from a12_lib import *
from a12_t2b import relc, ST5
cnt=0; ex=[]
for qa in range(3):
  for qb in range(3):
    for b in itertools.product((0,1),repeat=4):
      R=relc(qa,qb,b,ST5)
      for sb in ST5:
        A={sa for (sa,t) in R if t==sb}
        if A and 'O' not in A and 'L' in A and 'R' in A: cnt+=1; ex.append((qa,qb,b,sb,sorted(A)))
print('binary slices excluding O but allowing L and R:',cnt, ex[:3])
