import sys, itertools
import numpy as np
sys.argv=['x','5']
from w11_maj5 import *
mem=[(1,0,3),(2,0,4),(2,1,5)]
hp=[(0,5),(1,4),(2,3),(3,4),(3,5),(4,5)]
hard={p:0 for p in hp}
R=rel5(mem,6,hard,())
names=['L','R','M','4','5']
print('size',len(R))
for r in R.tolist(): print(''.join(names[x] for x in r), end='  ')
print()
# defect-free part
ze=[r for r in R.tolist() if all(x<3 for x in r)]
print('ZE states', [''.join(names[x] for x in r) for r in ze])
one=[r for r in R.tolist() if sum(1 for x in r if x>=3)==1]
print('one defect', [''.join(names[x] for x in r) for r in one])
