import numpy as np
from w7_lib import *
E=[(0,3),(0,5),(1,2),(1,4),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)]
A=np.zeros((6,6),bool)
for u,v in E: A[u,v]=A[v,u]=True
H=[(0,5,4),(1,2,5),(0,3,2),(1,4,3)]
I=Inst(A,H); print(I.valid(), I.cmax(), I.sigma(), I.soft_triangles())
s=I.sigma(); print(closed_mmaj_fast(s), two_decomposable(s))
