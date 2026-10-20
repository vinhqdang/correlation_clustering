import itertools, io, contextlib
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from w7_tri import sigmas
from w7_lib import two_decomposable
from w7_rel7 import closed_all_majorities
nd=0; ncl=0; tot=0; ex=None; seen=set()
for key,sig in sigmas.items():
    fs=frozenset(sig)
    if fs in seen: continue
    seen.add(fs)
    tot+=1
    v=[''.join(s) for s in sig]
    if not two_decomposable(v): nd+=1; ex=ex or (key,v)
    if not closed_all_majorities(v): ncl+=1
print('distinct triangle Sigma sets:',tot,'non-2-decomposable:',nd,'not closed under all majorities:',ncl,'example non-closed/non2dec:',ex)
