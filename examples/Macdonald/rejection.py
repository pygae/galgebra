from GA import *

enhance_print()

(e1,e2,e3) = MV.setup('e_1 e_2 e_3','[1,1,1]')
(v1,v2,v3,u1,u2,u3) = symbols('v__1 v__2 v__3 u__1 u__2 u__3')
v = MV('v', 'vector')
u = MV('u', 'vector')
uv = u^v
print uv
uv_sq = (uv*uv).scalar()
print uv_sq
print expand(uv_sq +((u1*v2-u2*v1)**2+(u1*v3-u3*v1)**2+(u2*v3-u3*v2)**2))
