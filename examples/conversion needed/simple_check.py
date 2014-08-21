from sympy import symbols,sin,cos,pi
from ga import Ga
from printer import Eprint,xpdf,Format
from metric import linear_expand

def main():
    #Eprint()
    Format()

    o3d = Ga.preset('o3d')
    print 'I =',o3d.I
    x,y,z = o3d.coords

    (ex,ey,ez) = o3d.mv()

    u, v, w = o3d.mv('u v w', 'vector')

    print '# 3D Orthogonal Rectangular Coordinates'

    print 'u =',u
    print 'v =',v
    print 'w =',w

    uv = u^v
    print 'u^v =',uv
    uvw = u^v^w
    print 'u^v^w =',uvw

    print '(e_x+e_y)(e_y+e_z) =',(ex+ey)*(ey+ez)

    print r'\nabla =',o3d.grad
    print r'\dot{\nabla} =',o3d.grad

    f= o3d.mv('F', 'vector',f=True)
    print 'F =',f
    print r'\nabla F =',o3d.grad*f
    print r'F \dot{\nabla} =',f*o3d.rgrad
    print r'F \nabla =',f*o3d.grad
    print r'\dot{\nabla} F =',o3d.rgrad*f
    print r'F | \nabla =',f | o3d.grad

    G = y * z * ex + y ** 2 * ey + z ** 3 * ez

    print 'G =',G
    print r'\nabla G =', o3d.grad * G
    print r'\nabla ^ G =', o3d.grad ^ G
    print '#3D Orthogonal Spherical Coordinates'

    sp3d = Ga.preset('sph3d')
    r, th, phi = sp3d.coords
    er, eth, ephi = sp3d.mv()
    f, F = sp3d.mv('f F','scalar vector',f=True)

    print 'f =',f
    print r'\nabla =',sp3d.grad
    print r'\nabla f =',sp3d.grad*f
    print 'F =',F
    print r'\nabla | F =',sp3d.grad|F
    print r'-I\nabla ^ F =',-sp3d.I*(sp3d.grad^F)

    Lap = sp3d.grad|sp3d.grad

    print r'%\nabla^{2} =',Lap
    print r'%\nabla^{2} f =',Lap*f

    print '#Parabolic submanifold $[u,v,u^{2}+v^{2}]$ of o3d'

    sm_coords = (u,v) = symbols('u,v',real=True)

    Ga.set_simp(1,1,1)  # Additionally turn on combining fractions

    para = o3d.sm([u,v,u**2+v**2],sm_coords)

    print 'g_{parabolic} =', para.g

    (eu,ev) = para.mv()
    para_grad = para.grad

    print r'\nabla_{parabolic} =',para_grad

    f,F = para.mv('f F','scalar vector',f=True)

    print 'f = ', f
    print 'F = ', F

    print r'\nabla_{parabolic} f =',para_grad*f
    print r'\nabla_{parabolic} | F =',para_grad|F
    print r'\nabla_{parabolic} ^ F =',para_grad^F
    print r'\nabla_{parabolic}  F =',para_grad*F

    Ga.set_simp(1,1,0)  # Additionally turn off combining fractions

    sph = o3d.sm([sin(u)*cos(v),sin(u)*sin(v),cos(u)],sm_coords)

    print '#Spherical submanifold $[sin(u)*cos(v),sin(u)*sin(v),cos(u)]$ of o3d'

    (eu,ev) = sph.mv()
    sph_grad = sph.grad

    print 'g_{sphere} =', sph.g

    print r'\nabla_{spherical} =',sph_grad

    f,F = sph.mv('f F','scalar vector',f=True)

    print r'\nabla_{sphere} f =',sph_grad*f
    print r'\nabla_{sphere} | F =',sph_grad|F
    print r'\nabla_{sphere} ^ F =',sph_grad^F
    print r'\nabla_{sphere}  F =',sph_grad*F

    s = symbols('s',real=True)

    circ = sph.sm([pi/4,s],(s,))

    print '#Circular submanifold $[pi/4,s]$ of sphere (latitude 45 deg North)'

    es = circ.mv()
    circ_grad = circ.grad

    print 'g_{circ} =',circ.g

    print r'\nabla_{circle} =',circ_grad

    f,F = circ.mv('f F','scalar vector',f=True)

    print r'\nabla_{circle} f =',circ_grad*f
    print r'\nabla_{circle} F =',circ_grad*F

    print r'#Unit Great Circle Manifold at $\pi/4$ radians to Equator'

    #Calculate great circle at 45 degs to equator through x-axis

    a = pi/4
    u = ey ^ ez
    X = cos(s) * ex + sin(s) * ey
    R = cos(a/2) + u * sin(a/2)

    RXRrev = R*X*R.rev()
    print 'X(s) =',RXRrev

    M_gc = o3d.sm(RXRrev.list(),(s,))
    es = M_gc.mv()
    print 'g =',M_gc.g
    print r'\nabla_{gc} =',M_gc.grad
    xpdf()

    return

if __name__ == "__main__":
    main()
