
# Initialize gnuplot
import Gnuplot,numpy,math

def wait(str=None, prompt='Press return to show results...\n'):
    if str is not None:
        print str
    raw_input(prompt)

zero = numpy.array([0.0,0.0,0.0])
ex = numpy.array([1.0,0.0,0.0])
ey = numpy.array([0.0,1.0,0.0])
ez = numpy.array([0.0,0.0,1.0])

def triad(th,phi):
    th_rd = math.pi*th/180.0
    phi_rd = math.pi*phi/180.0
    sth = math.sin(th_rd)
    cth = math.cos(th_rd)
    sphi = math.sin(phi_rd)
    cphi = math.cos(phi_rd)
    u = sth*ez+cth*(cphi*ex+sphi*ey)
    n = cth*ez-sth*(cphi*ex+sphi*ey) 
    t = sphi*ey-cphi*ex
    return(u,n,t)   

def seg(v1,v2):
    v12 = numpy.array([v1[0],v1[1],v1[2],v2[0],v2[1],v2[2]])
    return(v12)

def vec(v1,v2):
    v12 = numpy.array([v1[0],v1[1],v1[2],v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]])
    return(v12)

def drawarc(center,start,finish,rad=1.0,arrow=4.0,npts=20,thang=15.0):
    r1 = start-center
    r1sq = numpy.dot(r1,r1)
    r1mag = math.sqrt(r1sq)
    r2 = finish-center
    r2sq = numpy.dot(r2,r2)
    r2mag = math.sqrt(r2sq)
    u1 = r1/r1mag
    u2 = r2/r2mag
    costh = numpy.dot(u1,u2)
    sinth = math.sqrt(costh)
    u2perp = u2-costh*u1
    u2 = u2perp/math.sqrt(numpy.dot(u2perp,u2perp))
    th = math.acos(costh)
    dth = th/float(npts-1)
    rpts = []
    for ith in range(npts):
        th = float(ith)*dth
        cth = math.cos(th)
        sth = math.sin(th)
        r = center+rad*r1mag*(cth*u1+sth*u2)
        rpts.append(r)
        
    #Draw Arrow
    v1 = rpts[-2]
    v2 = rpts[-1]
    delv = v2-v1
    delv_mag = math.sqrt(numpy.dot(delv,delv))
    delu = delv/delv_mag
    cth = numpy.dot(u1,delu)
    sth = numpy.dot(u2,delu)
    delu_perp = -sth*u1+cth*u2
    delv_perp = delv_mag*math.sin(math.pi*thang/180.0)*delu_perp
    arrow1 = 1.0-arrow
    rpts.append(arrow1*v2+arrow*(v1+delv_perp))
    rpts.append(v2)
    rpts.append(arrow1*v2+arrow*(v1-delv_perp))
    arc = Gnuplot.Data(rpts,with_='lines ls 255')
    return(arc)

gp = Gnuplot.Gnuplot(persist = 1)
gp('set notics')
gp('set nolabel')
gp('set nogrid')
gp('set noborder')
gp('set notitle')
gp('set style line 255 lt rgb "black"')
gp('set view 65,100,0.5,1.5')
# Generate the vector field data

(u,n,t) = triad(45.0,45.0)

pt1 = 0.5*n+0.5*t+u
pt2 = 0.5*n-0.5*t+u
pt3 = -0.5*n-0.5*t+u
pt4 = -0.5*n+0.5*t+u
pts = [pt1,pt2,pt3,pt4,pt1]

Ex = vec(zero,ex)
Ey = vec(zero,ey)
Ez = vec(zero,ez)
U  = vec(zero,u)
x = -0.45*n+u
X  = vec(zero,x)
Xpar = vec(u,x)
delta = math.pi*60.0/180.0
xrot = 0.45*(-math.cos(delta)*n+math.sin(delta)*t)+u
Xrot = vec(zero,xrot)
Xrotpar = vec(u,xrot)
arc = drawarc(u,x,xrot,0.75)

v = [Ex,Ey,Ez,U,X,Xpar,Xrot,Xrotpar]

gv = Gnuplot.Data(v,with_ = 'vectors ls 255')
gpts = Gnuplot.Data(pts,with_ = 'lines ls 255')
# Plot everything: note that the with argument is 'vectors'
plot = gp.splot(gpts,gv,arc)
wait('Testing')
gp('set terminal postscript solid lw 0.5')
gp.hardcopy('rotor.eps', mode='eps')



