import numpy,math
from enthought.mayavi.mlab import *

def timeline():
    def f(tau):
        return(tau**2/2.0)
    def g(tau):
        return(0.5*(tau*math.sqrt(1+tau**2)+math.asinh(tau)))
    def df(tau):
        return(tau)
    def dg(tau):
        return(math.sqrt(1+tau**2))
    return(f,g,df,dg)

def plot_tau_r(f,g,df,dg,tau0,th0,ntau,nr):
    deg2rad = numpy.pi/180.0
    th0_rad = deg2rad*th0
    dtau = tau0/float(ntau)
    tau = numpy.arange(0.0,tau0+dtau,dtau)
    th  = numpy.arange(0.0,th0_rad,th0_rad/float(nr-1))
    fv = numpy.vectorize(f)
    gv = numpy.vectorize(g)
    x = fv(tau)
    y = numpy.zeros(len(x))
    z = gv(tau)
    lines = [plot3d(x,y,z,line_width=2.0,color=(0.0,0.0,0.0),tube_radius=None)]
    x1 = f(tau0)*numpy.cos(th)
    y1 = f(tau0)*numpy.sin(th)
    z1 = g(tau0)*numpy.array(len(x1)*[1.0])
    lines.append(plot3d(x1,y1,z1,line_width=2.0,color=(0.0,0.0,0.0),tube_radius=None))
    x_pt = f(tau0)*numpy.cos(th0_rad)
    y_pt = f(tau0)*numpy.sin(th0_rad)
    z_pt = g(tau0)
    e_th = quiver3d([x_pt],[y_pt],[z_pt],\
                    [-numpy.sin(th0_rad)],[numpy.cos(th0_rad)],[0.0],\
                    line_width=2.0,color=(0.0,1.0,0.0))
    e_tau = quiver3d([x_pt],[y_pt],[z_pt],\
                    [df(tau0)*numpy.cos(th0_rad)],[df(tau0)*numpy.sin(th0_rad)],[dg(tau0)],\
                    line_width=2.0,color=(0.0,0.0,1.0))
    lines.append(e_th)
    lines.append(e_tau)
    return(lines)

def coords(x0,y0,z0):
    xyz  = quiver3d([x0,x0,x0],[y0,y0,y0],[z0,z0,z0],\
                    [5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0],\
                    line_width=2.0,color=(1.0,0.0,0.0))
    return(xyz)

def surf_rev(f,g,tau,nth):
    pi = numpy.pi
    cos = numpy.cos
    sin = numpy.sin
    dth = pi/float(nth)
    th = numpy.arange(0.0,2*pi+1.5*dth,dth)
    (t,th) = numpy.meshgrid(tau,th)
    fv = numpy.vectorize(f)
    gv = numpy.vectorize(g)
    x = cos(th)*fv(t)
    y = sin(th)*fv(t)
    z = gv(t)
    return mesh(x, y, z,colormap='gist_earth')

tau = numpy.arange(0.01,2.015,0.01)

(f,g,df,dg) = timeline()

surf_rev(f,g,tau,250)
plot_tau_r(f,g,df,dg,1.5,45.0,50,50)
coords(0.5,-1.0,0.0)
#CR = raw_input()
