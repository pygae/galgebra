import numpy,math
from enthought.mayavi.mlab import * 
    
def sphere(x1,x2):
    (X1,X2) = numpy.meshgrid(x1,x2)
    cosx1 = numpy.cos(X1)
    sinx1 = numpy.sin(X1)
    cosx2 = numpy.cos(X2)
    sinx2 = numpy.sin(X2)
    xu1   = cosx1*cosx2
    xu2   = sinx1*cosx2
    xu3   = sinx2
    return mesh(xu1,xu2,xu3,colormap='gist_earth')

def cylinder(x1,x2,u2shft):
    (X1,X2) = numpy.meshgrid(x1,x2)
    xu1 = numpy.cos(X1)
    xu2 = numpy.sin(X1)-u2shft
    xu3   = numpy.tan(X2)
    return mesh(xu1,xu2,xu3,colormap='gist_earth')

def grids(x1,x2,x1_lim,x2_lim,nx,u2shft):
    x1_range = numpy.arange(x1_lim[0],x1_lim[1]+0.01,(x1_lim[1]-x1_lim[0])/float(nx))
    x2_range = numpy.arange(x2_lim[0],x2_lim[1]+0.01,(x2_lim[1]-x2_lim[0])/float(nx))
    cosx1_range = numpy.cos(x1_range)
    sinx1_range = numpy.sin(x1_range)
    cosx2_range = numpy.cos(x2_range)
    sinx2_range = numpy.sin(x2_range)
    lines = []
    sf = 1.001
    for x in x2:
        cosx2 = math.cos(x)
        sinx2 = math.sin(x)
        tanx2 = math.tan(x)
        u1 = sf*cosx2*cosx1_range
        u2 = sf*cosx2*sinx1_range
        u3 = sf*sinx2*numpy.array(len(u1)*[1.0])
        lines.append(plot3d(u1,u2,u3,line_width=0.5,color=(0.0,0.0,0.0),tube_radius=None))
        up1 = sf*cosx1_range
        up2 = sf*sinx1_range-u2shft
        up3 = tanx2*numpy.array(len(up1)*[1.0])
        lines.append(plot3d(up1,up2,up3,line_width=0.5,color=(0.0,0.0,0.0),tube_radius=None))
    for x in x1:
        cosx1 = math.cos(x)
        sinx1 = math.sin(x)
        u1 = sf*cosx1*cosx2_range
        u2 = sf*sinx1*cosx2_range
        u3 = sf*sinx2_range
        lines.append(plot3d(u1,u2,u3,line_width=0.5,color=(0.0,0.0,0.0),tube_radius=None))
        up1 = sf*cosx1*numpy.array([1.0,1.0])
        up2 = sf*sinx1*numpy.array([1.0,1.0])-u2shft
        up3 = numpy.array([math.tan(x2_range[0]),math.tan(x2_range[-1])])
        lines.append(plot3d(up1,up2,up3,line_width=0.5,color=(0.0,0.0,0.0),tube_radius=None))
    return(lines)

def tangents(x1,x2,u2shft):
    cosx1 = math.cos(x1)
    sinx1 = math.sin(x1)
    cosx2 = math.cos(x2)
    sinx2 = math.sin(x2)
    tanx2 = math.tan(x2)
    e1 = [-sinx1,cosx1,0.0]
    e2 = [-sinx2*cosx1,-sinx2*sinx1,cosx2]
    x  = [cosx1*cosx2,sinx1*cosx2,sinx2]
    ep1 = e1
    ep2 = [0.0,0.0,1.0/cosx2**2]
    xp  = [cosx1,sinx1-u2shft,tanx2]
    t  = quiver3d([x[0],x[0],xp[0],xp[0]],[x[1],x[1],xp[1],xp[1]],[x[2],x[2],xp[2],xp[2]],\
                  [e1[0],e2[0],ep1[0],ep2[0]],[e1[1],e2[1],ep1[1],ep2[1]],[e1[2],e2[2],ep1[2],ep2[2]],\
                  line_width=2.0,color=(1.0,0.0,0.0)) 
    return(t)
   
x1 = numpy.arange(0.0,2*math.pi+0.01,math.pi/50.0)
x2 = numpy.arange(0.5-math.pi/2.0,math.pi/2.0-0.49,(math.pi-1.0)/50.0)

X1 = numpy.arange(x1[0],x1[-1]+0.01,(x1[-1]-x1[0])/10.0)
X2 = numpy.arange(x2[0],x2[-1]+0.01,(x2[-1]-x2[0])/10.0)

s1 = sphere(x1,x2)
c1 = cylinder(x1,x2,3.0)
g1 = grids(X1,X2,[0.0,2.0*math.pi],[x2[0],x2[-1]],50,3.0)
t1 = tangents(X1[6],X2[8],3.0)

#CR = raw_input()
#savefig('mercator3.png', size=(2048,2048))


