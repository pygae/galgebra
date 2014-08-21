import numpy,math,random,sys
from enthought.mayavi.mlab import * 
from shapely.geometry import Polygon
from shapely.geometry import Point
import delaunay.core as core
  
def sphere(x1,x2):
    (X1,X2) = numpy.meshgrid(x1,x2)
    cosx1 = numpy.cos(X1)
    sinx1 = numpy.sin(X1)
    cosx2 = numpy.cos(X2)
    sinx2 = numpy.sin(X2)
    xu1   = cosx1*cosx2
    xu2   = sinx1*cosx2
    xu3   = sinx2
    return mesh(xu1,xu2,xu3,colormap='gist_earth',opacity=0.25)

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
        lines.append(plot3d(u1,u2,u3,line_width=1.5,color=(0.0,0.0,0.0),tube_radius=None))
    for x in x1:
        cosx1 = math.cos(x)
        sinx1 = math.sin(x)
        u1 = sf*cosx1*cosx2_range
        u2 = sf*sinx1*cosx2_range
        u3 = sf*sinx2_range
        lines.append(plot3d(u1,u2,u3,line_width=1.5,color=(0.0,0.0,0.0),tube_radius=None))
    return(lines)

def manifold_curve(X1,X2,delta):
    R = 1.0+delta
    npts = len(X1)
    u1 = numpy.zeros(npts)
    u2 = numpy.zeros(npts)
    u3 = numpy.zeros(npts)
    for i in range(npts):
        x1 = X1[i]
        x2 = X2[i]
        cos_x1 = math.cos(x1)
        sin_x1 = math.sin(x1)
        sin_x2 = math.sin(x2)
        cos_x2 = math.cos(x2)
        u1[i] = R*cos_x1*cos_x2
        u2[i] = R*sin_x1*cos_x2
        u3[i] = R*sin_x2
    return(plot3d(u1,u2,u3,line_width=0.5,color=(0.0,0.0,0.0),tube_radius=None))
    
def ellipse(xmin,xmaj,npts,mode=True):
    xmin /= 2.0
    xmaj /= 2.0
    if mode:
        X1 = numpy.zeros(npts)
        X2 = numpy.zeros(npts)
    else:
        pts = []
    dth = 2.0*math.pi/float(npts-1)
    for i in range(npts):
        th = float(i)*dth
        cos_th = math.cos(th)
        if mode:
            X1[i] = xmaj*cos_th
            X2[i] = xmin*(0.25+cos_th**2)*math.sin(th)
        else:
            pts.append((xmaj*cos_th,xmin*(0.25+cos_th**2)*math.sin(th)))
    if mode:
        return(X1,X2)
    else:
        return(pts)

def vector_map(pts,delta=0.0):
    R = 1.0+delta
    npts = len(pts)
    u1 = numpy.zeros(npts)
    u2 = numpy.zeros(npts)
    u3 = numpy.zeros(npts)
    i = 0
    for pt in pts:
        x1 = pt[0]
        x2 = pt[1]
        cos_x1 = math.cos(x1)
        sin_x1 = math.sin(x1)
        sin_x2 = math.sin(x2)
        cos_x2 = math.cos(x2)
        u1[i] = R*cos_x1*cos_x2
        u2[i] = R*sin_x1*cos_x2
        u3[i] = R*sin_x2
        i += 1
    return(u1,u2,u3)

def shrink_triangle(triangle,delta=0.001):
    u = []
    for pt in triangle:
        u.append(numpy.array(pt))
    u_cg = (u[0]+u[1]+u[2])/3.0
    delta1 = 1.0-delta
    new_triangle = []
    for ui in u:
        new_triangle.append(tuple(u_cg+delta1*(ui-u_cg)))
    return(new_triangle)
    
def triangulate_area_in_closed_curve(bndry_pts,xmax,xmin,ymax,ymin,N,vector_map):
    polygon = Polygon(bndry_pts)
    interior_pts = []
    dx = (xmax-xmin)/float(N+1)
    dy = (ymax-ymin)/float(N+1)
    x0 = xmin+dx
    y0 = ymin+dy
    for i in range(N):
        x = x0+float(i)*dx+random.uniform(-0.1,0.1)
        for j in range(N):
            y = y0+float(j)*dy+random.uniform(-0.1,0.1)
            pt = Point(x,y)
            if polygon.contains(pt):
                interior_pts.append((x,y))
    surface_pts = bndry_pts+interior_pts
    simplices = core.Triangulation(surface_pts)
    triangles = simplices.get_elements()
    pts = simplices.get_set()
    (u1,u2,u3) = vector_map(pts)
    indices = simplices.get_elements_indices()
    reduced_indices = []
    for index in indices:
        x1 = numpy.array(pts[index[0]])
        x2 = numpy.array(pts[index[1]])
        x3 = numpy.array(pts[index[2]])
        xcg = (x1+x2+x3)/3.0
        x1 = xcg+0.999*(x1-xcg)
        x2 = xcg+0.999*(x2-xcg)
        x3 = xcg+0.999*(x3-xcg)
        tmp = Polygon((x1,x2,x3))
        if polygon.contains(tmp):
            reduced_indices.append(index)
    return(triangular_mesh(u1,u2,u3,reduced_indices,line_width=0.25,tube_radius=0.005,representation='fancymesh'))

x1 = numpy.arange(0.0,2*math.pi+0.01,math.pi/50.0)
x2 = numpy.arange(0.5-math.pi/2.0,math.pi/2.0-0.49,(math.pi-1.0)/50.0)

X1 = numpy.arange(x1[0],x1[-1]+0.01,(x1[-1]-x1[0])/10.0)
X2 = numpy.arange(x2[0],x2[-1]+0.01,(x2[-1]-x2[0])/10.0)

#s1 = sphere(x1,x2)
g1 = grids(X1,X2,[0.0,2.0*math.pi],[x2[0],x2[-1]],50,3.0)

xmaj = 1.5
xmin = 2.0

(X1,X2) = ellipse(xmin,xmaj,150)
pts = ellipse(2.0,1.5,30,False)
triangulate_area_in_closed_curve(pts,xmaj/2.0,-xmaj/2.0,xmin/2.0,-xmin/2.0,7,vector_map)
manifold_curve(X1,X2,0.01)

CR = raw_input()
#savefig('manifold_simplex.png', size=(2048,2048))
