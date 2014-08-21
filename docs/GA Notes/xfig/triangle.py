import numpy
import delaunay.core as core

set = [(0,0),(2,1),(1,2),(4,0),(0,4),(4,4)]
T = core.Triangulation(set)

print T.get_elements()
print T.get_elements_indices()
print T.get_neighbours()
print T.get_set()
