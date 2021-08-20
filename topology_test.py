# -*- coding: utf-8 -*-
"""
graph-algorithm - geometry test

"""

import numpy as np
import graphtheoretical_DDM_algo as graph_algo
import approximate_DDM_algo as vert_enum_inc
import cdd

# plotting
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import pyplot  as plt
from scipy.spatial import ConvexHull

################  functions ########################
def simple_plot(poly_exact,poly_approx_incidence, poly_approx_graph, name, epsilon):
    ##### plot convex hull with pyplotlib #####
    # works only with number_type = float #
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title("{0}_eps_{1}".format(name,epsilon))
    # for cube, color in zip([poly_approx_incidence, poly_approx_graph,poly_exact], ['y','g','b']):
    for cube, color in zip([poly_approx_incidence,poly_exact], ['y','b']):
    # for cube, color in zip([poly_approx_graph], ['y']):
        hull = ConvexHull(cube)
        # draw the polygons of the convex hull
        for s in hull.simplices:
            tri = Poly3DCollection(cube[s])
            tri.set_color(color)
            tri.set_alpha(0.4)
            ax.add_collection3d(tri)
        # draw the vertices
        ax.scatter(cube[:, 0], cube[:, 1], cube[:, 2], marker='.', color='b')
    plt.show()


def dcel_to_off_file(dcel, name = ""):
    ##### OFF-File- Output #####
    number_of_vertices = len(dcel.vertices)
    number_of_faces = len(dcel.faces)
    number_of_edges = len(dcel.halfedges)//2
    f = open(r'C:\Users\David\iCloudDrive\Desktop\Bachelorarbeit\Algorithmen\off_files\{}.off'.format(name), 'w')
    f.write("OFF\n\n")
    f.write("{0} {1} {2}\n".format(number_of_vertices,number_of_faces,number_of_edges))
    for vertex in dcel.vertices:
        f.write("{0} {1} {2}\n".format(vertex.coord[0],vertex.coord[1],vertex.coord[2]))
    face_vertices_list = dcel.face_vertices()
    for face in face_vertices_list:
        f.write("{} ".format(len(face)))
        for vertex in face:
            f.write("{} ".format(vertex.ind))
        f.write("\n")

    
################## test #####################
''' files of certain h-reps of polytopes computed with BENSOLVETOOLS '''

# ball
# cube 
# cube+ball 
# cube+ball_polar
# cube+ball_polar+cube
# cube+ball_polar+cube_polar
# bensolvehedron(3,1)
# bensolvehedron(3,1)_polar+ball
# bensolvehedron(3,1)_polar+cube+ball
# bensolvehedron(3,2)

name = "bensolvehedron(3,1)"
vertexnumber_exact = []
vertexnumber_graph = []
vertexnumber_graph_filtered = []
vertexnumber_incidence = []
vertexnumber_incidence_filtered = []

def graph_algorithm_topology_test():
    for epsilon in [10,5,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,1e-1,5e-2,1e-2,5e-3,1e-3,1e-4,1e-5,1e-10,1e-11,1e-12,1e-13,5e-14,3e-14,1e-14,5e-15,3e-15,1e-15,1e-16,1e-17,1e-18,1e-19,0]:
        print("epsilon = ", epsilon)
        
        ### csv-data input ###
        A = np.genfromtxt('A_{}.csv'.format(name), delimiter=' ')
        b = np.genfromtxt('b_{}.csv'.format(name), delimiter=' ')
        A = A/b[:,None]
        b = np.ones(A.shape[0])
        
        # a = [1,-1]
        # A = np.array([list(x) for x in product(a,a,a)])
        # b = np.ones(A.shape[0])
        # print(A)
        
        ### Ball ###
        # a = [1,-1]
        # A = np.array([list(x) for x in product(a,a,a)])
        # b = np.ones(A.shape[0])
        # print(A)
        
        ### Cube ###
        # A = np.array([[-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]])
        # b = np.ones(A.shape[0])
    
        ### simplex ###
        # A = np.array([[-1,0,0],[0,-1,0],[0,0,-1],[1,1,1]])
        # b = np.ones(A.shape[0])
    
        # transform to input-format for cddlib
        b_eps = (1+epsilon)*b
        A_b = np.zeros((A.shape[0],A.shape[1]+1))
        A_b_eps = np.zeros((A.shape[0],A.shape[1]+1))
        A_b[:,0] = b
        A_b_eps[:,0] = b_eps
        A_b[:,1:] = -A
        A_b_eps[:,1:] = -A
        
        
        mat = cdd.Matrix(A_b,number_type = "float")
        mat.rep_type = cdd.RepType.INEQUALITY
        poly = cdd.Polyhedron(mat)
        mat_eps = cdd.Matrix(A_b_eps,number_type = "float")
        mat_eps.rep_type = cdd.RepType.INEQUALITY
        poly_eps = cdd.Polyhedron(mat_eps)
        
        # results computed with cdd:
        # exact polytope P:
        poly_exact = poly.get_generators()
        poly_exact = np.array([list(x)[1:] for x in poly_exact]) #transform cdd output to np.array
        print("poly_exact",poly_exact)
        # epsilon polytope (1+eps)*P:
        poly_eps = poly_eps.get_generators()
        poly_eps = np.array([list(x)[1:] for x in poly_eps]) #transform cdd output to np.array
        print("poly_eps",poly_eps)
        # result computed with graph algorithm:
        poly_approx_graph, poly_approx_graph_dcel = graph_algo.compute_vertex_enumeration(A,epsilon)
        print("poly_approx_graph",poly_approx_graph)
        # results computet with initial algortihm:
        poly_approx_incidence = np.array(vert_enum_inc.compute_vertex_enumeration(A,epsilon))
        print("poly_approx_incidence",poly_approx_incidence)
        
        # counting the number of vertices for each resulting polytope
        vertexnumber_exact.append(len(poly_exact))
        vertexnumber_graph.append(len(poly_approx_graph))
        vertexnumber_graph_filtered.append(len(np.unique(poly_approx_graph,axis=0)))
        vertexnumber_incidence.append(len(poly_approx_incidence))
        vertexnumber_incidence_filtered.append(len(np.unique(poly_approx_incidence,axis=0)))
        
        # visual output
        # dcel_to_off_file(poly_approx_graph_dcel,"{0}_eps_{1}".format(name,epsilon))
        simple_plot(poly_exact,poly_approx_incidence, poly_approx_graph, name, epsilon)
        
        
    print("exact",vertexnumber_exact)
    print("incidence",vertexnumber_incidence)
    print("incidence_filtered",vertexnumber_incidence_filtered)
    print("graph",vertexnumber_graph)
    print("graph_filtered",vertexnumber_graph_filtered)
    
        
graph_algorithm_topology_test()
    

    
    




