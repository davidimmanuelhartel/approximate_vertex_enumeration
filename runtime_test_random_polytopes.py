# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 10:58:34 2021

@author: David
"""
import numpy as np
import cdd
import graphtheoretical_DDM_algo as graph_algo
import approximate_DDM_algo as vert_enum_inc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import time

from matplotlib import pyplot  as plt
from scipy.spatial import ConvexHull

################## functions #########################
def simple_plot(poly_exact,poly_approx_incidence, poly_approx_graph):
    ##### plot convex hull with pyplotlib #####
    # works only with number_type = float #
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title("random polytope with {0} vertices".format(len(poly_exact)))
    for cube, color in zip([poly_approx_incidence, poly_approx_graph,poly_exact], ["b",'y','g']):
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
    
 
def sample_spherical(n, dim=3):
    #computes n random points on the sphere
    vec = np.random.randn(dim, n)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

def random_polytope_runtime(n):
    np.random.seed(seed=2)
    a = sample_spherical(n,3)
    ones = np.ones(np.shape(a)[1])
    vrep = np.append([ones],a,axis=0).T
    
    mat = cdd.Matrix(vrep,number_type = "float")
    mat.rep_type = cdd.RepType.GENERATOR
    poly = cdd.Polyhedron(mat)
    
    hrep = poly.get_inequalities()
    poly_exact = poly.get_generators()
    poly_exact = (np.array([list(x)[1:] for x in poly_exact]))

    A = np.array([list(x)[1:] for x in hrep]) #transform cdd output to np.array
    b = np.array([list(x)[0] for x in hrep]) 

    start = time.time()
    for epsilon in [10,5,1,1e-1,1e-2,1e-5,0]:
        # print("n=",n,"eps=",epsilon)
        # print("eps=",epsilon)
    
        ### results computed with cdd ###
        
        # b_eps = (1+epsilon)*np.ones(A.shape[0])
        # A_b_eps = np.zeros((A.shape[0],A.shape[1]+1))
        # A_b_eps[:,0] = b_eps
        # A_b_eps[:,1:] = -A
        # mat_eps = cdd.Matrix(A_b_eps,number_type = "float")
        # mat_eps.rep_type = cdd.RepType.INEQUALITY
        # poly_eps = cdd.Polyhedron(mat_eps)
        # poly_eps = poly_eps.get_generators()
        # poly_eps = np.array([list(x)[1:] for x in poly_eps]) #transform cdd output to np.array
        # print("poly_eps",poly_eps)
        
        # print("poly_exact",poly_exact)
    
        ### graphtheoretical_DDM_algo ###
        poly_approx_graph,_ = graph_algo.compute_vertex_enumeration(-A/b[:,None],epsilon)
        # print("poly_approx_graph",poly_approx_graph)
        
        ### approximate_DDM_algo ###
        # poly_approx_incidence = np.array(vert_enum_inc.compute_vertex_enumeration(A,epsilon))
        # print("poly_approx_incidence",poly_approx_incidence)
        
        # print("number of vertices:", len(poly_exact))
        
        
        ### plott output ###
        # simple_plot(poly_exact,poly_approx_incidence, poly_approx_graph)
    end = time.time()
    return end - start

def runtime(n):
    runtime = []
    for i in range(10,n,10):
        t = random_polytope_runtime(i)
        runtime.append([i,t])
    return runtime

# rt = runtime(500)
# np.savetxt("graph_algo.csv", rt, delimiter=",")
