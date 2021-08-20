# -*- coding: utf-8 -*-

'''
Implementation of the graptheoretical algorithm for 
approximate vertex enumeration
input: 
    Polytope P = {x in IR^d | Ax <= e}, given by a m*d-matrix A - e.g. an array of dimension (m,d)
    epsilon > 0
output:
    epsilon-approximate V-representation V of P 
    planar embedding of the polytope as halfedge data structure
'''

import numpy as np
import halfedge_data_structure as heds
from scipy.optimize import linprog #LP solver

def compute_vertex_enumeration(A, epsilon):
    V, A_independent, A_rest = startsolution(A,epsilon)
    # assign faces of vertex
    # since it is a 3-simplex, the order of the vertices is uniquely defined
    v0 = heds.Vertex(0,V[0])
    v1 = heds.Vertex(1,V[1])
    v2 = heds.Vertex(2,V[2])
    v3 = heds.Vertex(3,V[3])
    faces =  [[v1,v2,v3],[v3,v2,v0],[v3,v0,v1]] #only three of four faces are needed by implementation of build_heds
    #[v0,v2,v1] is not needed
    G = heds.heds()
    G.build_heds(faces)
    

    for h in A_rest:
        heds.cut(G,h,epsilon)

    vertice_coords = []
    for vertex in G.vertices:
        vertice_coords.append(vertex.coord)
        
    G.update_indexes()
    return np.array(vertice_coords),G


def startsolution(A,epsilon):
    dimension = A.shape[1]
    number_of_ineqs = A.shape[0]
    V = []
    # we determine the first approximative-rep by the means of d linearly independent
    # rows in A. 

    # find d linearly independent rows in A
    A_independent, A_rest = linearly_independent_rows(A)
    
    # Create a new hyperplane by takinng the negative of the sum of the d linearly
    # independent rows
    # "push" this hyperplane outward until A is included in this system of d+1 inequalities
    # There a LP must be solved:
    # max_x: hx
    # st:    Ax <= e
    h = -A_independent.sum(axis=0)   
    b = np.ones(number_of_ineqs)
    variable_bounds = [(None,None)]*dimension
    x_opt = linprog(-h, A_ub=A, b_ub=b, bounds = variable_bounds, method='revised simplex').x 
    h_new = h/h.dot(x_opt) #newe Hyperebene
    ##print("h_new=",h_new)
    A_independent = np.append(A_independent, [h_new],axis=0)

    
    # compute intersection between d of the d+1 hyperplanes in A_independent
    # for all d possible arrangementsa 
    maske = np.ones(dimension+1, dtype=bool) #Mask, to hide a row
    for i in range(dimension+1):
        maske[i] = False
        intersect = np.linalg.solve(A_independent[maske,], np.ones(dimension))
        V.append((1+epsilon/2)*intersect)
        maske[i] = True

    V = np.array(V)
    return V, A_independent, A_rest


def linearly_independent_rows(A):
    dimension = A.shape[1]
    number_of_ineqs = A.shape[0]
    for i in range(number_of_ineqs-2):
        for j in range(i+1,number_of_ineqs-1):
            for k in range(i+j+1,number_of_ineqs):
                rank = np.linalg.matrix_rank(A[[i,j,k],])
                if(rank == dimension): 
                    independent_inds = (i,j,k)
                    mask = np.invert(np.ones(number_of_ineqs, dtype=bool)) #array of only FALSE 
                    mask[list(independent_inds)] = True
                    A_independent = A[mask,]
                    A_rest = A[[not elm for elm in mask],]
                    return A_independent, A_rest

