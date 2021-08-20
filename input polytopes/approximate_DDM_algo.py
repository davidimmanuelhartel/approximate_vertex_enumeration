#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import linprog

'''
Implementation of Algorithmus 2 - "Approximate Vertex Enumeration" [LÖHNE 2020]
input: 
    Polytope P = {x in IR^d | Ax <= e}, given by m*d-matrix A - array of dimension (m,d)
    epsilon > 0
output:
    epsilon-approximate V-representation V of P

THIS VERSION WORKS ONLY FOR DIMENSION = 3 
(as the function linearly_independent_rows is only implemented for d = 3)
'''

def compute_vertex_enumeration(A,epsilon):
    # compute the start solution
    V, A_independent, A_rest = startsolution(A, epsilon)
    # A_rearranged = np.concatenate((A_independent, A_rest), axis = 0)
    incidencelist = []
    
    # initialize the incidencelist:
    for v in V:
        v_inequalities = A_independent.dot(v) >= 1  +1e-10
        # colums represent the inequalities, rows the vertices
        incidencelist.append(v_inequalities.tolist())
    
    for row in A_rest:
        # print("row",row)
        V = cut(incidencelist, epsilon,row,V)
        
    return V
    
    
   
def cut(incidencelist,epsilon,h,V):
    dimension = len(h)
    
    V_plus = []  # in H_plus = {x in IR^d | hx > 1 + epsilon }
    V_minus = [] # in H_minus = {x in IR^d | hx < 1 }
    V_null = []          # in H_null = {x in IR^d | hx < 1 }
    # divide up all vertices and 
    # change incidencelist according to the newly added inequality h:
    # print("len(V)=",len(V))
    for i in range(len(V)):
        if (np.dot(h,V[i]) > 1 + epsilon):
            V_plus.append(i)
            # incidencelist[i].append(True)
            incidencelist[i].append(False)
        elif (np.dot(h,V[i]) < 1):
            V_minus.append(i)
            incidencelist[i].append(False)
        else:
            incidencelist[i].append(True)
            V_null.append(i)
    
    # Every pair of vertices which is divided by h will be checked for a common edge
    # that is, the corresponding arrays in the incidence-list will be checked for 
    # |J_>=(v_1) \cap J_>=(v_2)| >= d-1
    # in this case - a new vertex will be added
    # print("V_plus", V_plus,"\nV_minus",V_minus)
    for ind_p1 in V_minus:
        for ind_p2 in V_plus:
            # print(ind_p1,ind_p2)
            I_point1 = incidencelist[ind_p1]
            I_point2 = incidencelist[ind_p2]
            I_point1_AND_I_point2 = [a and b for a, b in zip(I_point1, I_point2)]
            if I_point1_AND_I_point2.count(True) >= dimension-1:
                point1 = V[ind_p1]
                point2 = V[ind_p2]
                # print(point1,point2)
                c = 1 + epsilon/2
                point_new = intersection(point1,point2,h,c)
                V.append(point_new)
                row_new = [a and b for a, b in zip(incidencelist[ind_p1], incidencelist[ind_p2])] 
                row_new[-1] = True
                incidencelist.append(row_new)
        
    # remove points in V_plus from V and the corresponding row in incidencelist
    for index in sorted(V_plus, reverse=True):
        del incidencelist[index]
        del V[index]
        
    return V

def remove_point(V,point):
    ind = 0
    size = len(V)
    print(point)
    while ind != size and not (V[ind] == point).all():
        ind += 1
    if ind != size:
        V.pop(ind)
    else:
        raise ValueError('the list does not contain the point.')



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
    A_independent = np.append(A_independent, [h_new],axis=0)

    
    # compute intersection between d of the d+1 hyperplanes in A_independent
    # for all d possible arrangementsa 
    maske = np.ones(dimension+1, dtype=bool) #Mask, to hide a row
    for i in range(dimension+1):
        maske[i] = False
        intersect = np.linalg.solve(A_independent[maske,], np.ones(dimension))
        V.append((1+epsilon/2)*intersect)
        maske[i] = True

    return V, A_independent, A_rest


#********************************* functions *********************************#

#works only for dimension = 3
def linearly_independent_rows(A):
    # determines 3 linear independent rows for a (n x 3) matrix
    # --> needs to be adapted for arbitrary dimensions
    # idea: start with the first three rows of A. If those are not linearly 
    # independent, take another row instead. Do this with all possible arrangements 
    # of three rows in A, until 3 independent rows are chosen.
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




def intersection(point1,point2,h,c):
    '''
    calculates the intersection between the line between point1 and point 2 and
    the hyperplane H = {x in IR^d | hx = c}, giben by h in IR^d and c in IR
    '''
    
    # Einsetzen der Geraden point1 + s(point2 - point1) in hx=c und Umstellen nach s ergibt:
    s = (c-np.dot(h,point1))/(np.dot(h,(point2-point1)))
    # Resubstitution ergibt den SP:
    intersection = point1 + s*(point2 - point1)
    #print("intersection=",intersection)
    return intersection

#### old functions ###
def cut_v0(incidencelist,epsilon,h,V):
    dimension = len(h)
    
    V_plus = []  # in H_plus = {x in IR^d | hx > 1 + epsilon }
    V_minus = [] # in H_minus = {x in IR^d | hx < 1 }
    V_null = []  # in H_null = {x in IR^d | hx < 1 }
    # divide up all vertices and 
    # change incidencelist accordingly for the newly added inequality h:
    for i in range(len(V)):
        if (np.dot(h,V[i]) > 1 + epsilon):
            V_plus.append(V[i])
            incidencelist[i].append(True)
        if (np.dot(h,V[i]) < 1):
            V_minus.append(V[i])
            incidencelist[i].append(False)
        if (np.dot(h,V[i]) <= 1 + epsilon) and (h.dot(V[i]) >= 1):
            V_null.append(V[i])
            incidencelist[i].append(True)
    
    # Every pair of vertices which is divided by h will be checked for a common edge
    # that is, the corresponding arrays in the incidence-list will be checked for 
    # |J_>=(v_1) \cap J_>=(v_2)| >= d-1
    # in this case - a new vertex will be added
    
    for point1 in V_minus:
        for point2 in V_plus:
            ind_p1 = find_row_indexes(V,point1)
            ind_p2 = find_row_indexes(V,point2)
            I_point1 = incidencelist[ind_p1]
            I_point2 = incidencelist[ind_p2]
            I_point1_AND_I_point2 = [a and b for a, b in zip(I_point1, I_point2)]
            if I_point1_AND_I_point2.count(True) >= dimension-1:
                c = 1 + epsilon/2
                point_new = intersection(point1,point2,h,c)
                V.append(point_new)
                row_new = [a and b for a, b in zip(incidencelist[ind_p1], incidencelist[ind_p2])] 
                row_new[-1] = True
                incidencelist.append(row_new)
            
    # remove points in V_plus from V and the corresponding row in incidencelist
    for point in V_plus:
        ind = find_row_indexes(V,point)
        del incidencelist[ind]
        remove_point(V,point)
        
    return V


def find_row_indexes(matrix,vektor):
    '''
    Given a 2-dim np.array (matrix) and a 1-dim np.array (vektor) the function 
    returns the index of the vektor within the matrix
    '''
    result = (np.where(np.all(matrix == vektor, axis = 1)))
    # print(result)
    return result[0][0]


