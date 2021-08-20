# -*- coding: utf-8 -*-
"""
Implementation of the halfedge data structure and some functions for the graptheoretica algorithm for 
approximate vertex enumeration

"""
import numpy as np

class HalfEdge:
    """Implements a halfedge of a heds"""
    def __init__(self, origin, target):
        self.origin = origin #reference to the origin vertex
        self.target = target #reference to the target vertex
        self.twin = None #reference to the twin halfedge
        self.face = None #reference to the incident face
        self.next = None #reference to the next halfedge along the boundary of the face
        self.prev = None #reference to the preceding halfedge along the boundary of the face
        
    def __eq__(self, other):
        """returns True if equal to other"""
        if isinstance(other, HalfEdge):
            return (self.origin == other.origin) & (self.target == other.target)    
        return False
    def __hash__(self):
        """returns the hash value (integer)"""
        return hash((self.origin,self.target))
    def __str__(self):
        return 'h({0},{1})'.format(self.origin,self.target)
    def __repr__(self):
        return 'h({0},{1})'.format(self.origin,self.target)

        
class Vertex:
    """Implements a vertex of a heds"""
    def __init__(self,ind,coord=None):
        self.ind = ind
        self.halfedge = None# reference to arbitrary halfedge with vertex as origin 
        self.coord = coord  # store the IR^3 coordinates
        self.in_Z = False #flag if Vertex will be cut out

    def __eq__(self, other):
        """returns True if equal to other"""
        if isinstance(other, Vertex):
            return self.ind == other.ind
        return False
    def __hash__(self):
        """returns the hash value (integer)"""
        return hash(self.ind)
    def __str__(self):
        return 'v{}'.format(self.ind)
    def __repr__(self):
        return 'v{}'.format(self.ind)
    
    def adjacent_vertices(self):
    # returns all adjacent vertices for one given vertex
        vertexlist = []
        s = self.halfedge
        h = s
        while True:
            vertexlist.append(h.target)
            h = h.twin.next
            if h == s:
                break
        return vertexlist
    
    def halfedges(self):
    # returns a list all adjacent edges with vertex as target
        edgelist = []
        s = self.halfedge.twin
        h = s
        while True:
            edgelist.append(h)
            h = h.next.twin
            if h == s:
                break
        return edgelist
    
    def faces(self):
    # returns all adjacent faces with vertex on the boundary
        facelist = []
        s = self.halfedge.twin
        h = s
        while True:
            facelist.append(h.face)
            h = h.next.twin
            if h == s:
                break
        return facelist
    

    
class Face:
    """Implements a face of a heds"""
    def __init__(self,facenumber):  
        self.halfedge = None           #reference to an arbitrary halfedge of the face
        self.facenumber = facenumber 
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Face):
            return self.facenumber == other.facenumber
        return False
    def __str__(self):
        return 'f{}'.format(self.facenumber)
    def __repr__(self):
        return 'f{}'.format(self.facenumber)
    def __hash__(self):
        """returns the hash value (integer)"""
        return hash(self.facenumber)
    def vertices(self):
        # returns list containing all vertices of the Face in clockwise order 
        vertices = []
        s = self.halfedge
        h = s
        while True:
            vertices.append(h.target)
            h = h.next
            if h == s:
                break
        return vertices
        
    def halfedges(self):
        # returns list containing all halfedges of the Face  
        halfedges = []
        s = self.halfedge
        h = s
        while True:
            halfedges.append(h)
            h = h.next
            if h == s:
                break
        return halfedges

        

      
class heds:
    """
    Implements a doubly-connected edge list         
    given a list of basis cycles [i1,i2,...,ik,i1] with indices ij of 
    verxes vj in counter clockwise order.
    basis cycles do not contain any sub-cycles)
    """
    def __init__(self,faces = []):
        self.vertices = []      #list of all vertices
        self.halfedges = []     #list of all halfedges
        self.faces = []         #list of all faces
        self.facenumber = 0     #amount of added faces (!= actual amount of faces)
        self.verticenumber = 0
        #self.outer_face = None  #face, containing all others --> UNPRÄZISE!!
        self.cut_vertices = [] 
        self.adjacency = self.vertex_adjacencies        # dict of adjacent vertices for each vertex
        
        if self.faces != []:
            self.build_heds(faces)


    def __str__(self):
        return 'Vertices:{0}\nHalfedges:{1}\nFaces:{2})'.format(self.vertices,
                                                                self.halfedges,
                                                                self.faces)
    def __repr__(self):
        return 'Vertices:{0}\nHalfedges:{1}\nFaces:{2})'.format(self.vertices,
                                                                self.halfedges,
                                                                self.faces)
        # some functions:
    def vertex_adjacencies(self):
        # returns a dict with all adjacent vertices in ccw order for each vertex 
        adj = {}
        for vertex in self.vertices:
            vertices = vertex.adjacent_vertices
            adj["{}".format(vertex)] = vertices
        return adj
    
    def face_vertices(self):
        #returns a list containg vertices in cw order for every face 
        face_vertices_list = []
        for face in self.faces:
            face_vertices_list.append(face.vertices())
        return face_vertices_list
            
    def find_vertex(self,ind):
        #returns vertex object given the index number
        for vertex in self.vertices:
            if ind == vertex.ind:
                return vertex
        return None
        
    def find_halfedge(self, origin, target):
        # returns Halfedge, given origin and target vertices
        for edge in self.halfedges:
            if (edge.prev.target == origin and edge.target == target):
                return edge
        return None
    
    def update_indexes(self):
        # updates the vertex indexes, so that the maximum index
        # matches the number of vertices
        ind = 0
        for vertex in self.vertices:
            vertex.ind = ind
            ind = ind+1


        
    def heds_to_off_file(self,n):
        number_of_vertices = len(self.vertices)
        number_of_faces = len(self.faces)
        number_of_edges = len(self.halfedges)//2
        f = open("off_output_{}.off".format(n), 'w')
        f.write("OFF\n\n")
        f.write("{0} {1} {2}\n".format(number_of_vertices,number_of_faces,number_of_edges))
        for vertex in self.vertices:
            f.write("{0} {1} {2}\n".format(vertex.coord[0],vertex.coord[1],vertex.coord[2]))
        face_vertices_list = self.face_vertices()
        for face in face_vertices_list:
            f.write("{} ".format(len(face)))
            for vertex in face:
                f.write("{} ".format(vertex.ind))
            f.write("\n")
    """
    Euler-Operators
    """
    def kill_edge_face(self,edge1):
        '''
        removes both h1 and h1.twin and thus merges both faces into h1.face 
        '''
        # if edge1 not in self.halfedges:
        #     return
    
        edge2 = edge1.twin
        face = edge1.face #face to keep
        face_to_remove = edge2.face
        #assign face to the edges of the face to remove:
        edge = edge2
        startedge = edge
        while(True):
            edge.face = face
            edge = edge.next
            if(edge == startedge):
                break
        #reassign prev and next 
        edge1.prev.next = edge2.next
        edge2.next.prev = edge1.prev
        edge1.next.prev = edge2.prev 
        edge2.prev.next = edge1.next 
        
        if edge1.origin.halfedge == edge1:
            edge1.origin.halfedge = edge2.next
        if edge2.origin.halfedge == edge2:
            edge2.origin.halfedge = edge1.next
        if edge1.face.halfedge == edge1:
            edge1.face.halfedge = edge1.prev
        if edge2.face.halfedge == edge2:
            edge2.face.halfedge = edge2.next
        
        # if (face_to_remove in self.faces) and (face != face_to_remove):
        if (face != face_to_remove):
            try:
                self.faces.remove(face_to_remove) #in case the face already has been removed
            except:
                pass
        self.halfedges.remove(edge1)
        self.halfedges.remove(edge1.twin)
        
    
    def kill_edge_vertex(self,edge1): #argument is edge adjacent to face to be removed
        '''
        removes edge1,  edge1.twin and edge1.target
        '''
        #reassing prev and next
        edge1.prev.next = edge1.twin.next
        edge1.twin.next.prev = edge1.prev
        #check if edge.origin.halfedge == edge:
        if edge1.origin.halfedge == edge1:
            edge1.origin.halfedge = edge1.twin.next
        if (edge1.face.halfedge == edge1) or (edge1.face.halfedge == edge1.twin):
            edge1.face.halfedge = edge1.twin.next
        self.vertices.remove(edge1.target)
        
        #in case the halfedges already have been removed
        try:
            self.halfedges.remove(edge1)
            self.halfedges.remove(edge1.twin)
        except:
            pass


        

    def split_edge(self,edge1,hyperplane, epsilon, new_coord):
        '''
        adds a new vertex on edge1 and connects it with edge1.origin and edge1.target
        then, the edge will be replaced by the new halfedges. If coordinates
        are used, the cutting point of the straight line between
        edge1.origin and edge1.target and the hyperplane will be assigned
        '''
        v1 = edge1.origin
        v2 = edge1.target
        edge2 = edge1.twin
        # vertexnumber = len(self.vertices)
        # v_new = Vertex(vertexnumber+1)
        v_new = Vertex(self.verticenumber)
        self.verticenumber += 1
        e1 = HalfEdge(v1, v_new)
        e2 = HalfEdge(v_new, v1)
        f1 = HalfEdge(v_new, v2)
        f2 = HalfEdge(v2, v_new)
        # if v1 or v2 store referenes to halfedges which we want to remove, 
        # assign them to the new halfedges
        if v1.halfedge.target == v2:
            v1.halfedge = e1
        if v2.halfedge.target == v1:
            v2.halfedge = f2
        v_new.halfedge = e2
        # assign twin, face, next, prev for all halfedges
        e1.twin = e2
        e2.twin = e1
        f1.twin = f2
        f2.twin = f1
        e1.face = f1.face = edge1.face
        e2.face = f2.face = edge2.face
        edge1.face.halfedge = e1  # in case edge1.face.halfedge == edge1
        edge2.face.halfedge = e2
        e1.next =f1
        if edge1.next != edge2: #check, if edge1.target is external/leaf node
            f1.next = edge1.next
        else: f1.next = f2
        f2.next = e2
        if edge2.next != edge1: #similarily check, if edge2.target is external/leaf node
            e2.next = edge2.next
        else: e2.next = e1
        f1.prev = e1
        e1.prev = edge1.prev 
        e2.prev = f2
        f2.prev = edge2.prev
        edge1.prev.next = e1
        edge2.prev.next = f2
        edge1.next.prev = f1
        edge2.next.prev = e2
        #assign new coordinates, dependent if v1 in V_minus or V_null
        if new_coord:
            point1 = v1.coord
            point2 = v2.coord
            c = 1 + epsilon/2
            s = (c-np.dot(hyperplane,point1))/(np.dot(hyperplane,(point2-point1)))
            cutting_point = point1 + s*(point2 - point1)
            v_new.coord = cutting_point
        else:
            v_new.coord = v1.coord
        self.halfedges.extend([e1,e2,f1,f2])
        self.vertices.append(v_new)
        try: #in case the halfedges already have been removed
            self.halfedges.remove(edge1)
            self.halfedges.remove(edge2)
        except:
            pass
        
    def make_edge_vertex(self,edge1):
        '''
        adds a new vertex to a face and connects it with edge1.target
        '''
        face = edge1.face
        vertexnumber = len(self.vertices)
        vertex_new = Vertex(vertexnumber+1)
        h1 = HalfEdge(edge1.target, vertex_new)
        h2 = HalfEdge(vertex_new, edge1.target)
        vertex_new.halfedge = h2
        h1.twin = h2
        h2.twin = h1
        h1.face = h2.face = face
        h1.next = h2
        h2.next = edge1.next
        h1.prev = edge1
        h2.prev = h1
        edge1.next = h1
        h2.next.prev = h2
        self.vertices.append(vertex_new)
        self.halfedges.extend([h1,h2])
        return vertex_new

    def make_edge_face(self,vertex1,vertex2,face):
        '''
        connects vertex1 and vertex2, provided vertex1 and vertex2 are
        are incident to the given face. The face will be assigned to 
        h(vertex1,vertex2). A new face will be assigned to h(vertex2, vertex1)
        '''
        if (vertex1 == vertex2):
            return
        # in order to assign next and prev for new halfedges, we need 
        # the halfedges in the face with vertex1, vertex2 as target
        for halfedge in face.halfedges():
            if halfedge.target == vertex1:
                e1 = halfedge
            elif halfedge.target == vertex2:
                e2 = halfedge
                
        h1 = HalfEdge(vertex1,vertex2)
        h2 = HalfEdge(vertex2,vertex1)
        # the new halfedges split the facee into two new faces
        self.facenumber +=1 
        f_new = Face(self.facenumber)
        face.halfedge = h1
        f_new.halfedge = h2
    
        h1.twin = h2
        h2.twin = h1
        
        h1.next = e2.next
        h2.next = e1.next
        h1.next.prev = h1
        h2.next.prev = h2
        e1.next = h1
        e2.next = h2
        h1.prev = e1
        h2.prev = e2
        #assign new faces:
        e = h1
        while(True):
            e.face = face
            if e.target == vertex1:
                break
            e = e.next
        e = h2
        while(True):
            e.face = f_new
            if e.target == vertex2:
                break
            e = e.next
        self.faces.append(f_new)
        self.halfedges.extend([h1,h2])
        
    def make_edge_face1(self,e1,e2):
        '''
        connects e1.target and e2.target. The existing face will be assigned to all edges of the face of e1, 
        whereas a new face will be created and reassigned to all other edges of the face of e2
        '''
        face = e1.face
        if (e1.target == e2.target) or (e1 == e2):
            return
        h1 = HalfEdge(e1.target,e2.target)
        h2 = HalfEdge(e2.target,e1.target)
        # the new halfedges split the facee into two new faces
        self.facenumber +=1 
        f_new = Face(self.facenumber)
        face.halfedge = h1
        f_new.halfedge = h2
    
        h1.twin = h2
        h2.twin = h1
        
        h1.next = e2.next
        h2.next = e1.next
        h1.next.prev = h1
        h2.next.prev = h2
        e1.next = h1
        e2.next = h2
        h1.prev = e1
        h2.prev = e2
        #assign new faces:
        e = h1
        while(True):
            e.face = face
            if e.target == e1.target:
                break
            e = e.next
        e = h2
        while(True):
            e.face = f_new
            if e.target == e2.target:
                break
            e = e.next
        self.faces.append(f_new)
        self.halfedges.extend([h1,h2])
        return h2
        
    def build_heds(self,faces):
        """
        Implements a doubly-connected edge list         
        given a list of faces [v1,v2,...,vk,v1] with integer indices 
        or instances of vertices vj (coord = True) in counter clockwise order.
        basis cycles do not contain any sub-cycles)
        """
        for face in faces:
            vertexlist = []
            edgelist = []
            l = len(face)
            self.facenumber += 1 
            new_face = Face(self.facenumber)
            
            #assign face vertices:
            if type(face[0]) == int:  #indexes are given
                for i in range(l):
                    v = Vertex(face[i])
                    vertexlist.append(v)
            else:   # Vertex instances are given
                vertexlist = face
            
            #assign face edges:
            for i in range(l):
                if i == l-1:
                    v1 = vertexlist[i]
                    v2 = vertexlist[0]
                else:
                    v1 = vertexlist[i]
                    v2 = vertexlist[i+1]
                h = HalfEdge(v1,v2)
                h.face = new_face
                edgelist.append(h)
                # assign halfedge for vertex:
                v1.halfedge = h
            # assign prev and next for every edge:
            for i in range(len(edgelist)):
                if i == l-1:
                    e1 = edgelist[i]
                    e2 = edgelist[0]
                else:
                    e1 = edgelist[i]
                    e2 = edgelist[i+1]
                e1.next = e2
                e2.prev = e1
            # assign an arbitrary halfedge to the face :
            new_face.halfedge = edgelist[0]
            self.vertices.extend(vertexlist)
            self.halfedges.extend(edgelist)
            self.faces.append(new_face)

            
        # assign twins for every edge 
        new_halfedges = []
        self.facenumber += 1
        self.outer_face = Face(self.facenumber)

        for halfedge in self.halfedges:
            v1 = halfedge.origin
            v2 = halfedge.target
            twin = self.find_halfedge(v2,v1)
            # assign twin, if it is not already assigned
            if twin != None:
                halfedge.twin = twin
                twin.twin = halfedge
            # if there is no twin yet, since it is an outer edge
            # add halfedge and assign to outer face
            else:
                new_halfedge = HalfEdge(v2,v1)
                halfedge.twin =new_halfedge
                new_halfedge.twin = halfedge
                new_halfedge.face = self.outer_face
                new_halfedges.append(new_halfedge)
                
        # assign next and prev for outer new added halfedges:
        for halfedge1 in new_halfedges:
            for halfedge2 in new_halfedges:
                if (halfedge1.target == halfedge2.twin.target):
                    halfedge1.next = halfedge2
                    halfedge2.prev = halfedge1
        
        # assign an arbitrary halfedge to the face 
        self.outer_face.halfedge = new_halfedges[0]
        # finally, append halfedges and face to heds:
        self.halfedges.extend(new_halfedges)
        self.faces.append(self.outer_face)
        self.vertices = list(set(self.vertices)) #remove duplicates in the list
        self.verticenumber = len(self.vertices) # to count the total number of vertives
      

def cut(G,h,epsilon=0):
    """ Implementation of the cutting procedure from Loehne's 
    "Approximate Vertex Enumeration", p. 14 
    input: heds of Graph G, list Z of vertices to be cut, and epsilon to 
    determine the coordinates of the new added points.
    output: Graph G` produced by the cut
    """

    # cut_vertices = Z = H_plus = {x in IR^d | hx > 1 + epsilon }
    V_minus = [] # H_minus = {x in IR^d | hx < 1 }
    # V_null = [] # H_null = {x in IR^d | hx < 1 }
    
    Z = []
    cut_faces = [] #remember all faces affected by the cut
    for v in G.vertices: 
        if (np.dot(h,v.coord) > 1 + epsilon):
            v.in_Z = True #flag determinating if vertex will be cut
            Z.append(v)
            for face in v.faces():
                if face not in cut_faces:
                    cut_faces.append(face)
        elif (np.dot(h,v.coord) < 1):
            V_minus.append(v)
        # if (np.dot(h,v.coord) <= 1 + epsilon) and (h.dot(v.coord) >= 1):
        #     V_null.append(v)
    
    if Z == []:
        return
    
    # find all halfedges, which will be cut
    # that is, all halfedges in G with one and only one vertex in Z 
    
    cut = [] #edges, which will be cut
    
    for face in cut_faces:
        edge = face.halfedge
        start_vertex = edge.target
        while(True):
            if (edge.target.in_Z == True) and (edge.origin.in_Z == False):
                cut.append(edge)
                # create new vertices and assign the coordinates by split_edge
                if edge.origin in V_minus:
                    G.split_edge(edge,h,epsilon,True)
                else: 
                    G.split_edge(edge,h,epsilon,False)
            edge = edge.next
            
            if edge.target == start_vertex:
                break   

    # find all neighbors of the new vertices for every face
    face_neighbors = []
    for face in cut_faces:
        neighbors = []
        neighbor1 = neighbor2 = None
        start_edge = face.halfedge # start with arbitrary halfedge an then pass all face-edges
        next_start_edge = None # ... and remember a halfedge with target as next start_edge
        edge = start_edge
        # find halfedge with target vertex in Z:
        while(True):
            if (edge.target.in_Z == True) and (edge.origin.in_Z == False):
                next_start_edge = edge
                break
            edge = edge.next
            if edge==start_edge:
                break
        if next_start_edge != None: #that is, there exists a cutting edge in the face
            edge = next_start_edge
            while(True):
                if (edge.target.in_Z == True) and (edge.origin.in_Z == False):
                    neighbor1 = edge.prev
                if (edge.target.in_Z == False) and (edge.origin.in_Z == True):
                    neighbor2 = edge
                    if neighbor1 != None:
                        neighbors.append([neighbor1,neighbor2])
                        neighbor1 = None
                        neighbor2 = None
                edge = edge.next
                if edge == next_start_edge:
                    break
        face_neighbors.append(neighbors)
    # connect all neighbors:
    connected_vertices = []
    for neighbors in face_neighbors:
        for nbr_pair in neighbors:
            edge1 = nbr_pair[0]
            edge2 = nbr_pair[1]
            v1 = edge1.target
            v2 = edge2.target

            # if there is already an edge between v1 and v2, create a new 
            # vertex to maintain a simple graph
            if HalfEdge(v1,v2) in connected_vertices or HalfEdge(v2,v1) in connected_vertices: 
                G.make_edge_vertex(edge1)
                G.make_edge_face1(edge1.next,edge2)
            else:
                new_halfedge = G.make_edge_face1(edge1,edge2)
                connected_vertices.append(new_halfedge)
     
    # remove vertices
    for vertex in Z:
        halfedges = vertex.halfedges()
        l = len(halfedges)
        for i in range(l-1):
            G.kill_edge_face(halfedges[i]) #first, remove edges until one
        G.kill_edge_vertex(halfedges[l-1]) #second, remove the last edge
        if vertex in G.vertices: # in case all halfedges have been removed already
            G.vertices.remove(vertex)
            
    # in case the by the removal of the vertices merged face has been assigned 
    # to a deleted halfedge, we reassign the reference to a existing halfedge.
    # as we might have cut the graph at more than one distinct areas,
    # we need to reassign assign all faces. as we don't know the exact amount of faces
    # we do it via brute force: assign it for every edge new.
    for edge in connected_vertices:
        edge.face.halfedge = edge 



