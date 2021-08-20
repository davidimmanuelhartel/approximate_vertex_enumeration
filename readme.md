## Approximate vertex enumeration

This is the practical part of my bachelor's thesis *About a graph-theoretical algorithm for approximate vertex-enumeration* consisting of implementing an algorithm\
The thesis refers to [Andreas Löhne’s Approximate Vertex Enumeration](https://arxiv.org/abs/2007.06325), published in 2020. 
There, Löhne discusses the *approximate vertex enumeration*, which is the problem of computing a V-polytope that is close to a given H-polytope P and developes an approximate variant of Motzkin’s Double Description Method.  
For dimension three, a graph-theoretical version for approximate vertex enumeration emerges that is discussed and implemented in my thesis. 

As implementing the algorithm required the implementation of an efficient data structure for planar graphs, part of the practical part is also the implementation of the *halfedge data structure (HEDS)*.  As the implementation needs to be compared with Löhne’s original algorithm for dimension 3, the original algorithm is implemented as well.
