## Approximate vertex enumeration

This is the practical part of my bachelor's thesis *About a graph-theoretical algorithm for approximate vertex-enumeration*.\
The thesis refers to [Andreas Löhne’s Approximate Vertex Enumeration](https://arxiv.org/abs/2007.06325), published in 2020. 
There, the *approximate vertex enumeration*, which is the problem of computing a V-polytope that is close to a given H-polytope P is discussed and an approximate variant of Motzkin’s Double Description Method is developed.  
For dimension three, a graph-theoretical version for approximate vertex enumeration emerges.

The implementation uses *the halfedge data structure (HEDS)*, which allows efficient manipulation of planar graphs.  Finally, the imple-mentation will be compared with an implementation of L ̈ohne’s original algorithm fordimension 3
