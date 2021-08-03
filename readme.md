# Approximate vertex enumeration

This thesis refers to [Andreas Löhne’sApproximate Vertex Enumeration](https://arxiv.org/abs/2007.06325), published in 2020. 
There, the problem of computing a V-polytope that is close to a given H-polytope P is addressed and an approximate variant of Motzkin’s Double Description Method is developed.  
For dimension three a graphtheorical algorithm emerges fromthe correctness proof. Using this algorithm, a graph-theoretical version for approximate vertex enumerationthat relies more on the polytopes combinatorics and less on numerical computationswill be presented and implemented.  The implementation uses an edge-centered datastructure, which allows efficient manipulation of planar graphs.  Finally, the imple-mentation will be compared with an implementation of L ̈ohne’s original algorithm fordimension 3
