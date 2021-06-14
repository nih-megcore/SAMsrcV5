// definition of cortical mesh elements for coregistration of MEG with MRI

#ifndef H_MESH
#define H_MESH

typedef struct {
    double      Vp[3];      // mesh vertex position
    double      Vu[3];      // weakly-constrained source moment vector
    double      Vc[3];      // strongly-constrained source moment vector
    double      z2;         // z^2 source metric (from strong constraints)
    double      angle;      // acos|dot| product of constrained & unconstrained moment vectors
    double      span;       // eigenvalue span of moment vector solution
} MESHINFO;

#endif  // H_MESH
