#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <lfulib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

int ORDER;                              // spherical harmonic expansion order

void GetCoeffs(
    FIELD           *detector,          // detector[S] -- array of detection points
    int             S,                  // number of detection points
    HULL            *Hull,              // conducting volume vertices & normals
    double          **coeffs            // coeffs[BASIS][S] -- compensation coefficients for each field integration point
) {
    gsl_matrix      *CtC;               // C'C with index 0
    gsl_matrix      *CtB;               // C'b with index 0
    gsl_permutation *indx;              // LU decomposition permutations
    gsl_vector      *achan;             // row of C'b & solution vector
    double          **gradbasis;        // gradbasis[BASIS+1][3] -- gradient of basis functions
    double          *gradbasisn;        // gradbasisn[BASIS] -- gradients
    double          DipLoc[3];          // vertex location in hull origin fram
    double          SensLoc[3];         // sensor location in hull origin frame
    double          DipOri[3];          // vertex normal orienation
    double          SensOri[3];         // sensor normal orientation
    double          Bd[3][3];           // solution for dipole in spherical conductor
    double          x, tmp;
    int             i, j, k;
    int             v, w;
    int             signum;
    int             BASIS;              // number of basis functions

    BASIS = ((ORDER+1) * (ORDER+1) - 1);
    gradbasis = new_matrixE(BASIS + 1, 3, "gradbasis[][]");
    gradbasisn = new_arrayE(double, BASIS, "gradbasisn[]");

    // allocate gsl matices & vectors -- initializing C'C & C'b to zero
    CtC = gsl_matrix_alloc(BASIS, BASIS);
    CtB = gsl_matrix_alloc(BASIS, S);
    indx = gsl_permutation_alloc(BASIS);
    achan = gsl_vector_alloc(BASIS);

    // accumulate C'C & C'b over all conducting volume vertices
    gsl_matrix_set_zero(CtC);
    gsl_matrix_set_zero(CtB);
    for (i=0; i<Hull->nv; i++) {

        // translate vertex location to hull origin
        // - use inflated hull
        for (v=X_; v<=Z_; v++) {
            DipLoc[v] = Hull->infvertex[i].p[v] - Hull->Vo[v];
            DipOri[v] = Hull->infvertex[i].v[v];
        }

        // compute basis function for vertex of hull
        GetBasis(DipLoc, Hull->scale, gradbasis);
        for (j=0; j<BASIS; j++)
            for (v=X_, gradbasisn[j]=0.; v<=Z_; v++)    // gradbasisn = gradbasis <dot> vertex orientation
                gradbasisn[j] += gradbasis[j+1][v] * DipOri[v];

        // accumulate C'C upper triangle
        for (j=0; j<BASIS; j++)
            for (k=0; k<=j; k++) {
                tmp = gsl_matrix_get(CtC, j, k);
                tmp += gradbasisn[j] * gradbasisn[k];
                gsl_matrix_set(CtC, j, k, tmp);
            }

        // accumulate C'b over all detection locations
        for (j=0; j<S; j++) {

            // translate detector location to hull origin
            for (v=X_; v<=Z_; v++) {
                SensLoc[v] = detector[j].p[v] - Hull->Vo[v];
                SensOri[v] = detector[j].v[v];
            }

            // compute forward solution for dipole in spherical conductor
            xbd(DipLoc, SensLoc, Bd);                   // solve Sarvas forward solution for vertex & detector locations
            for (v=X_, x=0.; v<=Z_; v++)                // x is dot product of vertex normal field with detector orientation
                for (w=X_; w<=Z_; w++)
                    x += Bd[v][w] * DipOri[v] * SensOri[w];

            // C'b
            for (k=0; k<BASIS; k++) {
                tmp = gsl_matrix_get(CtB, k, j);
                tmp +=  x * gradbasisn[k];
                gsl_matrix_set(CtB, k, j, tmp);
            }
        }
    }

    // mirror C'C to lower triangle
    for (j=0; j<BASIS; j++)
        for (k=j+1; k<BASIS; k++)
            gsl_matrix_set(CtC, j, k, gsl_matrix_get(CtC, k, j));

    // compute coefficients a = [C'C]^-1 [C'b]
    if (gsl_linalg_LU_decomp(CtC, indx, &signum) != GSL_SUCCESS)
        cleanup("gsl_linalg_LU_decomp() failed");
    for (j=0; j<S; j++) {
        for (i=0; i<BASIS; i++)
            gsl_vector_set(achan, i, gsl_matrix_get(CtB, i, j));
        gsl_linalg_LU_svx(CtC, indx, achan);        // achan replaced by solution
        for (i=0; i<BASIS; i++)
            coeffs[i][j] = gsl_vector_get(achan, i);
    }

    // free gsl matrices & vectors
    gsl_matrix_free(CtC);
    gsl_matrix_free(CtB);
    gsl_permutation_free(indx);
    gsl_vector_free(achan);

    free_matrix(gradbasis);
    free(gradbasisn);
}
