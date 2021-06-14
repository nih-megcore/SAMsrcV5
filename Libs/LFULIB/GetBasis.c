#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lfulib.h>
#include <samlib.h>

void GetBasis(
    double  Vc[3],
    double  scale,
    double  **gradbasis         // gradbasis[BASIS+1][3]
) {
    double  ***fabrnm;          // fabrnm[ORDER+1][ORDER+1][3]
    double  ***fabinm;          // fabinm[ORDER+1][ORDER+1][3]
    double  **harmnorm;         // harmnorm[ORDER+1][ORDER+1]
    double  *fak;               // fak[2*ORDER+1]
    double  *scalevec;          // scalevec[ORDER+1];
    double  Vs[3];              // spherical coordinates
    double  d;
    int     i, j, k;
    int     count;

    // allocations
    fabrnm = new_arrayE(double **, ORDER+1, "fabrnm[]");
    for (i=0; i<=ORDER; i++) {
        fabrnm[i] = new_matrixE(ORDER+1, 3, "fabrnm[][]");
    }
    fabinm = new_arrayE(double **, ORDER+1, "fabinm[]");
    for (i=0; i<=ORDER; i++) {
        fabinm[i] = new_matrixE(ORDER+1, 3, "fabinm[][]");
    }
    harmnorm = new_matrixE(ORDER+1, ORDER+1, "harmnorm[][]");
    fak = new_arrayE(double, 2*ORDER+1, "fak[]");
    scalevec = new_arrayE(double, ORDER+1, "scalevec[]");

    // computation of harmomic norm depends only on hull scale

    for (i=1, fak[0]=1.; i<=2*ORDER; i++)
        fak[i] = (double)i * fak[i-1];
    for (i=1, scalevec[0]=1.; i<=ORDER; i++)
        scalevec[i] = scalevec[i-1] * scale;
    for (i=0; i<=ORDER; i++)
        for (j=0; j<=i; j++) {
            d = fak[i+j] / fak[i-j] / (double)(2 * i + 1);
            if (j == 0)
                d *= 2.;
            harmnorm[i][j] = scalevec[i] * (double)i * sqrt(d);
        }

    // convert cartesian to spherical coordinates
    CtoS(Vc, Vs);

    rnynm(Vs, fabrnm, fabinm);

    for (i=1; i<=ORDER; i++)
        for (j=0; j<=i; j++) {
            for (k=0; k<3; k++) {
                fabrnm[i][j][k] /= harmnorm[i][j];
                fabinm[i][j][k] /= harmnorm[i][j];
            }
        }

    for (i=1, count=0; i<=ORDER; i++) {
        for (j=0; j<=i; j++) {
            count++;
            for (k=0; k<3; k++)
                gradbasis[count][k] = fabrnm[i][j][k];
        }
        for (j=1; j<=i; j++) {
            count++;
            for (k=0; k<3; k++)
                gradbasis[count][k] = fabinm[i][j][k];
        }
    }

    for (i=0; i<=ORDER; i++) {
        free_matrix(fabrnm[i]);
        free_matrix(fabinm[i]);
    }
    free(fabrnm);
    free(fabinm);
    free_matrix(harmnorm);
    free(fak);
    free(scalevec);
}
