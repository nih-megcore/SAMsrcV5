#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lfulib.h>

void rnynm(
    double      Vs[3],
    double      ***fabrnm,      // fabrnm[ORDER+1][ORDER+1][3]
    double      ***fabinm       // fabinm[ORDER+1][ORDER+1][3]
) {
    double      **frnm;         // frnm[ORDER+1][ORDER+1]
    double      **finm;         // finm[ORDER+1][ORDER+1]
    double      **pnm;          // pnm[ORDER+1][ORDER+1]
    double      *ff;            // ff[2*(ORDER+1)+2]
    double      *sinn;          // sinn[ORDER+1]
    double      *cmphi;         // cmphi[ORDER+1]
    double      *smphi;         // smphi[ORDER+1]
    double      *rn;            // rn[ORDER+1];
    double      s, c;
    double      cphi;
    double      sphi;
    int         i, j, k;

    pnm = new_matrixE(ORDER+1, ORDER+1, "pnm[]");
    ff = new_arrayE(double, 2*(ORDER+1)+2, "ff[]");
    sinn = new_arrayE(double, ORDER+1, "sinn[]");
    cmphi = new_arrayE(double, ORDER+1, "cmphi[]");
    smphi = new_arrayE(double, ORDER+1, "smphi[]");
    rn = new_arrayE(double, ORDER+1, "rn[]");
    frnm = new_matrixE(ORDER+1, ORDER+1, "frnm[][]");
    finm = new_matrixE(ORDER+1, ORDER+1, "finm[][]");

    for (i=0; i<=ORDER; i++)
        for (j=0; j<=ORDER; j++) {
            frnm[i][j] = 0.;
            finm[i][j] = 0.;
            for (k=0; k<3; k++) {
                fabrnm[i][j][k] = 0.;
                fabinm[i][j][k] = 0.;
            }
        }

    c = cos(Vs[THETA_]);
    s = sin(Vs[THETA_]);

    for(i=1, ff[1]=1.; i<=ORDER; i++)
        ff[2*i+1] = ff[2*i-1] * (double)(2 * i + 1);
    for(i=2, sinn[1]=s; i<=ORDER; i++)
        sinn[i] = sinn[i-1] * s;

    for(i=1, pnm[0][0]=1.;i<=ORDER; i++)
        pnm[i][i] = ff[2*i-1] * sinn[i];
    for(i=0; i<=ORDER-1; i++)
        pnm[i+1][i] = c * (double)(2 * i + 1) * pnm[i][i];

    for(i=0;i <=ORDER; i++)
        for(j=2; j<=(ORDER-i); j++)
            pnm[i+j][i] = (c * (double)(2 * (i + j) - 1) * pnm[i+j-1][i] - (double)(2 * i + j - 1) * pnm[i+j-2][i]) / (double)j;

    cphi = cos(Vs[PHI_]);
    sphi = sin(Vs[PHI_]);
    cmphi[0] = 1.;
    smphi[0] = 0.;
    cmphi[1] = cphi;
    smphi[1] = sphi;
    for(i=2; i<=ORDER; i++) {
        cmphi[i] = cmphi[i-1] * cphi - smphi[i-1] * sphi;
        smphi[i] = cmphi[i-1] * sphi + smphi[i-1] * cphi;
    }

    for(i=1, rn[0]=1.;i<=ORDER; i++)
        rn[i] = rn[i-1] * Vs[RHO_];

    for(i=0; i<=ORDER; i++)
        for(j=0; j<=i; j++) {
            frnm[i][j] = pnm[i][j] * cmphi[j] * rn[i];
            finm[i][j] = pnm[i][j] * smphi[j] * rn[i];
        }

    for(i=0; i<=1; i++)
        for(j=0; j<=i; j++)
            for(k=0; k<3; k++) {
                fabrnm[i][j][k] = 0.;
                fabinm[i][j][k] = 0.;
            }

    for(i=1; i<=ORDER; i++)
        for(j=0; j<=i; j++) {
            fabrnm[i][j][2] = (double)(i + j) * frnm[i-1][j];
            fabinm[i][j][2] = (double)(i + j) * finm[i-1][j];
        }

    for(i=2; i<=ORDER; i++) {
        fabrnm[i][0][0] = -frnm[i-1][1];
        fabinm[i][0][0] = 0.;
        fabrnm[i][0][1] = -finm[i-1][1];
        fabinm[i][0][1] = 0.;
    }

    for(i=1; i<=ORDER; i++)
        for(j=1; j<=i; j++) {
            fabrnm[i][j][0] = 0.5 * (double)(i + j - 1) * (double)(i + j) * frnm[i-1][j-1];
            fabinm[i][j][0] = 0.5 * (double)(i + j - 1) * (double)(i + j) * finm[i-1][j-1];
            fabrnm[i][j][1] = -0.5 * (double)(i + j - 1) * (double)(i + j) * finm[i-1][j-1];
            fabinm[i][j][1] = 0.5 * (double)(i + j - 1) * (i + j) * frnm[i-1][j-1];
        }

    for(i=1; i<=ORDER; i++)
        for(j=1; j<=i-2; j++) {
            fabrnm[i][j][0] += -0.5 * frnm[i-1][j+1];
            fabinm[i][j][0] += -0.5 * finm[i-1][j+1];
            fabrnm[i][j][1] += -0.5 * finm[i-1][j+1];
            fabinm[i][j][1] += 0.5 * frnm[i-1][j+1];
        }

    free_matrix(pnm);
    free(ff);
    free(sinn);
    free(cmphi);
    free(smphi);
    free(rn);
    free_matrix(frnm);
    free_matrix(finm);
}
