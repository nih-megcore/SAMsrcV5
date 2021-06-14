// geig22() -- generalized eigensystem problem for 2x2 matrices
//  Exact solution from Mathematica
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <math.h>
#include <lfulib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


void	geig22(
    gsl_matrix  *A,         // 2x2 matrix A
    gsl_matrix  *B,         // 2x2 matrix B
    gsl_matrix  *eVec,      // 2x2 eigenvector matrix
	gsl_vector  *eVal       // 2x1 eigenvalues
)
{
	double	a11, a11_2;
	double	a12, a12_2;
	double	a22, a22_2;
	double	b11, b11_2;
	double	b12, b12_2;
	double	b22, b22_2;
	double	u1, u1_2;
	double	u2, u2_2;
	double	u3;
	double	u4;
	double	u5;
	double	u6;
	double	u7;
	double	u8;
	double	u9;
	double	u10;
	double	l1;
	double	l2;
	double	k1;
	double	k2;

    a11 = gsl_matrix_get(A, 0, 0);
    a12 = gsl_matrix_get(A, 0, 1);
    a22 = gsl_matrix_get(A, 1, 1);
    b11 = gsl_matrix_get(B, 0, 0);
    b12 = gsl_matrix_get(B, 0, 1);
    b22 = gsl_matrix_get(B, 1, 1);

	u1 = a12 * b11 - a11 * b12;
	u2 = a22 * b11 - 2. * a12 * b12 + a11 * b22;
	a11_2 = a11 * a11;
	a12_2 = a12 * a12;
	a22_2 = a22 * a22;
	b11_2 = b11 * b11;
	b12_2 = b12 * b12;
	b22_2 = b22 * b22;
	u1_2 = u1 * u1;
	u2_2 = u2 * u2;
	u3 = sqrt(-4. * (a12_2 - a11 * a22) * (b12_2 - b11 * b22) + u2_2);
	u4 = a22_2 * b11_2
		+ 2. * a12_2 * b11 * (b11 + b22)
		- 2. * a12 * b12 * (a22 * b11 + a11 *(2. * b11 + b22))
		+ a22 * (2. * a11 * b12_2 - 2. * a11 * b11 * b22 + b11 * u3)
		+ a11 * (2. * a11 * b12_2 + a11 * b22_2 - b22 * u3);
	u5 = a22_2 * b11_2
		+ 2. * a11_2 * b12_2
		+ 2. * a12_2 * b11 * (b11 + b22) - 2. * a12 * b12 * (a22 * b11 + a11 * (2. * b11 + b22))
		+ a11 * b22 *(a11 * b22 + u3)
		+ a22 * (2. * a11 * (b12_2 - b11 * b22) - b11 * u3);
	u6 = a22 * (-2. * b12_2 + b11 * b22) + b22 * (2. * a12 * b12 - a11 * b22 + u3);
	u7 = a22 * (-2. * b12_2 + b11 * b22) + b22 * (2. * a12 * b12 - a11 * b22 - u3);
	u8 = a22 * b11 * b12 - 2. * a12 * b11 * b22 + b12 * (a11 * b22 - u3);
	u9 = a22 * b11 * b12 - 2. * a12 * b11 * b22 + b12 * (a11 * b22 + u3);
	u10 = 2. * (b12_2 - b11 * b22);
	l1 = (u3 - u2) / u10;
	l2 = (-u3 - u2) / u10;
	k1 = M_SQRT2 * sqrt(u1_2 / u4);
	k2 = M_SQRT2 * sqrt(u1_2 / u5);
	if(l1 > l2) {
		gsl_vector_set(eVal, 0, l1);
		gsl_vector_set(eVal, 1, l2);
		gsl_matrix_set(eVec, 0, 0, k1 * (u6 / u8));
		gsl_matrix_set(eVec, 0, 1, k1);
		gsl_matrix_set(eVec, 1, 0, k2 * (u7 / u9));
		gsl_matrix_set(eVec, 1, 1, k2);
	} else {
		gsl_vector_set(eVal, 0, l2);
		gsl_vector_set(eVal, 1, l1);
		gsl_matrix_set(eVec, 0, 0, k2 * (u7 / u9));
		gsl_matrix_set(eVec, 0, 1, k2);
		gsl_matrix_set(eVec, 1, 0, k1 * (u6 / u8));
		gsl_matrix_set(eVec, 1, 1, k1);
	}
}
