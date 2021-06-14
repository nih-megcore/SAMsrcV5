// LUsolve() -- solve Ax = b using LU decomposition
//	This is call-compatible GSL eplacement GJ.c
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <stdlib.h>
#include <math.h>
#include <lfulib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

void	LUsolve(
	double	**a,	// a[N][N] (on input) & its inverse on output)
	double	*da,	// input & solution vector
	int		N		// matrix rank
)
{
	gsl_matrix		*A;					// gsl matrix copy of a
	gsl_matrix		*Ainv;				// inverse A
	gsl_permutation	*indx;				// LU decomposition permutations
	gsl_vector		*Da;				// input vector
	gsl_vector		*x;					// solution vector
	int				i;
	int				j;
	int				signum;

	// allocate gsl matrices & vectors	
	A = gsl_matrix_alloc(N, N);
	Ainv = gsl_matrix_alloc(N, N);
	indx = gsl_permutation_alloc(N);
	Da = gsl_vector_alloc(N);
	x = gsl_vector_alloc(N);

	// copy a to gsl A
	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			gsl_matrix_set(A, i, j, a[i][j]);
	
	// copy da to gsl Da
	for(i=0; i<N; i++)
		gsl_vector_set(Da, i, da[i]);

	// solve Ax = b
	gsl_linalg_LU_decomp(A, indx, &signum);
	gsl_linalg_LU_solve(A, indx, x, Da);
	for(i=0; i<N; i++)
		da[i] = gsl_vector_get(x, i);

	// solve A^-1
	gsl_linalg_LU_invert(A, indx, Ainv);

	// copy A^1 to a
	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			a[i][i] = gsl_matrix_get(Ainv, i, j);
	
	// free gsl matrices & vectors
	gsl_matrix_free(A);
	gsl_matrix_free(Ainv);
	gsl_vector_free(Da);
	gsl_vector_free(x);
	gsl_permutation_free(indx);
}
