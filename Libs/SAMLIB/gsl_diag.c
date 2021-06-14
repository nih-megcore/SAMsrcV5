// gsl_diag() -- routines for adding or scaling the diagonal of a gsl matrix
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <lfulib.h>
#include <samlib.h>

void	gsl_matrix_add_diag(
	gsl_matrix	*MyMatrix,
	double		factor
)
{
	int i;
	int	dim;

	dim = MyMatrix->size1;
	if(MyMatrix->size2 != dim)
		cleanup("matrix must be square for diagonal operation");
	for(i=0; i<dim; i++)
		gsl_matrix_set(MyMatrix, i, i, gsl_matrix_get(MyMatrix, i, i) + factor);
}

void	gsl_matrix_mul_diag(
	gsl_matrix	*MyMatrix,
	double		factor
)
{
	int i;
	int	dim;

	dim = MyMatrix->size1;
	if(MyMatrix->size2 != dim)
		cleanup("matrix must be square for diagonal operation");
	for(i=0; i<dim; i++)
		gsl_matrix_set(MyMatrix, i, i, gsl_matrix_get(MyMatrix, i, i) * factor);
}
