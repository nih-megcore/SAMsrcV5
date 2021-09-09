// gsl_diag() -- routines for adding or scaling the diagonal of a gsl matrix
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <lfulib.h>
#include <samlib.h>

// add factor to diagonal
void    gsl_matrix_add_diag(
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

// multiply diagonal by factor
void    gsl_matrix_mul_diag(
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

void    gsl_matrix_plus_equal(
    gsl_matrix  *MyMatrix,
    int         i,
    int         j,
    double      factor
)
{
    gsl_matrix_set(MyMatrix, i, j, gsl_matrix_get(MyMatrix, i, j) + factor);
}

void    gsl_vector_plus_equal(
    gsl_vector  *MyVector,
    int         i,
    double      factor
)
{
    gsl_vector_set(MyVector, i, gsl_vector_get(MyVector, i) + factor);
}
