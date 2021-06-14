#define ERROR(msg) exit (fprintf (stderr, "%s\n", msg))
#define LEFT	0
#define RIGHT	1

// matrix functions
void	copy_matrix(double * a, double * b, int n);
void	identity_matrix(double ** m, int n);
void	transpose(double *matrix, int n);
void	multiply(double *m1, double *m2, double * new_matrix, int n);
void	print_matrix(double *m, int rows, int columns);
double	*gen_matrix(int n);
void	swap_rows(double * i_mat, int row1, int row2, int n);


// Jacobi method functions
void	jacobi(double * a, int n, double * s, double * u, double * v);
void	generate_composite_matrix(double ** matrix, int p, int q, int n, double angle, int side);
void	rotate(double * a, int n, double * u, double * v, int p, int q);
int		not_converged(double * a, int n);
void	order_a(double * a, int n, double * u, double * v);
void	generate_s(double * a, double * s, int n);
