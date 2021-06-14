#include <lfulib.h>

void	xbd(
	double	Vd[3],			// dipole coordinate
	double	Vs[3],			// sensor coordinate
	double	bd[3][3]		// field vector for 3 cardinal moments
)
{
	double	A[3];
	double	gradF[3];
	double	VdcrossVs[3];
	double	a2;
	double	a;
	double	VsdotA;
	double	rs2;
	double	rs;
	double	rd2;
	double	F;
	double	k0;
	double	k1;
	double	k2;
	double	kall1;
	double	kall2;
	int		v;
	int		w;

	for(v=X_, a2=VsdotA=rs2=rd2=0.; v<=Z_; v++) {
		A[v] = Vs[v] - Vd[v];
		a2 += A[v] * A[v];
		VsdotA += Vs[v] * A[v];
		rs2 += Vs[v] * Vs[v];
		rd2 += Vd[v] * Vd[v];
	}
	a = sqrt(a2);
	rs = sqrt(rs2);
	F = a * (rs * a + VsdotA);
	k0 = VsdotA / a + 2. * rs;
	k1 = a2 / rs + 2. * a + k0;
	k2 = a + k0;
	for(v=X_; v<=Z_; v++)
		gradF[v] = k1 * Vs[v] - k2 * Vd[v];

	kall1 = MU0_4PI / F;
	kall2 = kall1 / F;

	VdcrossVs[X_] = Vd[Y_] * Vs[Z_] - Vd[Z_] * Vs[Y_];
	VdcrossVs[Y_] = Vd[Z_] * Vs[X_] - Vd[X_] * Vs[Z_];
	VdcrossVs[Z_] = Vd[X_] * Vs[Y_] - Vd[Y_] * Vs[X_];
	
	for(v=0; v<3; v++)
		for(w=0; w<3; w++)
			bd[v][w] = -kall2 * VdcrossVs[v] * gradF[w];

	bd[Y_][X_] += kall1 * Vd[Z_];
	bd[Z_][X_] += -kall1 * Vd[Y_];
	bd[X_][Y_] += -kall1 * Vd[Z_];
	bd[Z_][Y_] += kall1 * Vd[X_];
	bd[X_][Z_] += kall1 * Vd[Y_];
	bd[Y_][Z_] += -kall1 * Vd[X_];
}
