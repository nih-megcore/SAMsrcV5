// FIELD's will be of the following types:
//
//			position		orientation
//			========		===========
//	Cartesian:	X, Y, Z			X, Y, Z
//
//	Spherical:	RHO, THETA, PHI		RHO, THETA, PHI
//
//	5-Parameter Q:	RHO, THETA, PHI		PSI, |Q|
//
//	Euler:		X, Y, Z			PSI, THETA, PHI
//
//
//	All parameters will be passed in fundamental SI units
//	such as meters, Ampere-meters, Tesla, volts, &
//	all angles will be in radians
//
//	Author: SE Robinson
//			MEG Core Group
//			NIMH
//

#ifndef H_GEOMS
#define H_GEOMS

// FIELD indices
#define X_			0
#define Y_			1
#define Z_			2

#define RHO_		0				// radial displacement
#define THETA_		1				// "latitude" -- 0 to pi
#define PHI_		2				// "longitude" -- 0 to2 pi

#define PSI_		0
#define Q_			1

// MRI indices
#define S_			0				// sagittal
#define	C_			1				// coronal
#define A_			2				// axial
#define NA_			0				// nasion
#define LE_			1				// left ear
#define RE_			2				// right ear

// sensor types
#define MEG			0
#define EEG			1

// logicals
#define TRUE		1
#define FALSE		0
#define GOOD		0
#define BAD			1

// model parameter
#define NOSPHERE	0				// infinite homogeneous conductor
#define SPHERE		1				// dipole in sphere
#define MESH		2				// surface integral
#define MDIPOLE		3				// magnetic dipole
#define MULTIPOLE	4				// current multipole expansion

// miscellaneous
#define MINUTE		1.0e-300		// minute, as in small - not minute, as in 60 seconds!
#define INTPTS		12				// maximum field evaluation points
#define	SCALP		0				// head surface type
#define CORTICAL	1				// hull surface type
#define INNERSKULL	2				// inner skull surface type

// 'FIELD' -- position & orientation vector (either Cartesian or spherical)
typedef struct {
	double		p[3];				// position vector of field point
	double		v[3];				// orientation vector of field point
} FIELD;

// 'OPMFIELD' -- position & spherical orientation vector with gain 
typedef struct {
	double		p[3];				// position vector of field point
	double		v[2];				// orientation angles of field point
	double		g;				// gain of sensor
} OPMFIELD;

// 'XFORM' -- coordinate transformation data
typedef struct {
	double		Translate[3];		// The translation or local origin
	double		Rotate[3][3];		// 3 x 3 rotation matrix
} XFORM;

// 'QUATERNION' -- quaternion components
typedef struct {
	double		x;
	double		y;
	double		z;
	double		w;
} QUATERNION;

// 'COIL' -- the PHYSICAL coil data
typedef struct {
	FIELD		origin;				// coil origin vectors
	double		radius;				// coil radius
	int			SenseTurns;			// sense & turns (close-spaced)
	FIELD		B[INTPTS];			// B evaluation points -- up to 12 point integration
	double		w[INTPTS];			// integration weights
	int			NumInt;				// number of integration points
	int			Evaluated;			// evaluation flag & number of integration points (0 if not evaluated)
} COIL;

// 'SENSOR' -- the PHYSICAL MEG sensor data
typedef struct {
	FIELD		origin;				// sensor origin vectors
	COIL		*Coil;				// Coil[] -- coil array for sensor
	int			NumCoils;			// number of coils in sensor
	double		LocalSphere[3];		// local sphere origin for sensor
} SENSOR;

// 'OPM' -- the PHYSICAL OPM sensor data
typedef struct {
    FIELD       origin;             // sensor origin vector
    double      LocalSphere[3];     // local sphere origin for sensor
} OPM;

// 'ELECTRODE' -- the PHYSICAL EEG electrode data
typedef struct {
	double		position[3];		// electrode position in headframe
	double		gain;				// sensitivity (Volts/Volt)
	short		RefNumber;			// reference channel number (from a/d list)
} ELECTRODE;

// 'PROBE' -- the PHYSICAL probe data
typedef struct {
	XFORM		origin;				// probe transformation vectors
	SENSOR		*Sensor;			// Sensor[] -- sensor array for probe
	int			NumSensors;			// number sensors in probe
	int			Xformed;			// transformation flag
} PROBE;

// 'SNODE' -- scan node - a linked list of scan information
typedef struct cnode {
	SENSOR		*Sensor;			// sensor structure
	double		Lsp[3];				// local sphere origin
	double		*Xm;				// Xm[] -- measured X at time
	double		*Xp;				// Xp[] -- predicted X at time
	double		*Xs;				// Xs[] -- standard deviation of X
	int			T;					// number of time samples
	struct cnode	*next;			// link to next scan point
} SNODE;

// 'BNODE' -- measurement node
typedef struct {
	SENSOR		*Sensor;			// sensor structure
	double		Xm;					// measured field
	double		Xp;					// predicted field
	double		Xn;					// measurement noise variance
	double		Lsp[3];				// local sphere origin
} BNODE;

// 'TGLS' -- triangles mapping surface of conductor
typedef struct	tgls {
	FIELD		c;					// center of triangle element
	double		v1[3];				// coordinate of vertex-1
	double		v2[3];				// coordinate of vertex-2
	double		v3[3];				// coordinate of vertex-3
	double		ds;					// surface area of triangle
	struct		tgls	*next;		// pointer to next element
} TGLS;

// 'SVD' -- structure holding singular value decomposition sub-matrices
typedef struct	{
	double		**U;				// M by N input/ouput matrix
	double		*W;					// singular values
	double		**V;				// N by N orthogonal matrix
	int			M;					// row dimension
	int			N;					// column dimension
	double		noise;				// noise variance (least-significant singular value)
	double		span;				// span of singular values ('condition number')  
} SVD;

// 'SVDL' -- structure holding singular value decomposition sub-matrices
typedef struct	{
	long double	*U;					// M*N input/ouput matrix
	long double	*W;					// N singular values
	long double	*V;					// N*N orthogonal matrix
	long double	*rv1;				// workspace
	int			M;					// row dimension
	int			N;					// column dimension
	double		noise;				// noise variance (least-significant singular value)
	double		span;				// span of singular values ('condition number')  
} SVDL;

// 'PROTO' -- prototype sensor
typedef struct {
	int			NumCoils;			// number of coils in sensor
	double		Offset[4];			// Offset[] -- coil baseline array
	double		Radius[4];			// Radius[] -- coil radius array
	int			SenseTurns[4];		// SenseTurns[] -- sense-turn array
	double		gain;				// sensor gain (T/V)
} PROTO;

// 'MODEL' -- model structure
typedef struct {
	int			SensorType;			// magnetic or electric
	int			model;				// model type
	double		LocalSphere[3];		// local sphere origin
	double		sigma;				// conductivity (S)
	TGLS		surface;			// surface elements
} MODEL;

// 'DIGPTS' -- headshape digitization solution structure
typedef struct {
	double		Observed[3];		// observed digitization point
	double		Predicted[3];		// Predicted digitization point
	short		Flag;				// flag points to be used
} DIGPTS;

// 'CART' -- cartesian vector
typedef struct {
	double		x1;
	double		x2;
	double		x3;
} CART;

// 'MULTISPHERE' -- structure for multiple local sphere coordinates
typedef struct {
	char		ChanName[32];		// sensor name - 32 characters, for use in weight files
	double		Origin[3];			// sphere origin for this sensor (m)
	double		Radius;				// sphere radius for this sensor (m)
} MULTISPHERE;

// 'HEAD_MODEL' -- structure for '*.hdm' file
typedef struct {
	char		MriName[256];		// MRI file name (or path)
	double		SphereOrigin[3];	// mean or single-sphere origin (m)
	double		SphereRadius;		// mean or single-sphere radius (m)
	double		MriRes[3];			// MRI voxel resolution
	int			Fiducial[3][3];		// fiducial coordinates (voxels)
	double		SearchRadius;		// multi-sphere search radius (m)
	char		HeadShapeName[256];	// head-shape file name
	int			SurfaceType;		// head or hull
	int			NumSensors;			// number of sensors (reference + primary)
	MULTISPHERE	MultiSphere[500];	// lots of room for bigger systems
} HEAD_MODEL;

// 'SHAPEINFO' -- structure for list of head-shape points
typedef struct {					// 'SHAPEINFO' -- headshape digitization structure
	double		Observed[3];		// observed digitization point
	int			flag;				// flag points to be used
} SHAPEINFO;

// 'VOXINFO' -- structure to hold voxel array
typedef struct {
	double		p[3];				// voxel position
	int			flag;				// voxel flag
} VOXINFO;

#endif	// H_GEOMS
