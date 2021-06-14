// Swaps.h -- collection of byte swapping routines for short int (int16),
//	long int (int32), & float (float32).
//
//	Author:	Stephen E. Robinson
//			Neuromagnetism Laboratory
//			Dept.of Neurology
//			Henry Ford Hospital
//			Detroit, MI
//			Copyright (c) 2009
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


// swap endian for short int
static union {
	short int		s;
	unsigned char	b[2];
} sdat1, sdat2;
#define SHORTSWAP(sin) (sdat1.s=sin,sdat2.b[0]=sdat1.b[1],sdat2.b[1]=sdat1.b[0],sdat2.s)


// swap endian for long int
static union {
	long int		l;
	unsigned char	b[4];
} ldat1, ldat2;
#define LONGSWAP(lin) (ldat1.l=lin,ldat2.b[0]=ldat1.b[3],ldat2.b[1]=ldat1.b[2],ldat2.b[2]=ldat1.b[1],ldat2.b[3]=ldat1.b[0],ldat2.l)


// swap endian for float
static union {
	float			f;
	unsigned char	b[4];
} fdat1, fdat2;
#define FLOATSWAP(fin) (fdat1.f=fin,fdat2.b[0]=fdat1.b[3],fdat2.b[1]=fdat1.b[2],fdat2.b[2]=fdat1.b[1],fdat2.b[3]=fdat1.b[0],fdat2.f)