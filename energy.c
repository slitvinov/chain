#include <stdlib.h>
#include <stdio.h>
#include "energy.h"

#define SMALL (0.001)
#define TOLERANCE (0.05)

static real Sqrt(real);
static real Acos(real);
real
ebond0(real kb, real r0, real x0, real y0, real z0, real x1, real y1, real z1)
{
	real dx;
	real dy;
	real dz;
	real rsq;
	real th;

	dx = x0 - x1;
	dy = y0 - y1;
	dz = z0 - z1;
	rsq = dx*dx + dy * dy + dz * dz;
	th = Sqrt(rsq) - r0;
	return kb * th * th;
}


real
ebond(struct params *C, real x0, real y0, real z0, real x1, real y1, real z1)
{
	return ebond0(C->kb, C->r0, x0, y0, z0, x1, y1, z1);
}


real
eangle0(real kth, real th0, real x0, real y0, real z0, real x1, real y1, real z1, real x2, real y2, real z2)
{
	real c;
	real dtheta;
	real dx1;
	real dx2;
	real dy1;
	real dy2;
	real dz1;
	real dz2;
	real r1;
	real r2;
	real rsq1;
	real rsq2;

	dx1 = x0 - x1;
	dy1 = y0 - y1;
	dz1 = z0 - z1;
	rsq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
	r1 = Sqrt(rsq1);

	dx2 = x2 - x1;
	dy2 = y2 - y1;
	dz2 = z2 - z1;
	rsq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
	r2 = Sqrt(rsq2);

	c = dx1*dx2 + dy1*dy2 + dz1*dz2;
	c /= r1*r2;

	if (c > 1.0) c = 1.0;
	if (c < -1.0) c = -1.0;

	dtheta = Acos(c) - th0;
	return kth * dtheta * dtheta;
}

real
eangle(struct params *C, real x0, real y0, real z0, real x1, real y1, real z1, real x2, real y2, real z2)
{
	return eangle0(C->kth, C->th0, x0, y0, z0, x1, y1, z1, x2, y2, z2);
}

real
edihedral(struct params *C, real x0, real y0, real z0,
real x1, real y1, real z1,
real x2, real y2, real z2,
real x3, real y3, real z3)
{
	real ax;
	real ay;
	real az;
	real bx;
	real by;
	real bz;
	real c;
	real ra2inv;
	real rabinv;
	real rasq;
	real rb2inv;
	real rbsq;
	real vb1x;
	real vb1y;
	real vb1z;
	real vb2x;
	real vb2xm;
	real vb2y;
	real vb2ym;
	real vb2z;
	real vb2zm;
	real vb3x;
	real vb3y;
	real vb3z;

	vb1x = x0 - x1;
	vb1y = y0 - y1;
	vb1z = z0 - z1;

	vb2x = x2 - x1;
	vb2y = y2 - y1;
	vb2z = z2 - z1;

	vb2xm = -vb2x;
	vb2ym = -vb2y;
	vb2zm = -vb2z;

	vb3x = x3 - x2;
	vb3y = y3 - y2;
	vb3z = z3 - z2;
	// c,s calculation

	ax = vb1y*vb2zm - vb1z*vb2ym;
	ay = vb1z*vb2xm - vb1x*vb2zm;
	az = vb1x*vb2ym - vb1y*vb2xm;
	bx = vb3y*vb2zm - vb3z*vb2ym;
	by = vb3z*vb2xm - vb3x*vb2zm;
	bz = vb3x*vb2ym - vb3y*vb2xm;

	rasq = ax*ax + ay*ay + az*az;
	rbsq = bx*bx + by*by + bz*bz;
	ra2inv = rb2inv = 0.0;
	if (rasq > 0) ra2inv = 1.0/rasq;
	if (rbsq > 0) rb2inv = 1.0/rbsq;
	rabinv = Sqrt(ra2inv*rb2inv);
	c = (ax*bx + ay*by + az*bz)*rabinv;
	if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
		fprintf(stderr, "dihedral: tolerance\n");
		exit(2);
	}
	if (c > 1.0) c = 1.0;
	if (c < -1.0) c = -1.0;

	return C->kd * c;
}

#include <math.h>
real
Sqrt(real x)
{
	return sqrt(x);
}

real
Acos(real x)
{
	return acos(x);
}
