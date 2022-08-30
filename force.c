#include <stdlib.h>
#include <stdio.h>
#include "energy.h"
#include "force.h"

#define SMALL (0.001)
#define TOLERANCE (0.05)

static real Sqrt(real);
static real Acos(real);
void
fbond0(real kb, real r0, real x0, real y0, real z0, real x1, real y1, real z1, real *fx0, real *fy0, real *fz0, real *fx1, real *fy1, real *fz1)
{
	real dr;
	real dx;
	real dy;
	real dz;
	real fbond;
	real r;
	real rk;
	real rsq;

	dx = x0 - x1;
	dy = y0 - y1;
	dz = z0 - z1;
	rsq = dx*dx + dy * dy + dz * dz;
	r = Sqrt(rsq);
	dr = r - r0;
	rk = kb * dr;
	fbond = r > 0 ? -2 * rk / r : 0;
	*fx0 += dx * fbond;
	*fy0 += dy * fbond;
	*fz0 += dz * fbond;
	*fx1 -= dx * fbond;
	*fy1 -= dy * fbond;
	*fz1 -= dz * fbond;
}

void
fangle(real kth, real th0, real x0, real y0, real z0, real x1, real y1, real z1, real x2, real y2, real z2,
real *fx0, real *fy0, real *fz0,
real *fx1, real *fy1, real *fz1,
real *fx2, real *fy2, real *fz2
)
{
	real a;
	real a11;
	real a12;
	real a22;
	real c;
	real dtheta;
	real dx1;
	real dx2;
	real dy1;
	real dy2;
	real dz1;
	real dz2;
	real f1[3];
	real f3[2];
	real r1;
	real r2;
	real rsq1;
	real rsq2;
	real s;
	real tk;

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

	s = Sqrt(1.0 - c * c);
	if (s < SMALL) s = SMALL;
	s = 1.0 / s;

	dtheta = Acos(c) - th0;
	tk = kth * dtheta;

	a = -2.0 * tk * s;
	a11 = a * c / rsq1;
	a12 = -a / (r1 * r2);
	a22 = a * c / rsq2;

	f1[0] = a11 * dx1 + a12 * dx2;
	f1[1] = a11 * dy1 + a12 * dy2;
	f1[2] = a11 * dz1 + a12 * dz2;
	f3[0] = a22 * dx2 + a12 * dx1;
	f3[1] = a22 * dy2 + a12 * dy1;
	f3[2] = a22 * dz2 + a12 * dz1;

	*fx0 += f1[0];
	*fy0 += f1[1];
	*fz0  += f1[2];

	*fx1 -= f1[0] + f3[0];
	*fy1 -= f1[1] + f3[1];
	*fz1-= f1[2] + f3[2];

	*fx2 += f3[0];
	*fy2 += f3[1];
	*fz2 += f3[2];
}

void
fdihedral(real kd, real x0, real y0, real z0,
real x1, real y1, real z1,
real x2, real y2, real z2,
real x3, real y3, real z3,
real *fx0, real *fy0, real *fz0,
real *fx1, real *fy1, real *fz1,
real *fx2, real *fy2, real *fz2,
real *fx3, real *fy3, real *fz3
)
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
