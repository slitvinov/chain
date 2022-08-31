#include <stdlib.h>
#include <stdio.h>
#include "energy.h"

#define SMALL (0.001)
#define SMALLER   (0.00001)
#define TOLERANCE (0.05)
#define MY_PI (3.14159265358979323846)

static real Fabs(real);
static real Acos(real);
static real Sin(real);
static real Max(real, real);
static real Sqrt(real);

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
edihedral0(real kd, real phi0, real x0, real y0, real z0,
real x1, real y1, real z1,
real x2, real y2, real z2,
real x3, real y3, real z3)
{
	real b1mag;
	real b2mag;
	real b2mag2;
	real b3mag;
	real b3mag2;
	real c;
	real c0;
	real c1mag;
	real c2mag;
	real ctmp;
	real cx;
	real cy;
	real cz;
	real dphi;
	real p;
	real phi;
	real r12c1;
	real r12c2;
	real rb1;
	real rb3;
	real s12;
	real sb1;
	real sb3;
	real sc1;
	real sc2;
	real sin2;
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
	real dx;
	real cmag;
	real si;
	real b1mag2;
	real vb1x;

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

	// c0 calculation

	sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
	sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

	rb1 = Sqrt(sb1);
	rb3 = Sqrt(sb3);

	c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

	// 1st and 2nd angle

	b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
	b1mag = Sqrt(b1mag2);
	b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
	b2mag = Sqrt(b2mag2);
	b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
	b3mag = Sqrt(b3mag2);

	ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
	r12c1 = 1.0 / (b1mag*b2mag);
	c1mag = ctmp * r12c1;

	ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
	r12c2 = 1.0 / (b2mag*b3mag);
	c2mag = ctmp * r12c2;

	// cos and sin of 2 angles and final c

	sin2 = Max(1.0 - c1mag*c1mag,0.0);
	sc1 = Sqrt(sin2);
	if (sc1 < SMALL) sc1 = SMALL;
	sc1 = 1.0/sc1;

	sin2 = Max(1.0 - c2mag*c2mag,0.0);
	sc2 = Sqrt(sin2);
	if (sc2 < SMALL) sc2 = SMALL;
	sc2 = 1.0/sc2;

	s12 = sc1 * sc2;
	c = (c0 + c1mag*c2mag) * s12;

	cx = vb1y*vb2z - vb1z*vb2y;
	cy = vb1z*vb2x - vb1x*vb2z;
	cz = vb1x*vb2y - vb1y*vb2x;
	cmag = Sqrt(cx*cx + cy*cy + cz*cz);
	dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag;

	// error check

	if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
		fprintf(stderr, "error tollerance\n");
		exit(1);
	}


	if (c > 1.0) c = 1.0;
	if (c < -1.0) c = -1.0;

	phi =Acos(c);
	if (dx > 0.0) phi *= -1.0;
	si = Sin(phi);
	if (Fabs(si) < SMALLER) si = SMALLER;

	dphi = phi - phi0;
	if (dphi > MY_PI) dphi -= 2*MY_PI;
	else if (dphi < -MY_PI) dphi += 2*MY_PI;
	p = kd*dphi;
	p = p * dphi;
	return p;
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

real
Sin(real x)
{
	return sin(x);
}

real
Max(real x, real y)
{
	return x > y ? x : y;
}

real
Fabs(real x)
{
	return x > 0 ? x : -x;
}