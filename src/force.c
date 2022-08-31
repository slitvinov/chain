#include <stdlib.h>
#include <stdio.h>
#include "energy.h"
#include "force.h"

#define SMALL (0.001)
#define SMALLER   (0.00001)
#define TOLERANCE (0.05)
#define MY_PI (3.14159265358979323846)

static real Fabs(real);
static real Acos(real);
static real Sin(real);
static real Max(real, real);
static real Sqrt(real);

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
fangle0(real kth, real th0, real x0, real y0, real z0, real x1, real y1, real z1, real x2, real y2, real z2,
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
fdihedral0(real kd, real phi0, real x0, real y0, real z0,
real x1, real y1, real z1,
real x2, real y2, real z2,
real x3, real y3, real z3,
real *fx0, real *fy0, real *fz0,
real *fx1, real *fy1, real *fz1,
real *fx2, real *fy2, real *fz2,
real *fx3, real *fy3, real *fz3
)
{
	real a;
	real a11;
	real a12;
	real a13;
	real a22;
	real a23;
	real a33;
	real b1mag;
	real b1mag2;
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
	real f1[3];
	real f2[3];
	real f3[3];
	real f4[3];
	real p;
	real pd;
	real phi;
	real r12c1;
	real r12c2;
	real rb1;
	real rb3;
	real s1;
	real s12;
	real sb1;
	real sb2;
	real sb3;
	real sc1;
	real sc2;
	real sin2;
	real sx2;
	real sy2;
	real sz2;
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
	real dx;
	real cmag;
	real s2;
	real si;
	real siinv;

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
	sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
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

	s1 = sc1 * sc1;
	s2 = sc2 * sc2;
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
	siinv = 1.0/si;

	dphi = phi - phi0;
	if (dphi > MY_PI) dphi -= 2*MY_PI;
	else if (dphi < -MY_PI) dphi += 2*MY_PI;
	p = kd*dphi;
	pd = - 2.0 * p * siinv;
	p = p * dphi;

	a = pd;
	c = c * a;
	s12 = s12 * a;
	a11 = c*sb1*s1;
	a22 = -sb2 * (2.0*c0*s12 - c*(s1+s2));
	a33 = c*sb3*s2;
	a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12);
	a13 = -rb1*rb3*s12;
	a23 = r12c2 * (c2mag*c*s2 + c1mag*s12);

	sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
	sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
	sz2  = a12*vb1z + a22*vb2z + a23*vb3z;

	f1[0] = a11*vb1x + a12*vb2x + a13*vb3x;
	f1[1] = a11*vb1y + a12*vb2y + a13*vb3y;
	f1[2] = a11*vb1z + a12*vb2z + a13*vb3z;

	f2[0] = -sx2 - f1[0];
	f2[1] = -sy2 - f1[1];
	f2[2] = -sz2 - f1[2];

	f4[0] = a13*vb1x + a23*vb2x + a33*vb3x;
	f4[1] = a13*vb1y + a23*vb2y + a33*vb3y;
	f4[2] = a13*vb1z + a23*vb2z + a33*vb3z;

	f3[0] = sx2 - f4[0];
	f3[1] = sy2 - f4[1];
	f3[2] = sz2 - f4[2];

	*fx0 += f1[0];
	*fy0 += f1[1];
	*fz0 += f1[2];

	*fx1 += f2[0];
	*fy1 += f2[1];
	*fz1+= f2[2];

	*fx2+= f3[0];
	*fy2 += f3[1];
	*fz2 += f3[2];

	*fx3 += f4[0];
	*fy3 += f4[1];
	*fz3 += f4[2];
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