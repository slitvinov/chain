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

#include <math.h>
real
Sqrt(real x)
{
	return sqrt(x);
}
