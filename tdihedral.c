#include "energy.h"
#include "force.h"
#include <stdio.h>

int
main()
{
	real E0;
	real E1x;
	real E1y;
	real E1z;
	real fx0;
	real fx1;
	real fy0;
	real fy1;
	real fz0;
	real fz1;
	real fx2;
	real fy2;
	real fz2;
	real fx3;
	real fy3;
	real fz3;
	real h;
	real x0;
	real y0;
	real z0;
	real x1;
	real y1;
	real z1;
	real x2;
	real y2;
	real z2;
	real x3;
	real y3;
	real z3;

	h = 0.001;
	x0 = 1;
	y0 = 2;
	z0 = 3;
	x1 = 4;
	y1 = 9;
	z1 = 7;
	x2 = 10;
	y2 = 11;
	z2 = 14;
	x3 = 18;
	y3 = 19;
	z3 = 22;

	E0 = edihedral0(0.5, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

	E1x = edihedral0(0.5, x0 + h, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
	E1y = edihedral0(0.5, x0, y0 + h, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
	E1z = edihedral0(0.5, x0, y0, z0 + h, x1, y1, z1, x2, y2, z2, x3, y3, z3);
	printf("%.16e %.16e %.16e\n", -(E1x - E0)/h, -(E1y - E0)/h, -(E1z - E0)/h);
	fdihedral0(0.5, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, &fx0, &fy0, &fz0, &fx1, &fy1, &fz1, &fx2, &fy2, &fz2, &fx3, &fy3, &fz3);
	printf("%.16e %.16e %.16e\n", fx0, fy0, fz0);
	printf("\n");

	E1x = edihedral0(0.5, x0, y0, z0, x1 + h, y1, z1, x2, y2, z2, x3, y3, z3);
	E1y = edihedral0(0.5, x0, y0, z0, x1, y1 + h, z1, x2, y2, z2, x3, y3, z3);
	E1z = edihedral0(0.5, x0, y0, z0, x1, y1, z1 + h, x2, y2, z2, x3, y3, z3);
	printf("%.16e %.16e %.16e\n", -(E1x - E0)/h, -(E1y - E0)/h, -(E1z - E0)/h);
	printf("%.16e %.16e %.16e\n", fx1, fy1, fz1);
	printf("\n");

	E1x = edihedral0(0.5, x0, y0, z0, x1, y1, z1, x2 + h, y2, z2, x3, y3, z3);
	E1y = edihedral0(0.5, x0, y0, z0, x1, y1, z1, x2, y2 + h, z2, x3, y3, z3);
	E1z = edihedral0(0.5, x0, y0, z0, x1, y1, z1, x2, y2, z2 + h, x3, y3, z3);
	printf("%.16e %.16e %.16e\n", -(E1x - E0)/h, -(E1y - E0)/h, -(E1z - E0)/h);
	printf("%.16e %.16e %.16e\n", fx2, fy2, fz2);

	printf("%.16e %.16e %.16e\n", fx0 + fx1 + fx2 + fx3, fy0 + fy1 + fy2 + fy3, fz0 + fz1 + fz2 + fz3);
}
