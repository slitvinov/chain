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
	real h;
	real x0;
	real y0;
	real z0;
	real x1;
	real y1;
	real z1;
	
	h = 0.001;
	x0 = 1;
	y0 = 2;
	z0 = 3;
	x1 = 4;
	y1 = 8;
	z1 = 7;
	E0 = ebond0(1, 0.5, x0, y0, z0, x1, y1, z1);
	E1x = ebond0(1, 0.5, x0 + h, y0, z0, x1, y1, z1);
	E1y = ebond0(1, 0.5, x0, y0 + h, z0, x1, y1, z1);
	E1z = ebond0(1, 0.5, x0, y0, z0 + h, x1, y1, z1);
	
	printf("%.16e %.16e %.16e\n", -(E1x - E0)/h, -(E1y - E0)/h, -(E1z - E0)/h);
	fbond0(1, 0.5, x0, y0, z0, x1, y1, z1, &fx0, &fy0, &fz0, &fx1, &fy1, &fz1);
	printf("%.16e %.16e %.16e\n", fx0, fy0, fz0);
	printf("%.16e %.16e %.16e\n", fx0 + fx1, fy0 + fy1, fz0 + fz1);
}

