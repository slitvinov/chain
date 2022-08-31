#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "energy.h"
#include "force.h"

static void
force(struct params *C, int n, real *r, real *f)
{
	int i;
	int j;
	int k;
	int l;

	for (i = 0; i  < n - 1; i++) {
		j = i + 1;
		fbond0(C->kb, C->r0, r[3*i], r[3*i + 1], r[3*i + 2], r[3*j], r[3*j + 1], r[3*j + 2],
		    &f[3*i], &f[3*i + 1], &f[3*i + 2], &f[3*j], &f[3*j + 1], &f[3*j + 2]);
	}

	for (i = 0; i < n - 2; i++) {
		j = i + 1;
		k = i + 2;
		fangle0(C->kth, C->th0,  r[3*i], r[3*i + 1], r[3*i + 2],  r[3*j], r[3*j + 1], r[3*j + 2],  r[3*k], r[3*k + 1], r[3*k + 2],
		    &f[3*i], &f[3*i + 1], &f[3*i + 2], &f[3*j], &f[3*j + 1], &f[3*j + 2], &f[3*k], &f[3*k + 1], &f[3*k + 2]);
	}

	for (i = 0; i < n - 3; i++) {
		j = i + 1;
		k = i + 2;
		l = i + 3;
		fdihedral0(C->kd, C->phi0,
		    r[3*i], r[3*i + 1], r[3*i + 2], 
		    r[3*j], r[3*j + 1], r[3*j + 2], 
		    r[3*k], r[3*k + 1], r[3*k + 2],
		    r[3*l], r[3*l + 1], r[3*l + 2],
		    &f[3*i], &f[3*i + 1], &f[3*i + 2], &f[3*j], &f[3*j + 1], &f[3*j + 2], &f[3*k], &f[3*k + 1], &f[3*k + 2], &f[3*l], &f[3*l + 1], &f[3*l + 2]);
	}
}

static void
energy(struct params *C, int n, real *r, real *e)
{
	int i;
	int j;
	int k;
	int l;

	for (i = 0; i  < n - 1; i++) {
		j = i + 1;
		e[0] += ebond0(C->kb, C->r0, r[3*i], r[3*i + 1], r[3*i + 2], r[3*j], r[3*j + 1], r[3*j + 2]);
	}

	for (i = 0; i < n - 2; i++) {
		j = i + 1;
		k = i + 2;
		e[1] += eangle0(C->kth, C->th0,  r[3*i], r[3*i + 1], r[3*i + 2],  r[3*j], r[3*j + 1], r[3*j + 2],  r[3*k], r[3*k + 1], r[3*k + 2]);
	}

	for (i = 0; i < n - 3; i++) {
		j = i + 1;
		k = i + 2;
		l = i + 3;
		e[2] += edihedral0(C->kd, C->phi0,
		    r[3*i], r[3*i + 1], r[3*i + 2], 
		    r[3*j], r[3*j + 1], r[3*j + 2], 
		    r[3*k], r[3*k + 1], r[3*k + 2],
		    r[3*l], r[3*l + 1], r[3*l + 2]);
	}
}

int
main(void)
{
	int i;
	int j;
	int n;
	int N;
	real *r;
	real *rp;
	real *f;
	real t[3];
	real e[3];
	real dt;

	struct params C =  {
		.kb = 1.0, .r0 = 1.0, .kth = 0.0, .th0 = 165 *  3.141592 / 180, .kd = -1e-4, .phi0 = 170 * 3.141592 / 180			};
	r = NULL;
	n = N = 0;
	while (scanf("%lf %lf %lf", &t[0], &t[1], &t[2]) == 3) {
		if (n == N) {
			N = 2 * n + 1;
			r = realloc(r, 3 * N * sizeof *r);
			if (r == NULL) {
				fprintf(stderr, "realloc failed\n");
				exit(1);
			}
		}
		r[3 * n] = t[0];
		r[3 * n + 1] = t[1];
		r[3 * n + 2] = t[2];
		n++;
	}
	srand(1234);
	dt = 0.1;
	f = malloc(3 * n * sizeof *f);
	rp = malloc(3 * n * sizeof *f);
	if (f == NULL) {
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}
	for (i = 0; i < 3 * n; i++)
		f[i] = 0;
	for (j = 0; j < 100000000; j++) {
		for (i = 0; i < 3 * n; i++)
			rp[i] = r[i] + dt * f[i] / 2;
		for (i = 0; i < 3 * n; i++)
			f[i] = 0;
		force(&C,n, rp, f);
		for (i = 0; i < 3 * n; i++)
			r[i] += dt * f[i];
		if (j % 100000 == 0) {
			for (i = 0; i < n; i++)
				printf("%.16e %.16e %.16e\n", r[3*i], r[3*i + 1], r[3*i + 2]);
			printf("\n");
			e[0] = e[1] = e[2] = 0;
			energy(&C, n, r, e);
			fprintf(stderr, "%.16e %.16e %.16e %.16e\n", e[0], e[1], e[2], e[0] + e[1] + e[2]);
		}
	}
	free(r);
	free(rp);
	free(f);
}
