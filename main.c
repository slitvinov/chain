#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "energy.h"
#include "force.h"

static void
force(struct params *C, int m, real *r0, int *connect, real *r, real *f)
{
	int q;
	int i;
	int j;

	for (q = 0; q  < m; q++) {
		i = connect[2* q];
		j = connect[2 * q + 1];
		fbond0(C->kb, r0[q], r[3*i], r[3*i + 1], r[3*i + 2], r[3*j], r[3*j + 1], r[3*j + 2],
							&f[3*i], &f[3*i + 1], &f[3*i + 2], &f[3*j], &f[3*j + 1], &f[3*j + 2]);
	}

}

int
main(void)
{
	int c0;
	int c1;
	int i;
	int j;
	int n;
	int m;
	int N;
	int M;
	real *r;
	real *f;
	real t[3];
	real dt;
	real *r0;
	real r00;
	int *connect;
	FILE *file;
	char *path = "connect";

	struct params C =  {
		.kb = 1.0, .r0 = 1.0					};
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
	connect = NULL;
	r0 = NULL;
	file = fopen(path, "r");
	if (file == NULL) {
		fprintf(stderr, "fail to open '%s'\n", path);
		exit(1);
	}
	m = M = 0;
	while (fscanf(file, "%d %d %lf\n", &c0, &c1, &r00) == 3) {
		if (m == M) {
			M = 2 * n + 1;
			connect = realloc(connect, 2 * M * sizeof *connect);
			r0 = realloc(r0, M * sizeof *r0);
		}
		connect[2 * m] = c0;
		connect[2 * m + 1] =c1;
		r0[m] = r00;
		m++;
	}
	fclose(file);
	fprintf(stderr, "nm: %d %d\n", n, m);
	srand(1234);
	dt = 0.1;
	f = malloc(3 * n * sizeof *f);
	for (j = 0; j < 10000; j++) {
		for (i = 0; i < 3 * n; i++)
			f[i] = 0;
		force(&C, m, r0, connect, r, f);
		for (i = 0; i < 3 * n; i++)
			r[i] += dt * f[i];
		if (j % 100 == 0) {
			for (i = 0; i < n; i++)
				printf("%.16e %.16e %.16e\n", r[3*i], r[3*i + 1], r[3*i + 2]);
			printf("\n");
		}
	}
	free(r);
	free(f);
	free(connect);
}
