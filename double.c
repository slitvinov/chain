#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "energy.h"

static real
uniform0(void)
{
	return (real)rand()/(real)RAND_MAX;
}

static real
uniform(void)
{
	return 2 * uniform0() - 1;
}

static real
eng(struct params *C, int m, int *connect, real *r)
{
	int q;
	int i;
	int j;
	int k;
	real E;

	E = 0;
	for (q = 0; q  < m; q++) {
		i = connect[2* q];
		j = connect[2 * q + 1];
		E += ebond0(C->kb, C->r0, r[3*i], r[3*i + 1], r[3*i + 2], r[3*j], r[3*j + 1], r[3*j + 2]);
	}

	return E;
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
	real alpha;
	real E;
	real Es;
	real dE;
	real *r;
	real *s;
	real t[3];
	real Temp;
	real *tmp;
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
	file = fopen(path, "r");
	if (file == NULL) {
		fprintf(stderr, "fail to open '%s'\n", path);
		exit(1);
	}
	m = M = 0;
	while (fscanf(file, "%d %d\n", &c0, &c1) == 2) {
		if (m == M) {
			M = 2 * n + 1;
			connect = realloc(connect, 2 * M * sizeof *connect);
		}
		connect[2 * m] = c0;
		connect[2 * m + 1] =c1;
		m++;
	}
	fclose(file);
	fprintf(stderr, "nm: %d %d\n", n, m);
	srand(1234);
	alpha = 0.01;
	Temp = 0.0001;
	s = malloc(3 * n * sizeof *s);
	E = eng(&C, m, connect, r);
	for (j = 0; j < 10000000; j++) {
		for (i = 0; i < 3 * n; i++)
			s[i] = r[i] + alpha * uniform();
		Es = eng(&C, m, connect, s);
		dE = Es - E;
		if (dE < 0 || uniform0() < exp(-dE/Temp)  ) {
			tmp = s;
			s = r;
			r = tmp;
			E = Es;
		}
		if (j % 100000 == 0) {
			for (i = 0; i < n; i++)
				printf("%.16e %.16e %.16e\n", r[3*i], r[3*i + 1], r[3*i + 2]);
			printf("\n");
			fprintf(stderr, "Eng: %g %g\n", E, dE);
		}
	}
	free(r);
	free(s);
	free(connect);
}
