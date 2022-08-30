#include <stdio.h>
#include <stdlib.h>
#include "energy.h"

static real
uniform(void)
{
	return 2 * (real)rand()/(real)RAND_MAX - 1;
}

static real
eng(struct params *C, int n, real *r)
{
	int i;
	int j;
	real E;

	E = 0;
	for (i = 0; i + 1 < n; i++) {
		j = i + 1;
		E += ebond(C, r[3*i], r[3*i + 1], r[3*i + 2], r[3*j], r[3*j + 1], r[3*j + 2]);
	}
	return E;
}

int
main(void)
{
	int i;
	int j;
	int n;
	int N;
	real *r;
	real *s;
	real *tmp;
	real t[3];
	real alpha;
	real E;
	real Es;

	struct params C =  {
		.kb = 1.0, .r0 = 2.0, .kth = 3.0, .th0 = 4.0, .kd = 5.0 			};
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
	s = malloc(3 * n * sizeof *s);
	alpha = 0.1;
	E = eng(&C, n, r);

	for (j = 0; j < 10000; j++) {
		for (i = 0; i < 3 * n; i++)
			s[i] = r[i] + alpha * uniform();
		Es = eng(&C, n, s);
		if (Es < E) {
			tmp = s;
			s = r;
			r = tmp;
			E = Es;
		}
		if (j % 100 == 0) {
			for (i = 0; i < n; i++)
				printf("%.16e %.16e %.16e\n", r[3*i], r[3*i + 1], r[3*i + 2]);
			printf("\n");
		}
			
		fprintf(stderr, "%g\n", E);
	}
	free(r);
	free(s);
}

