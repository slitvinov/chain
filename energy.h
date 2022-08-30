typedef double real;
struct params {
	real kb;
	real kd;
	real kth;
	real r0;
	real th0;
};
real ebond(struct params*, real, real, real, real, real, real);
real eangle(struct params*, real, real, real, real, real, real, real, real, real);
real edihedral(struct params*, real, real, real,
real, real, real,
real, real, real,
real, real, real);