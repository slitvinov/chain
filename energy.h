typedef double real;
struct params {
	real kb;
	real kd;
	real kth;
	real r0;
	real th0;
};
real ebond0(real, real, real, real, real, real, real, real);
real ebond(struct params*, real, real, real, real, real, real);
real eangle(struct params*, real, real, real, real, real, real, real, real, real);
real eangle0(real, real, real, real, real, real, real, real, real, real, real);
real edihedral0(real, real, real, real,
real, real, real,
real, real, real,
real, real, real);