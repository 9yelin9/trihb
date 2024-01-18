#ifndef TRIHB_H
#define TRIHB_H

#define USE_MATH_DEFINES

#define BUF_SIZE 1024 // buffer size

#define DIM_L 2 // dimension of lattice
#define DIM_A 2 // dimension of angle
#define DIM_S 3 // dimension of spin

//#define N_THETA 512 // number of theta (0 - pi)
//#define N_PHI 1024 // number of phi (0 - 2*pi)

#define N_THETA 128
#define N_PHI 256

#define N_NN 6 // number of nearest-neighbors
#define RAD_NN 1 // radious of nearest-neighbors

#define ITV_M 5 // itv between monte carlo steps (to avoid autocorrelation)

#define LEN_S 3 // length of s_sub
#define LEN_P 2 // length of psi_coef
#define LEN_E 3 // length of e

#define N(L) (L * L) // size of lattice
#define RHO(N, T, rho1, rho2) ((-2 / (3 * sqrt(3) * N)) * (rho1 + rho2 / T)) // spin stiffness

#include <omp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <dirent.h>
#include <lapack.h>
#include <libgen.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

typedef struct Environment {
	int L; // length of lattice
	int N; // size of lattice
	int M; // num of Monte Carlo iteration
	double h; // magnetic field in z-direction
	double D; // strength of single-ion anisotropy
	double T; // temperature
} Env;

typedef struct Lattice {
	int alpha; // idx of sublattice (0:C, 1:B, 2:A)
	int angle[DIM_A]; // polar and azimuthal angle (theta, phi)
	int nn[N_NN]; // idx of nearest-neighbors
	double site[DIM_L]; // site
	double r_ij[N_NN][DIM_L]; // displacement vector of nearest-neighbors (with PBC)
	//struct Lat *nn[N_NN]; // addresses of nearest-neighbors
} Lat;

typedef struct Observables {
	int M; // num of Monte Carlo iteration
	double mz; // magnetization in z-direction
	double rho1; // dot term of spin stiffness
	double rho2; // cross term of spin stiffness
	double ozz; // C3 lattice symmetry order parameter
} Obs;

// values declared in "lib/libtrihb.c"
extern const double r_lat[DIM_L][DIM_L]; // lattice vectors
extern const double s_sup[DIM_L][DIM_L]; // triangular unit cell vectors
extern const double s_sub[LEN_S][DIM_L]; // triangular unit cell sublattice vectors
extern const double psi_coef[LEN_P][LEN_S]; // coefficients of psi
extern const double e[LEN_E][DIM_L]; // nearest-neighbor unit vectors
extern long seed; // seed for ran2

void ReadLatObs(Env env, char *fn_data, Lat *lat, Obs *obs); // read lattice and observables (in binary)
void WriteLatObs(Env env, char *fn_data, Lat *lat, Obs *obs); // write lattice and observables (in binary)
void InitLat(Env env, Lat *lat); // initialize lattice 
void CalcSpin(int *angle, double *spin); // calculate spin
double CalcEnergyL(Env env, Lat *lat); // calculate hamiltonian overall lattice
double CalcEnergyS(Env env, Lat *lat, int idx); // calculate hamiltonian on site
void MoveAngle(int *angle, int *angle_old); // move angle by ran2
void RestoreAngle(int *angle, int *angle_old); // restore angle (angle = angle_old)
void CalcMz(Env env, Lat *lat, double *mz); // calculate magnetization in z-direction
void CalcRho(Env env, Lat *lat, double *rho1, double *rho2); // calculate expectation values in spin stiffness
void CalcOzz(Env env, Lat *lat, double *ozz); // calculate C3 lattice symmetry order parameter
void MonteCarlo(Env env, Lat *lat); // Monte Carlo algorithm
void RunMonteCarlo(Env env, Lat *lat, Obs *obs); // run Monte Carlo simulation

#endif
