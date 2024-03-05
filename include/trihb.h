#ifndef TRIHB_H
#define TRIHB_H

#define USE_MATH_DEFINES

#define root2  1.41421356237309504880
#define root3  1.73205080756887729352 
#define root3h 0.86602540378443864676
#define root6  2.44948974278317809819

#define BUF_SIZE 1024 // buffer size

#define DIM_L 2 // dimension of lattice
#define DIM_A 2 // dimension of angle
#define DIM_S 3 // dimension of spin

#define N_SUB 3 // number of sites in sublattice

#define N_NN 6 // number of nearest-neighbors
#define D2_NN 1 // squared displacement of nearest-neighbors

#define DT 0.1 // interval between temperatures

#define ECS 10 // interval between exchange MC steps

#define RAND (2*ran2(&seed)-1) // random number between -1*RAND_LIMIT and 1*RAND_LIMIT
#define FLIP_LIMIT 0.3 // flip limitation

#define DELTA(lat, T, m) ((1/T[m] - 1/T[m+1]) * (CalcEnergy(env, lat[m]) - CalcEnergy(env, lat[m+1]))) // exchange probability factor
#define RHO(N, T, rho1, rho2) ((-2/(3*sqrt(3)*N)) * (rho1 + rho2/T)) // spin stiffness rho

#define MIN(x, y) (x < y ? x : y)

#include <omp.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <dirent.h>
#include <libgen.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

typedef struct Environment {
	char dir_save[BUF_SIZE]; // directory to save
	int L; // length of lattice
	int N; // size of lattice
	int M; // num of replicas
	int m; // index of replica at T=env.T
	int eqs; // num of equilibration steps
	int mcs; // num of MC steps
	double h; // magnetic field in z-direction
	double D; // strength of single-ion anisotropy
	double T; // temperature
} Env;

typedef struct Lattice {
	int alpha; // site in sublattice (0:C, 1:B, 2:A)
	int nn[N_NN]; // idx of nearest-neighbors
	double site[DIM_L]; // site
	double angle[DIM_A]; // polar and azimuthal angle (theta, phi)
	double r_ij[N_NN][DIM_L]; // displacement vector of nearest-neighbors in PBC
} Lat;

typedef struct Observable {
	double mz; // magnetization in z-direction
	double rho1; // dot term of spin stiffness
	double rho2; // cross term of spin stiffness
	double ozz; // C3 lattice symmetry order parameter
} Obs;

// values declared in "lib/libtrihb.c"
extern const double r_lat[DIM_L][DIM_L]; // lattice vectors
extern const double s_sup[DIM_L][DIM_L]; // triangular unit cell vectors
extern const double s_sub[N_SUB][DIM_L]; // triangular unit cell sublattice vectors
extern const double psi_coef[2][N_SUB]; // coefficients of psi
extern const double e[3][DIM_L]; // nearest-neighbor unit vectors
extern long seed; // seed for ran2

void SaveLatObs(Env *env, Lat **lat, Obs *obs, double *T); // write lattice and observables (in hdf5)
void CalcSpin(double *angle, double *spin); // calcusitee spin
double CalcEnergy(Env *env, Lat *lat); // calcusitee hamiltonian
double CalcEnergySite(Env *env, Lat *lat, int idx); // calcusitee hamiltonian on site
void FlipAngle(double *angle, double *angle_old); // flip angle
void UndoAngle(double *angle, double *angle_old); // undo angle (angle = angle_old)
void CalcMz(Env *env, Lat *lat, double *mz); // calcusitee magnetization in z-direction
void CalcRho(Env *env, Lat *lat, double *rho1, double *rho2); // calcusitee expectation values in spin stiffness
void CalcOzz(Env *env, Lat *lat, double *ozz); // calcusitee C3 sitetice symmetry order parameter
void InitLat(Env *env, Lat **lat); // initialize lattice
void InitT(Env *env, Lat **lat, double *T); // initialize T
void MonteCarlo(Env *env, Lat *lat, double T); // MC algorithm
void Exchange(Env *env, Lat **lat, double *T); // exchange MC algorithm
void RunEquilibration(Env *env, Lat **lat, double *T); // run equilibration
void RunMonteCarlo(Env *env, Lat **lat, Obs *obs, double *T); // run MC simulation

#endif
