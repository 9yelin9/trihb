#include "trihb.h"
#include "ran2.c"

const double r_lat[DIM_L][DIM_L] = {{1, 0}, {0.5, 0.86602540378443864676}};
const double s_sup[DIM_L][DIM_L] = {{1.5, 0.86602540378443864676}, {1.5, -0.86602540378443864676}};
const double s_sub[LEN_S][DIM_L] = {{0, 0}, {-0.5, 0.86602540378443864676}, {0.5, 0.86602540378443864676}};
const double psi_coef[LEN_P][LEN_S] = {{-4.89897948556635619638, 2.44948974278317809819, 2.44948974278317809819}, {0, 4.24264068711928514640, -4.24264068711928514640}};
const double e[LEN_E][DIM_L] = {{1, 0}, {0.5, 0.86602540378443864676}, {-0.5, 0.86602540378443864676}};
long seed = -1;

void ReadLatObs(Env env, char *fn_data, Lat *lat, Obs *obs) {
	FILE *f = fopen(fn_data, "rb");

	fread(lat, sizeof(Lat), env.N, f);
	fread(obs, sizeof(Obs), 1, f);

	fclose(f);
}

void WriteLatObs(Env env, char *fn_data, Lat *lat, Obs *obs) {
	FILE *f = fopen(fn_data, "wb");

	fwrite(lat, sizeof(Lat), env.N, f);
	fwrite(obs, sizeof(Obs), 1, f);

	fclose(f);
}

void InitLat(Env env, Lat *lat) {
	int i, j, k, a, b, n_nn;
	double norm, r_ij[DIM_L], site[DIM_L];

	// set angle and site
	for(i=0; i<env.N; i++) {
		lat[i].angle[0] = (int)(ran2(&seed) * N_THETA);
		lat[i].angle[1] = (int)(ran2(&seed) * N_PHI);

		for(j=0; j<DIM_L; j++) lat[i].site[j] = (i/env.L) * r_lat[0][j] + (i%env.L) * r_lat[1][j];
	}

	// find nearest-neighbors
	for(i=0; i<env.N; i++) {
		n_nn = 0;
		for(a=0; a<env.N*2; a++) {
			for(k=0; k<DIM_L; k++) site[k] = lat[i].site[k] + (a/env.L-env.L/2) * (env.L*r_lat[0][k]) + (a%env.L-env.L/2) * (env.L*r_lat[1][k]);
			for(j=0; j<env.N; j++) {
				if(i != j) {
					norm = 0;
					for(k=0; k<DIM_L; k++) {
						r_ij[k] = site[k] - lat[j].site[k];
						norm += pow(r_ij[k], 2);
					}
					norm = sqrt(norm);

					if(norm < RAD_NN + 1e-6) {
						lat[i].nn[n_nn] = j;
						for(k=0; k<DIM_L; k++) lat[i].r_ij[n_nn][k] = r_ij[k];
						n_nn++;
					}
				}
			}
		}
		if(n_nn != N_NN) {
			printf("site %d has %d nn\n", i, n_nn);
			exit(1);
		}
	}

	// find alpha
	for(i=0; i<env.N; i++) {
		lat[i].alpha = -1;
		for(a=0; a<LEN_S; a++) {
			for(b=0; b<env.N*2; b++) {
				norm = 0;
				for(k=0; k<DIM_L; k++) norm += fabs(lat[i].site[k] - s_sub[a][k] - (b/env.L-env.L/2) * s_sup[0][k] - (b%env.L-env.L/2) * s_sup[1][k]);
				if(norm < 1e-6) {
					lat[i].alpha = a;
					break;
				}
			}
			if(lat[i].alpha > -1) break;
		}
	}
}

void CalcSpin(int *angle, double *spin) {
	double theta, phi;
	theta = M_PI * angle[0] / N_THETA;
	phi = 2*M_PI * angle[1] / N_PHI;

	spin[0] = sin(theta)*cos(phi);
	spin[1] = sin(theta)*sin(phi);
	spin[2] = cos(theta);
}

double CalcEnergyL(Env env, Lat *lat) {
	int i, j, k;
	double spin_i[DIM_S], spin_j[DIM_S];
	double sum_J=0, sum_h=0, sum_D=0;

	for(i=0; i<env.N; i++) {
		CalcSpin(lat[i].angle, spin_i);
		for(j=0; j<N_NN; j++) {
			CalcSpin(lat[lat[i].nn[j]].angle, spin_j);
			for(k=0; k<DIM_S; k++) sum_J += spin_i[k] * spin_j[k];
		}
		sum_h += spin_i[2];
		sum_D += pow(spin_i[2], 2);
	}

	return sum_J / 2 - env.h * sum_h - env.D * sum_D;
}

double CalcEnergyS(Env env, Lat *lat, int idx) {
	int j, k;
	double spin_i[DIM_S], spin_j[DIM_S];
	double sum_J=0, sum_h=0, sum_D=0;

	CalcSpin(lat[idx].angle, spin_i);
	for(j=0; j<N_NN; j++) {
		CalcSpin(lat[lat[idx].nn[j]].angle, spin_j);
		for(k=0; k<DIM_S; k++) sum_J += spin_i[k] * spin_j[k];
	}
	sum_h += spin_i[2];
	sum_D += pow(spin_i[2], 2);

	return sum_J - env.h * sum_h - env.D * sum_D;
}

void MoveAngle(int *angle, int *angle_old) {
	int i;
	for(i=0; i<DIM_A; i++) angle_old[i] = angle[i];

	angle[0] = (int)(ran2(&seed) * N_THETA);
	angle[1] = (int)(ran2(&seed) * N_PHI);
}

void RestoreAngle(int *angle, int *angle_old) {
	int i;
	for(i=0; i<DIM_A; i++) angle[i] = angle_old[i];
}

void CalcMz(Env env, Lat *lat, double *mz) {
	int i;
	double spin_i[DIM_S];

	for(i=0; i<env.N; i++) {
		CalcSpin(lat[i].angle, spin_i);
		*mz += spin_i[2];
	}
	*mz /= env.N;
}

void CalcRho(Env env, Lat *lat, double *rho1, double *rho2) {
	int i, j, k, a;
	double spin_i[DIM_S], spin_j[DIM_S];
	double coef, dotNN, crsNN, rho1_a, rho2_a;

	for(a=0; a<LEN_E; a++) {
		rho1_a = rho2_a = 0;
		for(i=0; i<env.N; i++) {
			CalcSpin(lat[i].angle, spin_i);
			for(j=0; j<N_NN; j++) {
				coef = dotNN = 0;
				CalcSpin(lat[lat[i].nn[j]].angle, spin_j);
				for(k=0; k<DIM_L; k++) {
					coef += e[a][k] * lat[i].r_ij[j][k];
					dotNN += spin_i[k] * spin_j[k];
				}
				crsNN = spin_i[0] * spin_j[1] - spin_i[1] * spin_j[0];

				rho1_a += pow(coef, 2) * dotNN;
				rho2_a += coef * crsNN;
			}
		}
		*rho1 += rho1_a / 2;
		*rho2 += pow(rho2_a / 2, 2);
	}
}

void CalcOzz(Env env, Lat *lat, double *ozz) {
	int i, j;
	double spin_i[DIM_S], psi[LEN_P]={0};

	for(i=0; i<env.N; i++) {
		CalcSpin(lat[i].angle, spin_i);
		for(j=0; j<LEN_P; j++) {
			psi[j] += psi_coef[j][lat[i].alpha] * spin_i[2];
		}
	}
	for(i=0; i<LEN_P; i++) *ozz += pow(psi[i] / env.N, 2);
	*ozz = sqrt(*ozz);
}

void MonteCarlo(Env env, Lat *lat) {
	int i, idx, angle_old[DIM_A];
	double e0, e1;

	for(i=0; i<env.N; i++) {
		idx = (int)(env.N * ran2(&seed)) % env.N;

		e0 = CalcEnergyS(env, lat, idx);
		MoveAngle(lat[idx].angle, angle_old);
		e1 = CalcEnergyS(env, lat, idx);

		if(e1 - e0 > 0 && ran2(&seed) > exp((e0 - e1) / env.T)) RestoreAngle(lat[idx].angle, angle_old);
	}
}

void RunMonteCarlo(Env env, Lat *lat, Obs *obs) {
	int i;
	Obs obs0;

	for(i=0; i<env.M * ITV_M; i++) {
		MonteCarlo(env, lat);

		memset(&obs0, 0, sizeof(obs0));
		CalcMz(env, lat, &obs0.mz);
		CalcRho(env, lat, &obs0.rho1, &obs0.rho2);
		CalcOzz(env, lat, &obs0.ozz);

		if(i % ITV_M == 0) {
			obs->mz   += obs0.mz;
			obs->rho1 += obs0.rho1;
			obs->rho2 += obs0.rho2;
			obs->ozz  += obs0.ozz;
		}
	}
}
