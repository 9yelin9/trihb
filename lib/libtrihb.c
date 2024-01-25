#include "trihb.h"
#include "ran2.c"

const double r[DIM_L][DIM_L] = {{1, 0}, {0.5, root3h}}; // lattice vector
const double r_sup[DIM_L][DIM_L] = {{1.5, root3h}, {1.5, -root3h}}; // C3 supercell lattice vector
const double r_sub[N_SUB][DIM_L] = {{0, 0}, {-0.5, root3h}, {0.5, root3h}}; // C3 sublattice vector (C, A, B)

const double psi_coef[2][N_SUB] = {{-2*root6, root6, root6}, {0, -3*root2, 3*root2}}; // coefficient of C3 order parameter (C, A, B)
const double e[3][DIM_L] = {{1, 0}, {0.5, root3h}, {-0.5, root3h}}; // unit vector

long seed = -1;

void ReadLatObs(Env env, Lat *lat, Obs *obs, char *fn) {
	hid_t file_id, dataset_id;

	file_id = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT); 

	dataset_id = H5Dopen2(file_id, "/lat", H5P_DEFAULT);
	H5Dread(dataset_id, H5Dget_type(dataset_id), H5S_ALL, H5S_ALL, H5P_DEFAULT, lat);

	dataset_id = H5Dopen2(file_id, "/obs", H5P_DEFAULT);
	H5Dread(dataset_id, H5Dget_type(dataset_id), H5S_ALL, H5S_ALL, H5P_DEFAULT, obs);

	H5Dclose(dataset_id);
	H5Fclose(file_id);
}

void SaveLatObs(Env env, Lat *lat, Obs *obs, char *fn) {
	hid_t file_id, datatype_id, dataset_id, dataspace_id;

	file_id = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// Lat
	datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(Lat));
	H5Tinsert(datatype_id, "alpha", HOFFSET(Lat, alpha), H5T_NATIVE_INT);
	H5Tinsert(datatype_id, "angle", HOFFSET(Lat, angle), H5Tarray_create2(H5T_NATIVE_INT,    1, (hsize_t[1]){DIM_A}));
	H5Tinsert(datatype_id, "nn",    HOFFSET(Lat, nn),    H5Tarray_create2(H5T_NATIVE_INT,    1, (hsize_t[1]){N_NN}));
	H5Tinsert(datatype_id, "site",  HOFFSET(Lat, site),  H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, (hsize_t[1]){DIM_L}));
	H5Tinsert(datatype_id, "r_ij",  HOFFSET(Lat, r_ij),  H5Tarray_create2(H5T_NATIVE_DOUBLE, 2, (hsize_t[2]){N_NN, DIM_L}));

	dataspace_id = H5Screate_simple(1, (hsize_t[1]){env.N}, NULL);
	dataset_id = H5Dcreate2(file_id, "/lat", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, lat);

	// Obs
	datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(Obs));
	H5Tinsert(datatype_id, "mz",   HOFFSET(Obs, mz),   H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_id, "rho1", HOFFSET(Obs, rho1), H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_id, "rho2", HOFFSET(Obs, rho2), H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_id, "ozz",  HOFFSET(Obs, ozz),  H5T_NATIVE_DOUBLE);

	dataspace_id = H5Screate_simple(1, (hsize_t[1]){1}, NULL);
	dataset_id = H5Dcreate2(file_id, "/obs", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, lat);

	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Tclose(datatype_id);
	H5Fclose(file_id);
}

void InitLat(Env env, Lat *lat) {
	int i, j, k, a, b, cnt_alpha, cnt_nn;
	double d2, r_lat[DIM_L][DIM_L], r_ij[DIM_L];

	// set lattice displacement vector
	for(i=0; i<DIM_L; i++) for(j=0; j<DIM_L; j++)
		r_lat[i][j] = env.L * r[i][j];

	// set angle, site and alpha
	for(i=0; i<env.N; i++) {
		lat[i].angle[0] = (int)(ran2(&seed) * N_THETA);
		lat[i].angle[1] = (int)(ran2(&seed) * N_PHI);

		for(k=0; k<DIM_L; k++)
			lat[i].site[k] = (i/env.L) * r[0][k] + (i%env.L) * r[1][k];
		
		cnt_alpha = 0;
		for(j=0; j<N_SUB; j++) {
			for(a=-1; a<2; a++) for(b=-1; b<2; b++) {
				d2 = 0;
				for(k=0; k<DIM_L; k++) {
					r_ij[k] = a * r_sup[0][k] + b * r_sup[1][k] + r_sub[j][k] - lat[i].site[k];
					d2 += pow(r_ij[k], 2);
				}
				if(d2 < 1e-6) {
					lat[i].alpha = j;
					cnt_alpha++;
				}
			}
		}
		if(cnt_alpha != 1) {
			printf("Site %d has %d alpha\n", i, cnt_alpha);
			exit(1);
		}
	}

	// find nearest-neighbors
	for(i=0; i<env.N; i++) {
		cnt_nn = 0;
		for(j=0; j<env.N; j++) {
			if(i != j) {
				for(a=-1; a<2; a++) for(b=-1; b<2; b++) {
					d2 = 0;
					for(k=0; k<DIM_L; k++) {
						r_ij[k] = a * r_lat[0][k] + b * r_lat[1][k] + lat[j].site[k] - lat[i].site[k];
						d2 += pow(r_ij[k], 2);
					}
					if(d2 < D2_NN + 1e-6) {
						lat[i].nn[cnt_nn] = j;
						for(k=0; k<DIM_L; k++) lat[i].r_ij[cnt_nn][k] = r_ij[k];
						cnt_nn++;
					}
				}
			}
			if(cnt_nn != N_NN) {
				printf("Site %d has %d nn\n", i, cnt_nn);
				exit(1);
			}
		}
	}
}

void GetSpin(int *angle, double *spin) {
	double theta, phi;
	theta = M_PI * angle[0] / N_THETA;
	phi = 2*M_PI * angle[1] / N_PHI;

	spin[0] = sin(theta) * cos(phi);
	spin[1] = sin(theta) * sin(phi);
	spin[2] = cos(theta);
}

double CalcEnergy(Env env, Lat *lat, int idx) {
	int j, k;
	double spin_i[DIM_S], spin_j[DIM_S];
	double sum_J=0, sum_h=0, sum_D=0;

	GetSpin(lat[idx].angle, spin_i);
	for(j=0; j<N_NN; j++) {
		GetSpin(lat[lat[idx].nn[j]].angle, spin_j);
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

void RestAngle(int *angle, int *angle_old) {
	int i;
	for(i=0; i<DIM_A; i++) angle[i] = angle_old[i];
}

void CalcMz(Env env, Lat *lat, double *mz) {
	int i;
	double spin_i[DIM_S];

	for(i=0; i<env.N; i++) {
		GetSpin(lat[i].angle, spin_i);
		*mz += spin_i[2];
	}
	*mz /= env.N;
}

void CalcRho(Env env, Lat *lat, double *rho1, double *rho2) {
	int i, j, k, a;
	double spin_i[DIM_S], spin_j[DIM_S];
	double coef, dotNN, crsNN, rho1_a, rho2_a;

	for(a=0; a<3; a++) {
		rho1_a = rho2_a = 0;
		for(i=0; i<env.N; i++) {
			GetSpin(lat[i].angle, spin_i);
			for(j=0; j<N_NN; j++) {
				GetSpin(lat[lat[i].nn[j]].angle, spin_j);

				coef = dotNN = 0;
				for(k=0; k<DIM_L; k++) {
					coef += e[a][k] * lat[i].r_ij[j][k];
					dotNN += spin_i[k] * spin_j[k];
				}
				crsNN = spin_i[0] * spin_j[1] - spin_i[1] * spin_j[0];

				rho1_a += pow(coef, 2) * dotNN;
				rho2_a += coef * crsNN;
			}
		}
		*rho1 += rho1_a/2;
		*rho2 += pow(rho2_a/2, 2);
	}
}

void CalcOzz(Env env, Lat *lat, double *ozz) {
	int i, j;
	double spin_i[DIM_S], psi[2]={0};

	for(i=0; i<env.N; i++) {
		GetSpin(lat[i].angle, spin_i);
		for(j=0; j<2; j++) psi[j] += psi_coef[j][lat[i].alpha] * spin_i[2];
	}
	for(i=0; i<2; i++) *ozz += pow(psi[i]/env.N, 2);
	*ozz = sqrt(*ozz);
}

void MonteCarlo(Env env, Lat *lat) {
	int i, idx, angle_old[DIM_A];
	double e0, e1;

	for(i=0; i<env.N; i++) {
		idx = (int)(ran2(&seed) * env.N);

		e0 = CalcEnergy(env, lat, idx);
		MoveAngle(lat[idx].angle, angle_old);
		e1 = CalcEnergy(env, lat, idx);

		if(e1 - e0 > 0 && ran2(&seed) > exp((e0 - e1) / env.T)) RestAngle(lat[idx].angle, angle_old);
	}
}

void RunMonteCarlo(Env env, Lat *lat, Obs *obs) {
	int i;
	Obs obs0;

	for(i=0; i<env.Mmc * ITV; i++) {
		MonteCarlo(env, lat);

		memset(&obs0, 0, sizeof(obs0));
		CalcMz(env, lat, &obs0.mz);
		CalcRho(env, lat, &obs0.rho1, &obs0.rho2);
		CalcOzz(env, lat, &obs0.ozz);

		if(i % ITV == 0) {
			obs->mz   += obs0.mz;
			obs->rho1 += obs0.rho1;
			obs->rho2 += obs0.rho2;
			obs->ozz  += obs0.ozz;
		}
	}
}
