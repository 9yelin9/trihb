#include "trihb.h"
#include "ran2.c"

const double r[DIM_L][DIM_L] = {{1, 0}, {0.5, root3h}}; // lattice vector
const double r_sup[DIM_L][DIM_L] = {{1.5, root3h}, {1.5, -root3h}}; // C3 supercell lattice vector
const double r_sub[N_SUB][DIM_L] = {{0, 0}, {0.5, root3h}, {-0.5, root3h}}; // C3 sublattice vector (C, A, B)
const double psi_coef[2][N_SUB] = {{-2*root6, root6, root6}, {0, 3*root2, -3*root2}}; // coefficient of C3 order parameter (C, A, B)
const double e[3][DIM_L] = {{1, 0}, {0.5, root3h}, {-0.5, root3h}}; // unit vector
long seed = -1;
#pragma omp threadprivate(seed)

void TestLat(Env *env, Lat **lat, int m, int site) {
	int i;

	printf("=============== rep %d : site %d ===============\n", m, site);
	printf("alpha = %d\n", lat[m][site].alpha);
	printf("site = ");
	for(i=0; i<DIM_L; i++) printf("%f\t", lat[m][site].site[i]);
	printf("\nangle = ");
	for(i=0; i<DIM_L; i++) printf("%f\t", lat[m][site].angle[i]);
	printf("\nnn = ");
	for(i=0; i<N_NN; i++) printf("%d\t", lat[m][site].nn[i]);
	printf("\n\n");
}

void TestAlpha(Env *env, Lat **lat, int m) {
	int i, j;

	printf("=============== rep %d  ===============\n", m);
	for(i=0; i<env->L; i++) {
		for(j=0; j<i; j++) printf("  ");
		for(j=0; j<env->L; j++) {
			printf("%d   ", lat[m][i*env->L+j].alpha);
		}
		printf("\n");
	}
	printf("\n");
}

void GetFileName(Env *env, char *dn, char *ext, char *fn) {
	char dir_save[BUF_SIZE];
	sprintf(dir_save, "%s/%s", env->dir_save, dn);
	if(-access(dir_save, 0)) mkdir(dir_save, 0755);
	sprintf(fn, "%s/h%.4f_D%.4f_T%.4f.%s", dir_save, env->h, env->D, env->T, ext);
}

void SaveLatObs(Env *env, Lat **lat, Obs *obs, double *T) {
	char fn[BUF_SIZE];
	hid_t file_id, datatype_lat_id, datatype_obs_id, dataset_id, dataspace_id;

	GetFileName(env, "res", "h5", fn);
	file_id = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// lat
	datatype_lat_id = H5Tcreate(H5T_COMPOUND, sizeof(Lat));
	H5Tinsert(datatype_lat_id, "alpha", HOFFSET(Lat, alpha), H5T_NATIVE_INT);
	H5Tinsert(datatype_lat_id, "nn",    HOFFSET(Lat, nn),    H5Tarray_create2(H5T_NATIVE_INT,    1, (hsize_t[1]){N_NN}));
	H5Tinsert(datatype_lat_id, "site",  HOFFSET(Lat, site),  H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, (hsize_t[1]){DIM_L}));
	H5Tinsert(datatype_lat_id, "angle", HOFFSET(Lat, angle), H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, (hsize_t[1]){DIM_A}));
	H5Tinsert(datatype_lat_id, "r_ij",  HOFFSET(Lat, r_ij),  H5Tarray_create2(H5T_NATIVE_DOUBLE, 2, (hsize_t[2]){N_NN, DIM_L}));

	// obs
	datatype_obs_id = H5Tcreate(H5T_COMPOUND, sizeof(Obs));
	H5Tinsert(datatype_obs_id, "mz",   HOFFSET(Obs, mz),   H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_obs_id, "rho1", HOFFSET(Obs, rho1), H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_obs_id, "rho2", HOFFSET(Obs, rho2), H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_obs_id, "ozz",  HOFFSET(Obs, ozz),  H5T_NATIVE_DOUBLE);

	dataspace_id = H5Screate_simple(1, (hsize_t[1]){env->N}, NULL);
	dataset_id = H5Dcreate2(file_id, "/lat", datatype_lat_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_id, datatype_lat_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, lat[env->m]);

	dataspace_id = H5Screate_simple(1, (hsize_t[1]){1}, NULL);
	dataset_id = H5Dcreate2(file_id, "/obs", datatype_obs_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_id, datatype_obs_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &obs[env->m]);

	H5Tclose(datatype_lat_id);
	H5Tclose(datatype_obs_id);
	H5Sclose(dataspace_id);
	H5Dclose(dataset_id);
	H5Fclose(file_id);

	printf("Save as \"%s\"\n", fn);
}

void GetSpin(double *angle, double *spin) {
	spin[0] = sin(angle[0]) * cos(angle[1]);
	spin[1] = sin(angle[0]) * sin(angle[1]);
	spin[2] = cos(angle[0]);
}

double CalcEnergy(Env *env, Lat *lat) {
	int i, j, k;
	double spin_i[DIM_S], spin_j[DIM_S];
	double sum_J=0, sum_h=0, sum_D=0;

	for(i=0; i<env->N; i++) {
		GetSpin(lat[i].angle, spin_i);
		for(j=0; j<N_NN; j++) {
			GetSpin(lat[lat[i].nn[j]].angle, spin_j);
			for(k=0; k<DIM_S; k++) sum_J += spin_i[k] * spin_j[k];
		}
		sum_h += spin_i[2];
		sum_D += pow(spin_i[2], 2);
	}

	return sum_J/2 - env->h * sum_h - env->D * sum_D;
}

double CalcEnergySite(Env *env, Lat *lat, int idx) {
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

	return sum_J - env->h * sum_h - env->D * sum_D;
}

void FlipAngle(double *angle, double *angle_old) {
	int i;
	double dot, spin[DIM_S], spin_old[DIM_S];
	for(i=0; i<DIM_A; i++) angle_old[i] = angle[i];
	GetSpin(angle_old, spin_old);

	while(1) {
		angle[0] = acos(RAND);
		angle[1] = M_PI * RAND;
		GetSpin(angle, spin);

		dot = 0;
		for(i=0; i<DIM_S; i++) dot += spin[i] * spin_old[i];
		if(dot - (1 - 2*FLIP_LIMIT) > 0) break;
	}
}

void UndoAngle(double *angle, double *angle_old) {
	int i;
	for(i=0; i<DIM_A; i++) angle[i] = angle_old[i];
}

void CalcMz(Env *env, Lat *lat, double *mz) {
	int i;
	double spin_i[DIM_S];

	for(i=0; i<env->N; i++) {
		GetSpin(lat[i].angle, spin_i);
		*mz += spin_i[2];
	}
	*mz /= env->N;
}

void CalcRho(Env *env, Lat *lat, double *rho1, double *rho2) {
	int i, j, k, a;
	double spin_i[DIM_S], spin_j[DIM_S];
	double coef, dotNN, crsNN, rho1_a, rho2_a;

	for(a=0; a<3; a++) {
		rho1_a = rho2_a = 0;
		for(i=0; i<env->N; i++) {
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

void CalcOzz(Env *env, Lat *lat, double *ozz) {
	int i, j;
	double spin_i[DIM_S], psi[2]={0};

	for(i=0; i<env->N; i++) {
		GetSpin(lat[i].angle, spin_i);
		for(j=0; j<2; j++) psi[j] += psi_coef[j][lat[i].alpha] * spin_i[2];
	}
	for(i=0; i<2; i++) *ozz += pow(psi[i]/env->N, 2);
	*ozz = sqrt(*ozz);
}

void InitLat(Env *env, Lat **lat) {
	int i, j, k, a, b, m, cnt_alpha, cnt_nn;
	double d2, r_lat[DIM_L][DIM_L], r_ij[DIM_L];

	// set lattice displacement vector
	for(i=0; i<DIM_L; i++) for(j=0; j<DIM_L; j++)
		r_lat[i][j] = env->L * r[i][j];

	// set angle, site and alpha
	for(i=0; i<env->N; i++) {
		lat[0][i].angle[0] = acos(RAND);
		lat[0][i].angle[1] = M_PI * RAND;

		for(k=0; k<DIM_L; k++)
			lat[0][i].site[k] = (i/env->L) * r[0][k] + (i%env->L) * r[1][k];
		
		cnt_alpha = 0;
		for(j=0; j<N_SUB; j++) {
			for(a=-env->L; a<env->L; a++) for(b=-env->L; b<env->L; b++) {
				d2 = 0;
				for(k=0; k<DIM_L; k++) {
					r_ij[k] = a * r_sup[0][k] + b * r_sup[1][k] + r_sub[j][k] - lat[0][i].site[k];
					d2 += pow(r_ij[k], 2);
				}
				if(d2 < 1e-6) {
					lat[0][i].alpha = j;
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
	for(i=0; i<env->N; i++) {
		cnt_nn = 0;
		for(j=0; j<env->N; j++) {
			if(i != j) {
				for(a=-1; a<2; a++) for(b=-1; b<2; b++) {
					d2 = 0;
					for(k=0; k<DIM_L; k++) {
						r_ij[k] = a * r_lat[0][k] + b * r_lat[1][k] + lat[0][j].site[k] - lat[0][i].site[k];
						d2 += pow(r_ij[k], 2);
					}
					if(d2 < D2_NN + 1e-6) {
						lat[0][i].nn[cnt_nn] = j;
						for(k=0; k<DIM_L; k++) lat[0][i].r_ij[cnt_nn][k] = r_ij[k];
						cnt_nn++;
					}
				}
			}
		}
		if(cnt_nn != N_NN) {
			printf("Site %d has %d nn\n", i, cnt_nn);
			exit(1);
		}
	}

	// initialize other replicas
	for(m=1; m<env->M; m++) {
		for(i=0; i<env->N; i++) {
			lat[m][i].angle[0] = acos(RAND);
			lat[m][i].angle[1] = M_PI * RAND;

			lat[m][i].alpha = lat[0][i].alpha;
			for(j=0; j<DIM_L; j++) lat[m][i].site[j] = lat[0][i].site[j];
			for(j=0; j<N_NN; j++) {
				lat[m][i].nn[j] = lat[0][i].nn[j];
				for(k=0; k<DIM_L; k++) {
					lat[m][i].r_ij[j][k] = lat[0][i].r_ij[j][k];
				}
			}
		}
	}
}

void MonteCarlo(Env *env, Lat *lat, double T) {
	int i, idx;
	double e0, e1, angle_old[DIM_A];

	for(i=0; i<env->N; i++) {
		idx = (int)(ran2(&seed) * env->N);

		e0 = CalcEnergySite(env, lat, idx);
		FlipAngle(lat[idx].angle, angle_old);
		e1 = CalcEnergySite(env, lat, idx);

		if(e1 - e0 > 0 && ran2(&seed) > exp((e0 - e1) / T)) UndoAngle(lat[idx].angle, angle_old);
	}
}

void Exchange(Env *env, Lat **lat, double *T) {
	int m = (int)(ran2(&seed) * (env->M-1));
	Lat *lat_tmp;

	if(ran2(&seed) < exp(DELTA(lat, T, m))) {
		lat_tmp = &(*lat[m]);
		*lat[m] = *lat[m+1];
		*lat[m+1] = *lat_tmp;
	}
}

void RunEquilibration(Env *env, Lat **lat, double *T) {
	int i, m;

	for(i=0; i<env->eqs; i++) {
#pragma omp parallel for ordered num_threads(env->M)
		for(m=0; m<env->M; m++) MonteCarlo(env, lat[m], T[m]);
	}
}

void RunMonteCarlo(Env *env, Lat **lat, Obs *obs, double *T) {
	FILE *f; 
	char fn[BUF_SIZE];
	GetFileName(env, "log", "txt", fn);
	f = fopen(fn, "w");

	int i, m;
	Obs *obs0=(Obs*)malloc(sizeof(Obs) * env->M);

	fprintf(f, "%10s%12s%12s%12s%12s%12s\n", "itr", "e", "mz", "rho1", "rho2", "ozz");
	for(i=1; i<env->mcs+1; i++) {
		memset(obs0, 0, sizeof(Obs) * env->M);

#pragma omp parallel for ordered num_threads(env->M)
		for(m=0; m<env->M; m++) {
			MonteCarlo(env, lat[m], T[m]);

			CalcMz(env, lat[m], &obs0[m].mz);
			CalcRho(env, lat[m], &obs0[m].rho1, &obs0[m].rho2);
			CalcOzz(env, lat[m], &obs0[m].ozz);
		}

		for(m=0; m<env->M; m++) {
			obs[m].mz   += obs0[m].mz;
			obs[m].rho1 += obs0[m].rho1;
			obs[m].rho2 += obs0[m].rho2;
			obs[m].ozz  += obs0[m].ozz;
		}

		fprintf(f, "%10d%12f%12f%12f%12f%12f\n", i, CalcEnergy(env, lat[env->m]), obs0[env->m].mz, obs0[env->m].rho1, obs0[env->m].rho2, obs0[env->m].ozz);
		if(i % ECS == 0) Exchange(env, lat, T);
	}

	for(m=0; m<env->M; m++) {
		obs[m].mz   /= env->mcs;
		obs[m].rho1 /= env->mcs;
		obs[m].rho2 /= env->mcs;
		obs[m].ozz  /= env->mcs;
	}

	free(obs0);
	fclose(f);

	printf("Save as \"%s\"\n", fn);
}

