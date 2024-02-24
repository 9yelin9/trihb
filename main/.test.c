#include "trihb.h"

int main(int argc, char *argv[]) {
	if(argc < 1) {
		printf("%s\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1);

	/*
	Env env = {
		.mcs = 1e3,
		.L = 10,
		.N = env.L * env.L,
		.h = 1.2,
		.D = 0,
	};

	Lat lat[N_REP][env.N];
	Obs obs[N_REP];

	int m;
	double T_min=0.0001, T_max=2, T[N_REP];
	for(m=0; m<N_REP; m++) T[m] = T_min + (T_max - T_min) * m / N_REP;

	InitLat(&env, lat); InitT(&env, lat, T);

	double Mz, Rho, Ozz;
	RunMonteCarlo(&env, lat, obs, T, env.mcs, env.mcs+1);

	for(m=0; m<N_REP; m++) {
		Mz  = obs[m].mz / env.mcs;	
		Rho = (-2/(3*sqrt(3)*env.N)) * (obs[m].rho1/env.mcs + (obs[m].rho2/env.mcs)/T[m]);
		Ozz = obs[m].ozz / env.mcs;

		printf("%16f%16f%16f%16f%16f%16f\n", env.h, env.D, T[m], Mz, Rho, Ozz);
	}

	int i, j;
	double spin_i[DIM_S];
	for(i=0; i<env.N; i++) {
		//lat[i].angle[0] = 0;
		//lat[i].angle[0] = N_THETA / 2;
		lat[i].angle[0] = i % 2 ? N_THETA : 0;
		lat[i].angle[1] = 0;

		CalcSpin(lat[i].angle, spin_i);
		for(j=0; j<DIM_L; j++) printf("%f\t", lat[i].site[j]);
		for(j=0; j<DIM_S; j++) printf("%f\t", spin_i[j]);
		printf("\n");
	}

	CalcMz(env, lat, &obs.mz);
	CalcRho(env, lat, &obs.rho1, &obs.rho2);
	CalcOzz(env, lat, &obs.ozz);

	printf("L=%d : mz=%f rho1=%f rho2=%f rho=%f ozz=%f\n", env.L, obs.mz, obs.rho1, obs.rho2, RHO(env.N, env.T, obs.rho1, obs.rho2), obs.ozz);
	*/

	return 0;
}
