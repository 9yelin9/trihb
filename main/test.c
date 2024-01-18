#include "trihb.h"

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <L>\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1);

	Env env = {
		.L = atoi(argv[1]),
		.N = N(env.L),
		.M = 1e3,
		.h = 1.2,
		.D = 0,
		.T = -1,
	};
	Lat lat[env.N];
	Obs obs={0};

	InitLat(env, lat);

	int i;
	double Mz, Rho, Ozz;
	for(i=10; i>0; i--) {
		env.T = 0.05 * i;
		memset(&obs, 0, sizeof(obs));
		RunMonteCarlo(env, lat, &obs);

		Mz  = obs.mz / env.M;	
		Rho = RHO(env.N, env.T, obs.rho1 / env.M, obs.rho2 / env.M);
		Ozz = obs.ozz / env.M;

		printf("%16f%16f%16f%16f%16f%16f\n", env.h, env.D, env.T, Mz, Rho, Ozz);
	}

	/*
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

	printf("L=%d : e=%f mz=%f rho1=%f rho2=%f rho=%f ozz=%f\n", env.L, CalcEnergyL(env, lat), obs.mz, obs.rho1, obs.rho2, RHO(env.N, env.T, obs.rho1, obs.rho2), obs.ozz);
	*/

	return 0;
}
