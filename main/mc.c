#include "trihb.h"

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <L> <M> <eqs> <mcs> <h> <D> <T>\n", argv[0]);
		exit(1);
	}

	int m;
	Env env = {
		.L = atoi(argv[1]),
		.N = env.L * env.L,
		.M = atoi(argv[2]),
		.m = -1,
		.eqs = (int)atof(argv[3]),
		.mcs = (int)atof(argv[4]),
		.h = atof(argv[5]),
		.D = atof(argv[6]),
		.T = atof(argv[7]),
	};
	sprintf(env.dir_save, "data/L%d_M%d_eqs%.e_mcs%.e", env.L, env.M, (double)env.eqs, (double)env.mcs);
	if(-access(env.dir_save, 0)) mkdir(env.dir_save, 0755);

	double T[env.M], T_min=env.T-DT*(env.M/2);
	for(m=0; m<env.M; m++) {
		T[m] = T_min + DT*m;
		if(fabs(T[m] - env.T) < 1e-6) env.m = m;
	}
	if(env.m < 0) {
		printf("No replica exists at T=%.4f\n", env.T);
		exit(1);
	}

	Lat **lat=(Lat**)malloc(sizeof(Lat*) * env.M);
	for(m=0; m<env.M; m++) lat[m] = (Lat*)malloc(sizeof(Lat) * env.N);
	Obs *obs=(Obs*)malloc(sizeof(Obs) * env.M);

	time_t t0=time(NULL);

	InitLat(&env, lat);
	RunEquilibration(&env, lat, T);
	RunMonteCarlo(&env, lat, obs, T);
	SaveLatObs(&env, lat, obs, T);

	printf("%12s%12s%12s%12s%12s%12s\n", "T", "mz", "rho1", "rho2", "rho", "ozz");
	printf("%12f%12f%12f%12f%12f%12f\n", T[env.m], obs[env.m].mz, obs[env.m].rho1, obs[env.m].rho2, RHO(env.N, T[env.m], obs[env.m].rho1, obs[env.m].rho2), obs[env.m].ozz);
	for(m=0; m<env.M; m++) free(lat[m]);
	free(lat);
	free(obs);

	time_t t1=time(NULL);
	printf("mc done : %lds\n", t1-t0);

	return 0;
}
