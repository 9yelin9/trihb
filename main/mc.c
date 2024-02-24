#include "trihb.h"

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <L> <M> <eqs> <mcs> <h> <D> <T_min> <T_max>\n", argv[0]);
		exit(1);
	}

	int m;
	Env env = {
		.L = atoi(argv[1]),
		.N = env.L * env.L,
		.M = atoi(argv[2]),
		.eqs = (int)atof(argv[3]),
		.mcs = (int)atof(argv[4]),
		.h = atof(argv[5]),
		.D = atof(argv[6]),
	};
	double T_min=atof(argv[7]), T_max=atof(argv[8]), T[env.M];
	for(m=0; m<env.M; m++) T[m] = T_min + (T_max - T_min) * m / (env.M-1);

	sprintf(env.dir_save, "data/Nt%d_Np%d_L%d_M%d_eqs%.e_mcs%.e_Tmin%.4f_Tmax%.4f", N_THETA, N_PHI, env.L, env.M, (double)env.eqs, (double)env.mcs, T_min, T_max);
	if(-access(env.dir_save, 0)) mkdir(env.dir_save, 0755);
	printf("dir_save : %s\n\n", env.dir_save);

	Lat **lat=(Lat**)malloc(sizeof(Lat*) * env.M); for(m=0; m<env.M; m++) lat[m] = (Lat*)malloc(sizeof(Lat) * env.N);
	Obs *obs=(Obs*)malloc(sizeof(Obs) * env.M);

	time_t t0=time(NULL);

	InitLat(&env, lat); //InitT(&env, lat, T);
	RunMonteCarlo(&env, lat, obs, T, env.eqs); memset(obs, 0, sizeof(Obs) * env.M);
	RunMonteCarlo(&env, lat, obs, T, env.mcs);
	SaveLatObs(&env, lat, obs, T);

	time_t t1=time(NULL);
	printf("mc done : %lds\n", t1-t0);

	for(m=0; m<env.M; m++) free(lat[m]); free(lat);
	free(obs);

	return 0;
}
