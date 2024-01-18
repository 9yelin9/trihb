#include "trihb.h"

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <dir_init> <init_mode=eq,mc> <dir_save> <L> <M> <h> <D> <T>\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1);

	char *dir_init=argv[1], *init_mode=argv[2], *dir_save=argv[3];
	Env env = {
		.L = atoi(argv[4]),
		.N = N(env.L),
		.M = (int)atof(argv[5]),
		.h = atof(argv[6]),
		.D = atof(argv[7]),
		.T = atof(argv[8]),
	};
	Lat lat[env.N];
	Obs obs={0};

	char fn[BUF_SIZE];
	sprintf(fn, "h%.4f_D%.4f_T%.4f.bin", env.h, env.D, env.T);

	if(strstr(dir_init, "none")) InitLat(env, lat);
	else {
		char fn_init[BUF_SIZE];
		sprintf(fn_init, "%s/%s", dir_init, fn);
		ReadLatObs(env, fn_init, lat, &obs);
		if(strstr(init_mode, "eq")) memset(&obs, 0, sizeof(obs));
	}

	time_t t0=time(NULL);

	RunMonteCarlo(env, lat, &obs);

	char fn_save[BUF_SIZE];
	sprintf(fn_save, "%s/%s", dir_save, fn);
	WriteLatObs(env, fn_save, lat, &obs);
	//printf("%d\t%f\t%f\t%f\n", env.M, obs.mz/env.M, RHO(env.N, env.T, obs.rho1/env.M, obs.rho2/env.M), obs.ozz/env.M);

	time_t t1=time(NULL);
	printf("%s done : %lds\n", fn_save, t1-t0);

	return 0;
}
