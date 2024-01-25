#include "trihb.h"

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <dir_data> <init_mode=mc,eq> <Meq> <Mmc> <L> <h> <D> <T>\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1);

	char *dir_data=argv[1], *init_mode=argv[2];
	Env env = {
		.Meq = (int)atof(argv[3]),
		.Mmc = (int)atof(argv[4]),
		.L = atoi(argv[5]),
		.N = env.L * env.L,
		.h = atof(argv[6]),
		.D = atof(argv[7]),
		.T = atof(argv[8]),
	};
	Lat lat[env.N];
	Obs obs_eq={0}, obs_mc={0};

	char dir_save[BUF_SIZE], fn[BUF_SIZE], path_init[BUF_SIZE], path_save[BUF_SIZE];
	sprintf(fn, "h%.4f_D%.4f_T%.4f.h5", env.h, env.D, env.T);
	sprintf(path_init, "data/%s/L%d_Mmc%.e/%s", dir_data, env.L, (double)env.Meq, fn);

	time_t t0=time(NULL);

	if(strstr(init_mode, "mc")) {
		if(env.Meq == 0) {
			InitLat(env, lat);
			RunMonteCarlo(env, lat, &obs_mc);
		}
		else {
			ReadLatObs(env, lat, &obs_mc, path_init);
			RunMonteCarlo(env, lat, &obs_mc);
		}

		env.Mmc += env.Meq;
		sprintf(dir_save, "data/%s/L%d_Mmc%.e", dir_data, env.L, (double)env.Mmc);
		if(-access(dir_save, 0)) mkdir(dir_save, 0755);

		sprintf(path_save, "%s/%s", dir_save, fn);
		SaveLatObs(env, lat, &obs_mc, path_init);
	}
	else if(strstr(init_mode, "eq")) {
		ReadLatObs(env, lat, &obs_eq, path_init);
		RunMonteCarlo(env, lat, &obs_mc);

		sprintf(dir_save, "data/%s/L%d_Meq%.e_Mmc%.e", dir_data, env.L, (double)env.Meq, (double)env.Mmc);
		if(-access(dir_save, 0)) mkdir(dir_save, 0755);

		sprintf(path_save, "%s/%s", dir_save, fn);
		SaveLatObs(env, lat, &obs_mc, path_init);
	}
	else {
		printf("\"%s\" is wrong init_mode\n", init_mode);
		exit(1);
	}

	time_t t1=time(NULL);
	printf("%s done : %lds\n", path_save, t1-t0);

	return 0;
}
