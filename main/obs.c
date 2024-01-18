#include "trihb.h"

void ReadInfo(char *dir_data, Env *env) {
	FILE *f;
	char fn_info[BUF_SIZE], buf[BUF_SIZE];
	double M;

	sprintf(fn_info, "%s/info.txt", dir_data);
	f = fopen(fn_info, "r");

	while(!feof(f)) {
		fgets(buf, sizeof(buf), f);
		if (strstr(buf, "L")) sscanf(buf, "L %d", &env->L);
		else if(strstr(buf, "M")) {
			sscanf(buf, "M %lf", &M);
			env->M = (int)M;
		}
	}
	fclose(f);
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <dir_data>\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1);
	
	char *dir_data=argv[1];
	Env env;
	ReadInfo(dir_data, &env);
	env.N = N(env.L);

	char dir_data_par[BUF_SIZE];
	strcpy(dir_data_par, dir_data);
	dirname(dir_data_par);

	char fn_obs[BUF_SIZE], fn_script[BUF_SIZE];
	sprintf(fn_obs, "%s/obs.txt", dir_data);
	sprintf(fn_script, "%s/script.txt", dir_data_par);

	time_t t0=time(NULL);

	FILE *f_script=fopen(fn_script, "r"), *f_obs=fopen(fn_obs, "w");
	int M_tot;
	char fn_data[BUF_SIZE];
	double Mz, Rho, Ozz;
	Lat lat[env.N];
	Obs obs;

	while(fscanf(f_script, "%lf%lf%lf", &env.h, &env.D, &env.T) == 3) {
		sprintf(fn_data, "%s/h%.4f_D%.4f_T%.4f.bin", dir_data, env.h, env.D, env.T);
		ReadLatObs(env, fn_data, lat, &obs);

		M_tot = env.M + obs.M;
		Mz  = obs.mz / M_tot;	
		Rho = RHO(env.N, env.T, obs.rho1 / M_tot, obs.rho2 / M_tot);
		Ozz = obs.ozz / M_tot;

		fprintf(f_obs, "%16f%16f%16f%16f%16f%16f\n", env.h, env.D, env.T, Mz, Rho, Ozz);
	}	
	fclose(f_script);
	fclose(f_obs);

	time_t t1=time(NULL);
	printf("%s done : %lds\n", fn_obs, t1-t0);

	return 0;
}
