/********************************************************
 * Created by Yongjun ZHENG on 01 Nov 2016.             *
 * Email: yongjun.zheng@meteo.fr                        *
 * Copyright © 2016 Yongjun ZHENG. All rights reserved. *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include "utils.h"
#include "coll.h"
#include "shm.h"

#if defined(_NOMEM_)
#undef _COMPUTATION_
#undef _OUTPUT_
#endif

#define sstmac_app_name halo

int main(int argc, char *argv[]) {
#if defined(_SHM_)
  extern char *GO_MBLOCK;
#else
  char *GO_MBLOCK;
#endif
	int ierr, ncpu, icpu;
	double t_all, t_init, t_halo = 0.0L;
	double tt;

	/***************************************************************************
	 * 1 MPI initialization
	 **************************************************************************/
	ierr = MPI_Init(&argc, &argv);
	if (ierr != MPI_SUCCESS) {
		fprintf(stderr, "MPI_Init failed: %d\n", ierr);
		return -1;
	}
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
	if (ierr != MPI_SUCCESS) {
		fprintf(stderr, "MPI_Comm_size failed: %d\n", ierr);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &icpu);
	if (ierr != MPI_SUCCESS) {
		fprintf(stderr, "MPI_Comm_rank failed: %d\n", ierr);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// print the current time before collective communication
  print_time("Before collective communication", icpu);
  t_all = t_init = mysecond();

#if defined(_SHM_)
	//int shmid = atoi(getenv("GO_MBLOCK_SHMID"));
	//shm_malloc((void **)(&GO_MBLOCK), shmid);
	if (icpu == 0) fprintf(stderr, "Shared memory seen from virtual MPI: %p\nMPI_PROC_NULL=%d\n", GO_MBLOCK, MPI_PROC_NULL);
#else
	GO_MBLOCK = (char *)SMPI_SHARED_MALLOC(sizeof(struct t_config) + sizeof(struct t_domain) * ncpu);
	// shmc((void **)(&GO_MBLOCK), sizeof(struct t_config) + sizeof(struct t_domain) * ncpu + sizeof(struct t_stat) * ncpu, 0);
#endif

	/***************************************************************************
	 * 2 Configuration
	 **************************************************************************/
	// read in the configuration
	struct t_config *config;
	config = (struct t_config *)GO_MBLOCK;
#if !defined(_EXTERNAL_)
	if (icpu == 0) read_config("./config.xml", config);
	// broadcast
#if defined(_COLL_)
	ierr = bcast(config, sizeof(struct t_config), MPI_BYTE, MPI_COMM_WORLD);
#else
	ierr = MPI_Bcast(config, sizeof(struct t_config), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
	// print the current time after bcast config
	print_time("After bcast config", icpu);
#endif

#if defined(_WEAKSCALE_)
int nx, ny, nz;
for (nz = 32; nz <= 128; nz += 32) {
for (ny = 900; ny <= 3600; ny += 900) {
for (nx = 1800; nx <= 7200; nx += 1800) {
  config->nx = nx;
  config->ny = ny;
  config->nz = nz;
#endif

	/***************************************************************************
	 * 3 domain partition
	 **************************************************************************/
	struct t_domain *dm_z = NULL;
	struct t_domain *da_z = NULL;
	da_z = (struct t_domain *)((char *)GO_MBLOCK + sizeof(struct t_config));
#if defined(_EXTERNAL_)
	dm_z = da_z + icpu;
#else
	if (icpu == 0) {
		for (icpu = 0; icpu < ncpu; icpu++) {
			dm_z = da_z + icpu;
			// 3.1 domain partition
			dm_z->ntile_x = config->ntile_x, dm_z->ntile_y = config->ntile_y;
			domain_partition('Z', icpu, ncpu, config, dm_z);
		}
		icpu = 0;
	}
	dm_z = da_z + icpu;
	// broadcast
#if defined(_COLL_)
	ierr = bcast(da_z, sizeof(struct t_domain) * ncpu, MPI_BYTE, MPI_COMM_WORLD);
#else
	ierr = MPI_Bcast(da_z, sizeof(struct t_domain) * ncpu, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
	// print the current time after bcast partition
	print_time("After bcast partition", icpu);
#endif

	/***************************************************************************
	 * 4 test halo exchange and transpose
	 **************************************************************************/
	// 4.0 allocate memory
	double *bufr_halo = NULL;
	halo_exchange_workmem(1, dm_z, &bufr_halo);

	double *fldz, *fldz_x, *fldz_y;
#if defined(_NOMEM_)
    fldz = fldz_x = fldz_y = NULL;
#else
	double *dblock;
	if (dm_z->mcube > 0) {
		dblock = (double *)SMPI_SHARED_MALLOC(sizeof(double) * (dm_z->mcube * 3));
		// field
		fldz = dblock;
		// derivative in x direction
		fldz_x = dblock + dm_z->mcube;
		// derivative in x direction
		fldz_y = dblock + dm_z->mcube * 2;

		if (fldz == NULL || fldz_x == NULL || fldz_y == NULL) {
			fprintf(stderr, "Program terminated: can't allocate memory for fldz, fldz_x, or fldz_y!\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
#endif

	// 4.1 usually we scatter the global field to each processor
	//     but here we assign the local field analytically
	int i, j, k;
#if !defined(_NOMEM_)
	for (j = dm_z->jts; j <= dm_z->jte; j++) {
		double y = config->ys + (j - config->whalo_y) * config->dy;
		for (i = dm_z->its; i <= dm_z->ite; i++) {
			double x = config->xs + (i - config->whalo_x) * config->dx;
			double f = sin(x) * cos(2.0L*y);
			for (k = dm_z->kts; k <= dm_z->kte; k++) fldz[INDEX(i - dm_z->ims, j - dm_z->jms, k - dm_z->kms, dm_z->mx, dm_z->mslice)] = f;
		}
	}
#endif

	// 4.2 exchange halo region for calculation of the derivatives by finite difference method
	//     and tranpose the partition for calculation of the derivatives by spectral method
	double tavg_halo = 0.0L;
	double tmin_halo = 1.0E9;
	double tmax_halo = 0.0L;
	struct t_stat *stat = (struct t_stat*)((char *)GO_MBLOCK + sizeof(struct t_config) + sizeof(struct t_domain) * ncpu + sizeof(struct t_stat) * icpu);
	double t;
	int iter;

	FILE *fp;
	char fn[128];
	if (icpu == 0) {
		// create directory
		sprintf(fn, "%08d", ncpu);
		DIR *dp = opendir(fn);
		if (dp) {
	    closedir(dp);
		} else if (ENOENT == errno) {
	    mkdir(fn, 0750);
		} else {
	    fprintf(stderr, "Program terminated: can't open directory %s\n", fn);
	    MPI_Abort(MPI_COMM_WORLD, 1);
		}

		// create file
		sprintf(fn, "%08d/%06dx%06d.stat", ncpu, config->nx, config->ny);
		fp = fopen(fn, "wt");
		if (fp==NULL) {
			fprintf(stderr, "Program terminated: can't open file %s\n", fn);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	// print the current time
	t_init = mysecond() - t_init;
	if (icpu == 0) fprintf(stderr, "Initialized elapsed time: %lf seconds\n", t_init);

	for (iter = 1; iter <= config->niter; iter++) {
		// exchange halo region
		tt = mysecond();
		t = halo_exchange(icpu, da_z, fldz, stat, bufr_halo);
		t_halo += mysecond() - tt;

		tavg_halo += t;
		tmin_halo = tmin_halo > t ? t : tmin_halo;
		tmax_halo = tmax_halo < t ? t : tmax_halo;
		if (icpu == 0) fprintf(stderr, "H used %15lf seconds for %6d: %15lf %15lf %15lf\n", t_halo, iter, tmin_halo, tavg_halo/iter, tmax_halo);

		// computation by finite difference method
#if defined(_COMPUTATION_)
		t = MPI_Wtime();
		fd_derivative(config->ntrunc, 'x', config->dx, dm_z, fldz, fldz_x);
		fd_derivative(config->ntrunc, 'y', config->dy, dm_z, fldz, fldz_y);
		t = MPI_Wtime() - t;
#else
		t = 0.0L;
#endif
		stat->cal_time = t;
	}
#if defined(_STAT_)
	output_stat(icpu, ncpu, stat, fp);
#endif
	if (icpu == 0) fclose(fp);

	// 4.3 gather from each processor to the global field and output the calculated results
#if defined(_OUTPUT_)
	output_one(dm_z, fldz);
	output_all(icpu, da_z, fldz);
	output_all(icpu, da_z, fldz_x);
	output_all(icpu, da_z, fldz_y);
	// output analytical derivative in x direction
	for (j = dm_z->jts; j <= dm_z->jte; j++) {
		double y = config->ys + (j - config->whalo_y) * config->dy;
		for (i = dm_z->its; i <= dm_z->ite; i++) {
			double x = config->xs + (i - config->whalo_x) * config->dx;
			double fx = cos(x) * cos(2.0L*y);
			for (k = dm_z->kts; k <= dm_z->kte; k++) fldz_x[INDEX(i - dm_z->ims, j - dm_z->jms, k - dm_z->kms, dm_z->mx, dm_z->mslice)] = fx;
		}
	}
	output_all(icpu, da_z, fldz_x);
	// output analytical derivative in y direction
	for (j = dm_z->jts; j <= dm_z->jte; j++) {
		double y = config->ys + (j - config->whalo_y) * config->dy;
		for (i = dm_z->its; i <= dm_z->ite; i++) {
			double x = config->xs + (i - config->whalo_x) * config->dx;
			double fy = sin(x) * sin(2.0L*y) * (-2.0L);
			for (k = dm_z->kts; k <= dm_z->kte; k++) fldz_y[INDEX(i - dm_z->ims, j - dm_z->jms, k - dm_z->kms, dm_z->mx, dm_z->mslice)] = fy;
		}
	}
	output_all(icpu, da_z, fldz_y);
#endif

	// 4.4 deallocate the memory
	halo_exchange_workmem(0, dm_z, &bufr_halo);
#if defined(_NOMEM_)
	fldz = fldz_x = fldz_y = NULL;
#else
	fldz = fldz_x = fldz_y = NULL;
	if (dm_z->mcube > 0) SMPI_SHARED_FREE(dblock);
#endif
	config = NULL;
	dm_z = NULL;
	da_z = NULL;
#if defined(_SHM_)
	//shm_free((void **)(&GO_MBLOCK));
#else
	SMPI_SHARED_FREE((void **)(&GO_MBLOCK));
#endif


#if defined(_WEAKSCALE_)
}
}
}
#endif

	t_all = mysecond() - t_all;
	if (icpu == 0) {
		fprintf(stderr, "Total elapsed time: %lf seconds\n", t_all);
		fprintf(stderr, "****** %d TASKS SUCCEED ******\n", ncpu);
	}

	/***************************************************************************
	 * 5 MPI finalization
	 **************************************************************************/
	ierr = MPI_Finalize();
	if (ierr != MPI_SUCCESS) {
		fprintf(stderr, "MPI_Finalize failed: %d\n", ierr);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}else{
		return 0;
	}
}
