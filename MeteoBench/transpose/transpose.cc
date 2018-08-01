/********************************************************
 * Created by Yongjun ZHENG on 10 Sep 2016.             *
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
#if defined(_DUMPI_)
#include <dumpi/libdumpi/libdumpi.h>
#endif

#if defined(_NOMEM_)
#undef _COMPUTATION_
#undef _OUTPUT_
#endif

#if defined(_ALLTOALLV_) || defined(_ALLTOALL_)
#define _COMM_
#endif

#define sstmac_app_name transpose

int main(int argc, char *argv[]) {
#if defined(_SHM_)
    extern char *GO_MBLOCK;
#else
    char *GO_MBLOCK;
#endif

    int COMM_OFFSET = 0;
    int ierr, ncpu, icpu;
    double t_all, t_init, t_transpose;

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
#if defined(_DUMPI_)
    libdumpi_disable_profiling();
#endif

    // print the current time before collective communication
    print_time("Before collective communication", icpu);
    t_all = t_init = mysecond();

#if defined(_SHM_)
    //int shmid = atoi(getenv("GO_MBLOCK_SHMID"));
    //shm_malloc((void **)(&GO_MBLOCK), shmid);
    if (icpu == 0) fprintf(stderr, "Shared memory seen from virtual MPI: %p\nMPI_PROC_NULL=%d\n", GO_MBLOCK, MPI_PROC_NULL);
#else
    //GO_MBLOCK = (char *)malloc(sizeof(struct t_config) + sizeof(struct t_domain) * ncpu * 2);
    shmc((void **)(&GO_MBLOCK), sizeof(struct t_config) + sizeof(struct t_domain) * ncpu * 2 + sizeof(struct t_stat) * ncpu, 0);
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
    struct t_domain *dm_y = NULL, *dm_x = NULL;
    struct t_domain *da_y = NULL, *da_x = NULL;
    da_y = (struct t_domain *)((char *)GO_MBLOCK + sizeof(struct t_config));
    da_x = (struct t_domain *)((char *)GO_MBLOCK + sizeof(struct t_config) + sizeof(struct t_domain) * ncpu);
#if defined(_EXTERNAL_)
    dm_y = da_y + icpu;
    dm_x = da_x + icpu;
#else
    if (icpu == 0) {
        for (icpu = 0; icpu < ncpu; icpu++) {
            dm_y = da_y + icpu;
            dm_x = da_x + icpu;
            // 3.1 domain partition
            dm_y->ntile_x = config->ntile_x, dm_y->ntile_z = config->ntile_z;
            dm_x->ntile_y = config->ntile_y, dm_x->ntile_z = config->ntile_z;
            domain_partition('Y', icpu, ncpu, config, dm_y);
            dm_x->ntile_y = dm_y->ntile_x, dm_x->ntile_z = dm_y->ntile_z;
            domain_partition('X', icpu, ncpu, config, dm_x);
        }
        icpu = 0;
    }
    dm_y = da_y + icpu;
    dm_x = da_x + icpu;
    // broadcast
#if defined(_COLL_)
    ierr = bcast(da_y, sizeof(struct t_domain) * ncpu * 2, MPI_BYTE, MPI_COMM_WORLD);
#else
    ierr = MPI_Bcast(da_y, sizeof(struct t_domain) * ncpu * 2, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    // print the current time after bcast partition
    print_time("After bcast partition", icpu);
#endif
    if (0 && icpu==0) { // turn on if want to check shared memory
        FILE *fid = fopen("mblock.bin", "wb");
        fwrite(config, sizeof(struct t_config), 1, fid);
        fwrite(da_y, sizeof(struct t_domain), ncpu, fid);
        fwrite(da_x, sizeof(struct t_domain), ncpu, fid);
        fclose(fid);
    }

    // 3.1 split communicator
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;
#if !defined(_ALLTOALLV_) && !defined(_ALLTOALL_)
    COMM_OFFSET = (icpu / dm_y->ntile_x) * dm_y->ntile_x;
#endif
#if defined(_COMM_)
    //int grank[1][3];
    //grank[0][0] = COMM_OFFSET;
    //grank[0][1] = grank[0][0] + dm_y->ntile_x - 1;
    //grank[0][2] = 1;
    //MPI_Group ggrp,lgrp;
    //MPI_Comm_group(MPI_COMM_WORLD, &ggrp);
    //MPI_Group_range_incl(ggrp, 1, grank, &lgrp);
    //MPI_Comm_create_group(MPI_COMM_WORLD, lgrp, 0, &comm);
    MPI_Comm_split(MPI_COMM_WORLD, icpu / dm_y->ntile_x, icpu, &comm);
#endif

    /***************************************************************************
     * 4 transpose
     **************************************************************************/
    // 4.0 allocate memory
    int iflg = 1;
    double *bufr_x = NULL;
    int *count_offset_x = NULL;
    int max_x, nsend_x, nrecv_x;
#if defined(_RINGALL_) || defined(_RINGK_) || defined(_INDEX_) || defined(_ALLTOALL_)
    partition_transpose_workmem_ring(comm, iflg, da_y+COMM_OFFSET, da_x+COMM_OFFSET, &bufr_x, &max_x);
#elif defined(_ALLTOALLV_)
    partition_transpose_workmem_alltoallv(comm, iflg, dm_y, dm_x, &bufr_x, &count_offset_x);
#else
    partition_transpose_workmem(comm, iflg, da_y+COMM_OFFSET, da_x+COMM_OFFSET, &bufr_x, &max_x, &nsend_x, &nrecv_x);
#endif

    double *fldy, *fldy_y, *fldx, *fldx_x;
#if defined(_NOMEM_)
    fldy = fldy_y = fldx = fldx_x = NULL;
#else
    double *dblock;
    if (dm_y->mcube * 2 + dm_x->mcube * 2 > 0) {
        dblock = (double *)malloc(sizeof(double) * (dm_y->mcube * 2 + dm_x->mcube * 2));
        // for partition in xz
        fldy = dblock;
        fldy_y = dblock + dm_y->mcube;
        // for partition in yz
        fldx = dblock + dm_y->mcube * 2;
        fldx_x = dblock + dm_y->mcube * 2 + dm_x->mcube;

        if (fldy == NULL || fldy_y == NULL || fldx == NULL || fldx_x == NULL) {
            fprintf(stderr, "Program terminated: can't allocate memory for fldy, fldy_y, fldx, or fldx_x!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
#endif

    // 4.1 usually we scatter the global field to each processor
    //     but here we assign the local field analytically
    int i, j, k;
#if !defined(_NOMEM_)
    for (j = dm_y->jts; j <= dm_y->jte; j++) {
        double y = config->ys + (j - config->whalo_y) * config->dy;
        for (i = dm_y->its; i <= dm_y->ite; i++) {
            double x = config->xs + (i - config->whalo_x) * config->dx;
            double f = sin(x) * cos(2.0L*y);
            for (k = dm_y->kts; k <= dm_y->kte; k++) fldy[INDEX(i - dm_y->ims, j - dm_y->jms, k - dm_y->kms, dm_y->mx, dm_y->mslice)] = f;
        }
    }
#endif

    // 4.2 tranpose the partition for calculation of the derivatives by spectral method
    double tavg_tran = 0.0L;
    double tmin_tran = 1.0E9;
    double tmax_tran = 0.0L;
    struct t_stat *stat = (struct t_stat*)((char *)GO_MBLOCK + sizeof(struct t_config) + sizeof(struct t_domain) * ncpu *2 + sizeof(struct t_stat) * icpu);
    double t, tt;
    int iter;

    char fn[128];

    FILE *fp;
    if (icpu == 0) {
        sprintf(fn, "%08d", ncpu);
        // create directory
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

#if defined(_DUMPI_)
    libdumpi_enable_profiling();
#endif
    int itag;
    for (iter = 1; iter <= config->niter; iter++) {
        // transpose xz -> yz
        itag = 1;
        tt = mysecond();
#if defined(_RINGALL_)
        strcpy(fn, "transpose_ringall");
        t = partition_transpose_ringall(comm, itag, da_y+COMM_OFFSET, da_x+COMM_OFFSET, fldy, fldx, stat, bufr_x, max_x);
#elif defined(_RINGK_)
        strcpy(fn, "transpose_ringone");
        t = partition_transpose_ringone(comm, itag, da_y+COMM_OFFSET, da_x+COMM_OFFSET, fldy, fldx, stat, bufr_x, max_x, config->radix);
#elif defined(_INDEX_)
        strcpy(fn, "transpose_index");
        t = partition_transpose_index(comm, itag, da_y+COMM_OFFSET, da_x+COMM_OFFSET, fldy, fldx, stat, bufr_x, max_x, config->radix);
#elif defined(_ALLTOALL_)
        strcpy(fn, "transpose_alltoall");
        t = partition_transpose_alltoall(comm, itag, da_y+COMM_OFFSET, da_x+COMM_OFFSET, fldy, fldx, stat, bufr_x, max_x);
#elif defined(_ALLTOALLV_)
        strcpy(fn, "transpose_alltoallv");
        t = partition_transpose_alltoallv(comm, itag, da_y, da_x, fldy, fldx, stat, bufr_x, count_offset_x);
#else
        strcpy(fn, "transpose_default");
        t = partition_transpose(comm, itag, da_y+COMM_OFFSET, da_x+COMM_OFFSET, fldy, fldx, stat, bufr_x, max_x, nsend_x, nrecv_x);
#endif
        t_transpose += (tt = mysecond() - tt);

        tavg_tran += t;
        tmin_tran = tmin_tran > t ? t : tmin_tran;
        tmax_tran = tmax_tran < t ? t : tmax_tran;
        if (icpu == 0) fprintf(stderr, "T used %15lf seconds for %6d: %15lf %15lf %15lf using %s\n", t_transpose, iter, tmin_tran, tavg_tran/iter, tmax_tran, fn);

#if defined(_COMPUTATION_)
        t = MPI_Wtime();
        fft_derivative('x', dm_x, fldx, fldx_x, 2.0L*PI);
        t = MPI_Wtime() - t;
#else
        t = 0.0L;
#endif
        stat->cal_time = t;
    }
#if defined(_DUMPI_)
    libdumpi_disable_profiling();
#endif
#if defined(_STAT_)
    output_stat(icpu, ncpu, stat, fp);
#endif
    if (icpu == 0) fclose(fp);

    // 4.3 gather from each processor to the global field and output the calculated results
#if defined(_OUTPUT_)
    output_one(dm_y, fldy);
    output_all(icpu, da_y, fldy);
    output_all(icpu, da_x, fldx);
    // Fx
    output_all(icpu, da_x, fldx_x);
    // output analytical derivative in x direction
    for (j = dm_x->jts; j <= dm_x->jte; j++) {
        double y = config->ys + (j - config->whalo_y) * config->dy;
        for (i = dm_x->its; i <= dm_x->ite; i++) {
            double x = config->xs + (i - config->whalo_x) * config->dx;
            double fx = cos(x) * cos(2.0L*y);
            for (k = dm_x->kts; k <= dm_x->kte; k++) fldx_x[INDEX(i - dm_x->ims, j - dm_x->jms, k - dm_x->kms, dm_x->mx, dm_x->mslice)] = fx;
        }
    }
    output_all(icpu, da_x, fldx_x);
    // Fy
    fft_derivative('y', dm_y, fldy, fldy_y, PI);
    output_all(icpu, da_y, fldy_y);
    // output analytical derivative in y direction
    for (j = dm_y->jts; j <= dm_y->jte; j++) {
        double y = config->ys + (j - config->whalo_y) * config->dy;
        for (i = dm_y->its; i <= dm_y->ite; i++) {
            double x = config->xs + (i - config->whalo_x) * config->dx;
            double fy = -2.0L * sin(x) * sin(2.0L*y);
            for (k = dm_y->kts; k <= dm_y->kte; k++) fldy_y[INDEX(i - dm_y->ims, j - dm_y->jms, k - dm_y->kms, dm_y->mx, dm_y->mslice)] = fy;
        }
    }
    output_all(icpu, da_y, fldy_y);
#endif

    // 4.4 deallocate the memory
    iflg = 0;
#if defined(_RINGALL_) || defined(_RINGK_) || defined(_INDEX_) || defined(_ALLTOALL_)
    partition_transpose_workmem_ring(comm, iflg, da_y+COMM_OFFSET, da_x+COMM_OFFSET, &bufr_x, &max_x);
#elif defined(_ALLTOALLV_)
    partition_transpose_workmem_alltoallv(comm, iflg, dm_y, dm_x, &bufr_x, &count_offset_x);
#else
    partition_transpose_workmem(comm, iflg, da_y+COMM_OFFSET, da_x+COMM_OFFSET, &bufr_x, &max_x, &nsend_x, &nrecv_x);
#endif
#if defined(_NOMEM_)
    fldy = fldy_y = fldx = fldx_x = NULL;
#else
    fldy = fldy_y = fldx = fldx_x = NULL;
    if (dm_y->mcube * 2+ dm_x->mcube * 2 > 0) free(dblock);
#endif
    config = NULL;
    dm_y = dm_x = NULL;
    da_y = da_x = NULL;
#if defined(_SHM_)
    //shm_free((void **)(&GO_MBLOCK));
#else
    shmd((void **)(&GO_MBLOCK), sizeof(struct t_config) + sizeof(struct t_domain) * ncpu * 2 + sizeof(struct t_stat) * ncpu, 0);
#endif

#if defined(_COMM_)
    //MPI_Group_free(&ggrp);
    //MPI_Group_free(&lgrp);
    MPI_Comm_free(&comm);
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
#if defined(_DUMPI_)
    libdumpi_enable_profiling();
#endif
    ierr = MPI_Finalize();
    if (ierr != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Finalize failed: %d\n", ierr);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }else{
        return 0;
    }
}
