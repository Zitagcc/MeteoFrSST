/********************************************************
 * Created by Yongjun ZHENG on 01 Sep 2016.             *
 * Email: yongjun.zheng@meteo.fr                        *
 * Copyright © 2016 Yongjun ZHENG. All rights reserved. *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <libxml/xmlreader.h>
#include "utils.h"

void print_time(const char *prefix, int icpu) {
    time_t current_time;
    char *c_time_string;

    /* Obtain current time. */
    current_time = time(NULL);

    if (current_time == ((time_t)-1)) {
        (void) fprintf(stderr, "Failure to obtain the current time.\n");
        exit(EXIT_FAILURE);
    }

    /* Convert to local time format. */
    c_time_string = ctime(&current_time);

    if (c_time_string == NULL) {
        (void) fprintf(stderr, "Failure to convert the current time.\n");
        exit(EXIT_FAILURE);
    }

    /* Print, ctime() has already added a terminating newline character. */
    if (icpu == 0) {
        (void) fprintf(stderr, "%s: current time is %s", prefix, c_time_string);
        //exit(EXIT_SUCCESS);
    }
}

void read_config(const char *filename, struct t_config *config) {
    xmlTextReaderPtr reader;
    int ret;
    reader = xmlReaderForFile(filename, NULL, 0);
    if (reader != NULL) {
        ret = xmlTextReaderRead(reader);
        while (ret == 1) {
            const char *name, *value;

            name = (const char *)xmlTextReaderConstName(reader);
            if (name == NULL) name = (const char *)(BAD_CAST "--");
            value = (const char *)xmlTextReaderGetAttribute(reader, (const xmlChar *)"value");

            if(!strcmp(name, "xs")) config->xs = strtod(value, NULL);
            else if(!strcmp(name, "xe")) config->xe = strtod(value, NULL);
            else if(!strcmp(name, "ys")) config->ys = strtod(value, NULL);
            else if(!strcmp(name, "ye")) config->ye = strtod(value, NULL);
            else if(!strcmp(name, "zs")) config->zs = strtod(value, NULL);
            else if(!strcmp(name, "ze")) config->ze = strtod(value, NULL);
            else if(!strcmp(name, "nx")) config->nx = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "ny")) config->ny = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "nz")) config->nz = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "ibnd_x")) config->ibnd_x = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "ibnd_y")) config->ibnd_y = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "ibnd_z")) config->ibnd_z = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "ntile_x")) config->ntile_x = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "ntile_y")) config->ntile_y = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "ntile_z")) config->ntile_z = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "whalo_x")) config->whalo_x = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "whalo_y")) config->whalo_y = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "whalo_z")) config->whalo_z = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "niter")) config->niter = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "radix")) config->radix = (int)strtol(value, NULL, 0);
            else if(!strcmp(name, "ntrunc")) config->ntrunc = (int)strtol(value, NULL, 0);

            ret = xmlTextReaderRead(reader);
        }
        xmlFreeTextReader(reader);
        if (ret != 0) fprintf(stderr, "%s : failed to parse\n", filename);

        /***********************************************************************
        ** after reading
        ***********************************************************************/
        // from degree to radian
        config->xs *= PI/180.0L;
        config->xe *= PI/180.0L;
        config->ys *= PI/180.0L;
        config->ye *= PI/180.0L;

#if defined(_CELL_)
        // resolution
        config->dx = (config->xe - config->xs) / config->nx;
        config->dy = (config->ye - config->ys) / config->ny;
        config->dz = (config->ze - config->zs) / config->nz;

        config->xs += config->dx * 0.5L;
        config->xe -= config->dx * 0.5L;
        config->ys += config->dy * 0.5L;
        config->ye -= config->dy * 0.5L;
        config->zs += config->dz * 0.5L;
        config->ze -= config->dz * 0.5L;
#else
        // resolution
        config->dx = (config->xe - config->xs) / (config->nx - 1);
        config->dy = (config->ye - config->ys) / (config->ny - 1);
        config->dz = (config->ze - config->zs) / (config->nz - 1);
#endif
    }
}

void fd_derivative(int ntrunc, char dir, double dx, struct t_domain *dm, double *fld, double *deriv) {
#if !defined(_NOMEM_)
	int i, j, k;
	int ii, jj, kk;

	for (k = dm->kts; k <= dm->kte; k++) {
		kk = k - dm->kms;
		for (j = dm->jts; j <= dm->jte; j++) {
			jj = j - dm->jms;
			for (i = dm->its; i <= dm->ite; i++) {
				ii = i - dm->ims;
                if (ntrunc == 2) {
                    if (dir == 'x') {
                        deriv[INDEX(ii, jj, kk, dm->mx, dm->mslice)] = (-fld[INDEX(ii - 1, jj, kk, dm->mx, dm->mslice)] + fld[INDEX(ii + 1, jj, kk, dm->mx, dm->mslice)]) / (2.0L * dx);
                    }else if (dir == 'y') {
                        deriv[INDEX(ii, jj, kk, dm->mx, dm->mslice)] = (-fld[INDEX(ii, jj - 1, kk, dm->mx, dm->mslice)] + fld[INDEX(ii, jj + 1, kk, dm->mx, dm->mslice)]) / (2.0L * dx);
                    }
                }else if (ntrunc == 4) {
                    if (dir == 'x') {
                        deriv[INDEX(ii, jj, kk, dm->mx, dm->mslice)] = (fld[INDEX(ii - 2, jj, kk, dm->mx, dm->mslice)] - 8*fld[INDEX(ii - 1, jj, kk, dm->mx, dm->mslice)] + 8*fld[INDEX(ii + 1, jj, kk, dm->mx, dm->mslice)] - fld[INDEX(ii + 2, jj, kk, dm->mx, dm->mslice)]) / (12.0L * dx);
                    }else if (dir == 'y') {
                        deriv[INDEX(ii, jj, kk, dm->mx, dm->mslice)] = (fld[INDEX(ii, jj - 2, kk, dm->mx, dm->mslice)] - 8*fld[INDEX(ii, jj - 1, kk, dm->mx, dm->mslice)] + 8*fld[INDEX(ii, jj + 1, kk, dm->mx, dm->mslice)] - fld[INDEX(ii, jj + 2, kk, dm->mx, dm->mslice)]) / (12.0L * dx);
                    }
                }else if (ntrunc == 6) {
                    if (dir == 'x') {
                        deriv[INDEX(ii, jj, kk, dm->mx, dm->mslice)] = (-fld[INDEX(ii - 3, jj, kk, dm->mx, dm->mslice)] + 9*fld[INDEX(ii - 2, jj, kk, dm->mx, dm->mslice)] - 45*fld[INDEX(ii - 1, jj, kk, dm->mx, dm->mslice)] + 45*fld[INDEX(ii + 1, jj, kk, dm->mx, dm->mslice)] - 9*fld[INDEX(ii + 2, jj, kk, dm->mx, dm->mslice)] + fld[INDEX(ii + 3, jj, kk, dm->mx, dm->mslice)]) / (60.0L * dx);
                    }else if (dir == 'y') {
                        deriv[INDEX(ii, jj, kk, dm->mx, dm->mslice)] = (-fld[INDEX(ii, jj - 3, kk, dm->mx, dm->mslice)] + 9*fld[INDEX(ii, jj - 2, kk, dm->mx, dm->mslice)] - 45*fld[INDEX(ii, jj - 1, kk, dm->mx, dm->mslice)] + 45*fld[INDEX(ii, jj + 1, kk, dm->mx, dm->mslice)] - 9*fld[INDEX(ii, jj + 2, kk, dm->mx, dm->mslice)] + fld[INDEX(ii, jj + 3, kk, dm->mx, dm->mslice)]) / (60.0L * dx);
                    }
                }else if (ntrunc == 8) {
                    if (dir == 'x') {
                        deriv[INDEX(ii, jj, kk, dm->mx, dm->mslice)] = (3*fld[INDEX(ii - 4, jj, kk, dm->mx, dm->mslice)] - 32*fld[INDEX(ii - 3, jj, kk, dm->mx, dm->mslice)] + 168*fld[INDEX(ii - 2, jj, kk, dm->mx, dm->mslice)] - 672*fld[INDEX(ii - 1, jj, kk, dm->mx, dm->mslice)] + 672*fld[INDEX(ii + 1, jj, kk, dm->mx, dm->mslice)] - 168*fld[INDEX(ii + 2, jj, kk, dm->mx, dm->mslice)] + 32*fld[INDEX(ii + 3, jj, kk, dm->mx, dm->mslice)] - 3*fld[INDEX(ii + 4, jj, kk, dm->mx, dm->mslice)]) / (840.0L * dx);
                    }else if (dir == 'y') {
                        deriv[INDEX(ii, jj, kk, dm->mx, dm->mslice)] = (3*fld[INDEX(ii, jj - 4, kk, dm->mx, dm->mslice)] - 32*fld[INDEX(ii, jj - 3, kk, dm->mx, dm->mslice)] + 168*fld[INDEX(ii, jj - 2, kk, dm->mx, dm->mslice)] - 672*fld[INDEX(ii, jj - 1, kk, dm->mx, dm->mslice)] + 672*fld[INDEX(ii, jj + 1, kk, dm->mx, dm->mslice)] - 168*fld[INDEX(ii, jj + 2, kk, dm->mx, dm->mslice)] + 32*fld[INDEX(ii, jj + 3, kk, dm->mx, dm->mslice)] - 3*fld[INDEX(ii, jj + 4, kk, dm->mx, dm->mslice)]) / (840.0L * dx);
                    }
                }
			}
		}
	}
#endif
}

void output_one(struct t_domain *dm, double *fld) {
    char fn[128];
    sprintf(fn, "%08d/%06dx%06d.%03d", dm->ncpu, dm->nx_g, dm->ny_g, dm->icpu);
    FILE *fp = fopen(fn, "ab");
    fwrite(fld, dm->mx * dm->my * dm->mz * sizeof(double), 1, fp);
    fclose(fp);
}

void output_all(int icpu, struct t_domain *da, double *fld) {
#if !defined(_NOMEM_)
    int ierr;
    int m, n, i, j, k;
    int is, ie, js, je, ks, ke;
    struct t_domain *dm = da + icpu;

    // gather and unpack to global array
    if (icpu == 0) {
        // allocatte global array
        m = dm->nx_g * dm->ny_g * dm->nz_g;
        double *FG = (double *)malloc(m * sizeof(double));

        // copy from local array to global array
        double *G, *L;
        for (k = dm->kts; k <= dm->kte; k++) {
            G = FG + INDEX(dm->its - dm->whalo_x, dm->jts - dm->whalo_y, k - dm->whalo_z, dm->nx_g, dm->nx_g * dm->ny_g);
            L = fld + INDEX(dm->its - dm->ims, dm->jts - dm->jms, k - dm->kms, dm->mx, dm->mslice);
            for (j = dm->jts; j <= dm->jte; j++) {
                memcpy(G, L, dm->nx * sizeof(double));
                G += dm->nx_g;
                L += dm->mx;
            }
        }

        // receive from slave processors
        for (n = 1; n < dm->ncpu; n++) {
            da++;
            m = da->ncube;
            if (m <= 0) continue;
            double *FL = (double *)malloc(m * sizeof(double));

            // receive
            ierr = MPI_Recv(FL, m, MPI_DOUBLE, n, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // unpack to global array
            is = da->its - dm->whalo_x;
            ie = da->ite - dm->whalo_x;
            js = da->jts - dm->whalo_y;
            je = da->jte - dm->whalo_y;
            ks = da->kts - dm->whalo_z;
            ke = da->kte - dm->whalo_z;
            unpack(FL, FG, is, ie, js, je, ks, ke, dm->nx_g, dm->nx_g * dm->ny_g);

            free(FL);
        }

        // write out the result
        char fn[128];
        sprintf(fn, "%08d/%06dx%06d.dat", dm->ncpu, dm->nx_g, dm->ny_g);
        FILE *fp = fopen(fn, "ab");
        fwrite(FG, dm->nx_g * dm->ny_g * dm->nz_g * sizeof(double), 1, fp);
        fclose(fp);

        free(FG);
    }else{
        // allocatte local array
        m = dm->ncube;
        if (m <= 0) return;
        double *FL = (double *)malloc(m * sizeof(double));

        // pack local array
        is = dm->its - dm->ims;
        ie = dm->ite - dm->ims;
        js = dm->jts - dm->jms;
        je = dm->jte - dm->jms;
        ks = dm->kts - dm->kms;
        ke = dm->kte - dm->kms;
        pack(fld, FL, is, ie, js, je, ks, ke, dm->mx, dm->mslice);

        // send local array to master processor
        ierr = MPI_Send(FL, m, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

        free(FL);
    }
#endif
}

void do_stat(struct t_stat *stat_avg, struct t_stat *stat_max, struct t_stat *stat_min, struct t_stat *stat, int navg) {
    // avg
    stat_avg->cal_time += stat->cal_time;
    stat_avg->mpi_time += stat->mpi_time;
    stat_avg->pack_time += stat->pack_time;
    stat_avg->unpack_time += stat->unpack_time;
    stat_avg->total_time += stat->total_time;
    stat_avg->nbyte_pack += stat->nbyte_pack;
    stat_avg->nbyte_unpack += stat->nbyte_unpack;
    stat_avg->nbyte_send += stat->nbyte_send;
    stat_avg->nbyte_recv += stat->nbyte_recv;
    if (navg) {
        stat_avg->cal_time /= navg;
        stat_avg->mpi_time /= navg;
        stat_avg->pack_time /= navg;
        stat_avg->unpack_time /= navg;
        stat_avg->total_time /= navg;
        stat_avg->nbyte_pack /= navg;
        stat_avg->nbyte_unpack /= navg;
        stat_avg->nbyte_send /= navg;
        stat_avg->nbyte_recv /= navg;
    }
    // max
    if (stat->cal_time > stat_max->cal_time) stat_max->cal_time = stat->cal_time;
    if (stat->mpi_time > stat_max->mpi_time) stat_max->mpi_time = stat->mpi_time;
    if (stat->pack_time > stat_max->pack_time) stat_max->pack_time = stat->pack_time;
    if (stat->unpack_time > stat_max->unpack_time) stat_max->unpack_time = stat->unpack_time;
    if (stat->total_time > stat_max->total_time) stat_max->total_time = stat->total_time;
    if (stat->nbyte_pack > stat_max->nbyte_pack) stat_max->nbyte_pack = stat->nbyte_pack;
    if (stat->nbyte_unpack > stat_max->nbyte_unpack) stat_max->nbyte_unpack = stat->nbyte_unpack;
    if (stat->nbyte_send > stat_max->nbyte_send) stat_max->nbyte_send = stat->nbyte_send;
    if (stat->nbyte_recv > stat_max->nbyte_recv) stat_max->nbyte_recv = stat->nbyte_recv;
    // min
    if (stat->cal_time < stat_min->cal_time) stat_min->cal_time = stat->cal_time;
    if (stat->mpi_time < stat_min->mpi_time) stat_min->mpi_time = stat->mpi_time;
    if (stat->pack_time < stat_min->pack_time) stat_min->pack_time = stat->pack_time;
    if (stat->unpack_time < stat_min->unpack_time) stat_min->unpack_time = stat->unpack_time;
    if (stat->total_time < stat_min->total_time) stat_min->total_time = stat->total_time;
    if (stat->nbyte_pack < stat_min->nbyte_pack) stat_min->nbyte_pack = stat->nbyte_pack;
    if (stat->nbyte_unpack < stat_min->nbyte_unpack) stat_min->nbyte_unpack = stat->nbyte_unpack;
    if (stat->nbyte_send < stat_min->nbyte_send) stat_min->nbyte_send = stat->nbyte_send;
    if (stat->nbyte_recv < stat_min->nbyte_recv) stat_min->nbyte_recv = stat->nbyte_recv;
}

void output_stat(int icpu, int ncpu, struct t_stat *stat, FILE *fp) {
    int ierr;

#if defined(_SHM_)
    ierr = MPI_Gather(NULL, 0, MPI_CHAR, NULL, 0, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

    if (icpu == 0) {
        struct t_stat stat_avg, stat_max, stat_min;
        stat_avg = *stat;
        stat_max = *stat;
        stat_min = *stat;

        int k;
        for (k = 1; k < ncpu; k++) {
#if !defined(_SHM_)
            ierr = MPI_Recv(stat, sizeof(struct t_stat), MPI_BYTE, MPI_ANY_SOURCE, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
            stat++;
#endif

            do_stat(&stat_avg, &stat_max, &stat_min, stat, k+1==ncpu ? ncpu : 0);
        }
#define PRIu64 "llu"
        fprintf(fp, "%15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %10" PRIu64 " %10" PRIu64 " %10" PRIu64 " %10" PRIu64 " ",  stat_avg.cal_time, stat_avg.mpi_time, stat_avg.pack_time, stat_avg.unpack_time, stat_avg.total_time, stat_avg.nbyte_pack, stat_avg.nbyte_unpack, stat_avg.nbyte_send, stat_avg.nbyte_recv);
        fprintf(fp, "%15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %10" PRIu64 " %10" PRIu64 " %10" PRIu64 " %10" PRIu64 " ",  stat_max.cal_time, stat_max.mpi_time, stat_max.pack_time, stat_max.unpack_time, stat_max.total_time, stat_max.nbyte_pack, stat_max.nbyte_unpack, stat_max.nbyte_send, stat_max.nbyte_recv);
        fprintf(fp, "%15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %10" PRIu64 " %10" PRIu64 " %10" PRIu64 " %10" PRIu64 "\n", stat_min.cal_time, stat_min.mpi_time, stat_min.pack_time, stat_min.unpack_time, stat_min.total_time, stat_min.nbyte_pack, stat_min.nbyte_unpack, stat_min.nbyte_send, stat_min.nbyte_recv);
    }else{
#if !defined(_SHM_)
        ierr = MPI_Send(stat, sizeof(struct t_stat), MPI_BYTE, 0, 999, MPI_COMM_WORLD);
#endif
    }
}
