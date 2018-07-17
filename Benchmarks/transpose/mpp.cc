/********************************************************
 * Created by Yongjun ZHENG on 08 Sep 2016.             *
 * Email: yongjun.zheng@meteo.fr                        *
 * Copyright © 2016 Yongjun ZHENG. All rights reserved. *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "coll.h"
#include "mpp.h"
#include "shm.h"

int pack(double *src, double *dst, int is, int ie, int js, int je, int ks, int ke, int ncol, int nslice) {
    int i, j, k;
    int nx, nbyte, nbyte_total = 0;
    double *slice;

    nx = ie - is + 1;
    nbyte = nx * sizeof(double);
    src += INDEX(is, js, ks, ncol, nslice); // start address

    // pack
    for (k = ks; k <= ke; k++) {
        slice = src;          // start address of this slice
        for (j = js; j <= je; j++) {
            nbyte_total += nbyte;
            memcpy(dst, src, nbyte);
            src += ncol;      // go to next row
            dst += nx;        // go to next row
        }
        src = slice + nslice; // go to next slice
    }

    return nbyte_total;
}

int unpack(double *src, double *dst, int is, int ie, int js, int je, int ks, int ke, int ncol, int nslice) {
    int i, j, k;
    int nx, nbyte, nbyte_total = 0;
    double *slice;

    nx = ie - is + 1;
    nbyte = nx * sizeof(double);
    dst += INDEX(is, js, ks, ncol, nslice); // start address

    // unpack
    for (k = ks; k <= ke; k++) {
        slice = dst;          // start address of this slice
        for (j = js; j <= je; j++) {
            nbyte_total += nbyte;
            memcpy(dst, src, nbyte);
            src += nx;        // go to next row
            dst += ncol;      // go to next row
        }
        dst = slice + nslice; // go to next slice
    }

    return nbyte_total;
}

double mysecond() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1.0E-6;
#if defined(_NANO_)
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1.0E-9;
#endif
}

int prime_factors(int num, int m, int *factors) {
    int n = num;
    int i, j;

    // first check whether 2 is its prime factor
    i = 0;
    while (n % 2 == 0) {
        factors[i++] = 2;
        if (i >= m) {
            fprintf(stderr, "The number of prime factors of %d is greater than the size (=%d) of array of factors\n", num, m);
            return -1;
        }else{
            n /= 2;
        }
    }

    // find all odd prime factors
    j = 3;
    while (n > 1) {
        while (n % j == 0) {
            factors[i++] = j;
            if (i >= m) {
                fprintf(stderr, "The number of prime factors of %d is greater than the size (=%d) of array of factors\n", num, m);
                return -1;
            }else{
                n /= j;
            }
        }
        j += 2;
    }

    // if num is a prime number
    if (i == 0) factors[i++] = num;

    // return the number of prime factors
    return i;
}

void tile_setup_2d(int ncpu, int nx, int ny, int *ntile_x, int *ntile_y) {
    if (*ntile_x == 0 && *ntile_y == 0) {
        // aspect ratio of the grid
        double aspect_grid = (double)nx / ny;
        // find the prime factors of ncpu
        int factors[64];
        int nfactors = prime_factors(ncpu, 64, factors);

        double aspect, aspect_tile = 1.0E9, aspect_tile_opt = 1.0E9;
        int i, j, k, mx, my, ntile_x_opt = 0, ntile_y_opt = 0;
        // partition so that aspect ratio of tile is closed to that of grid as far as possible
        for (k = 1; k <= nfactors; k++) {
            for (j = 0; j <= nfactors - k; j++) {
                for (i = j, mx = 1; i < j + k; i++) {
                    mx *= factors[i];
                }
                my = ncpu / mx;
                aspect = (double)mx / my;

                // best closed to the aspect ratio of grid but mx <= nx and my <= ny is not guarantteed
                if ( fabs(aspect - aspect_grid) < fabs(aspect_tile - aspect_grid)) {
                    *ntile_x = mx;
                    *ntile_y = my;
                    aspect_tile = aspect;
                }
                // closed to the aspect ratio of grid also mx <= nx and my <= ny is satified
                if ( fabs(aspect - aspect_grid) < fabs(aspect_tile_opt - aspect_grid) &&  mx <= nx && my <= ny) {
                    ntile_x_opt = mx;
                    ntile_y_opt = my;
                    aspect_tile_opt = aspect;
                }
            }
        }
        // use optimal partition if found one
        if (ntile_x_opt > 0 && ntile_y_opt > 0) {
            *ntile_x = ntile_x_opt;
            *ntile_y = ntile_y_opt;
        }
    }else if (*ntile_x == 0) {
        *ntile_x = 1;
        *ntile_y = ncpu;
    }else if (*ntile_y == 0) {
        *ntile_x = ncpu;
        *ntile_y = 1;
    }else if (*ntile_x * *ntile_y != ncpu) {
        fprintf(stderr, "Please properly set ntile_x and ntile_y so that their product is equal to the number of CPUs: %d * %d <> %d\n", *ntile_x, *ntile_y, ncpu);
    }
}

void domain_partition_2d(int ncpu, int icpu, int nx, int ny, int ibnd_x, int ibnd_y, int ntile_x, int ntile_y, int whalo_x, int whalo_y, int *its, int *ite, int *jts, int *jte, int *w_neighbor, int *e_neighbor, int *s_neighbor, int *n_neighbor) {
    int i, j, k;
    // x partition
    int ix[ntile_x];
    for (i = 0, j = nx / ntile_x, k = nx % ntile_x; i < ntile_x; i++) {
        if (i < k) {
            ix[i] = j + 1;
        }else{
            ix[i] = j;
        }
    }
    // y partition
    int iy[ntile_y];
    for (i = 0, j = ny / ntile_y, k = ny % ntile_y; i < ntile_y; i++) {
        if (i < k) {
            iy[i] = j + 1;
        }else{
            iy[i] = j;
        }
    }

    // index of my tile
    i = icpu % ntile_x;
    j = icpu / ntile_x;

    // my domain
    for (k = 0, *its = whalo_x; k < i; k++) *its += ix[k];
    *ite = *its + ix[i] - 1;
    for (k = 0, *jts = whalo_y; k < j; k++) *jts += iy[k];
    *jte = *jts + iy[j] - 1;

    // neighbors in x direction
    if (ibnd_x == 0) {
        // specified boundary condition
        if (i == 0) {
            *w_neighbor = MPI_PROC_NULL;
            if (ntile_x > 1) {
                *e_neighbor = icpu + 1;
            }else{
                *e_neighbor = MPI_PROC_NULL;
            }
        }else if (i + 1 == ntile_x) {
            if (ntile_x > 1) {
                *w_neighbor = icpu - 1;
            }else{
                *w_neighbor = MPI_PROC_NULL;
            }
            *e_neighbor = MPI_PROC_NULL;
        }else{
            *w_neighbor = icpu - 1;
            *e_neighbor = icpu + 1;
        }
    }else if (ibnd_x == 1) {
        // periodic boundary condition
        *w_neighbor = j * ntile_x + (i - 1 + ntile_x) % ntile_x;
        *e_neighbor = j * ntile_x + (i + 1 + ntile_x) % ntile_x;
    }

    // neighbors in y direction
    if (ibnd_y == 0) {
        // specified boundary condition
        if (j == 0) {
            *s_neighbor = MPI_PROC_NULL;
            if (ntile_y > 1) {
                *n_neighbor = (j + 1) * ntile_x + i;
            }else{
                *n_neighbor = MPI_PROC_NULL;
            }
        }else if (j + 1 == ntile_y) {
            if (ntile_y > 1) {
                *s_neighbor = (j - 1) * ntile_x + i;
            }else{
                *s_neighbor = MPI_PROC_NULL;
            }
            *n_neighbor = MPI_PROC_NULL;
        }else{
            *s_neighbor = (j - 1) * ntile_x + i;
            *n_neighbor = (j + 1) * ntile_x + i;
        }
    }else if (ibnd_y == 1) {
        // periodic boundary condition
        *s_neighbor = ((j - 1 + ntile_y) % ntile_y) * ntile_x + i;
        *n_neighbor = ((j + 1 + ntile_y) % ntile_y) * ntile_x + i;
    }
}


void domain_partition(char dir, int icpu, int ncpu, struct t_config *config, struct t_domain *dm) {
    dm->icpu = icpu;
    dm->ncpu = ncpu;
    dm->nx_g = config->nx;
    dm->ny_g = config->ny;
    dm->nz_g = config->nz;
    dm->whalo_x = config->whalo_x;
    dm->whalo_y = config->whalo_y;
    dm->whalo_z = config->whalo_z;
    dm->dir = dir;

    if (dir == 'Z') {
        // no partition in z direction
        dm->ntile_z = 1;
        if (dm->ntile_x == 0 || dm->ntile_y == 0) {
            tile_setup_2d(ncpu, config->nx, config->ny, &dm->ntile_x, &dm->ntile_y); // setup the number of tiles in x and y directions
            if (icpu == 0) fprintf(stderr, "xy: nx=%d ny=%d ibnd_x=%d ibnd_y=%d ntile_x=%d ntile_y=%d whalo_x=%d whalo_y=%d\n", config->nx, config->ny, config->ibnd_x, config->ibnd_y, dm->ntile_x, dm->ntile_y, config->whalo_x, config->whalo_y);
        }

        dm->nz = config->nz;
        dm->mz = config->nz + 2 * config->whalo_z;
        dm->kts = config->whalo_z;
        dm->kte = config->whalo_z + config->nz - 1;
        dm->kms = dm->kts - config->whalo_z;
        dm->kme = dm->kte + config->whalo_z;
        dm->d_neighbor = MPI_PROC_NULL;
        dm->u_neighbor = MPI_PROC_NULL;

        domain_partition_2d(ncpu, icpu, config->nx, config->ny, config->ibnd_x, config->ibnd_y, dm->ntile_x, dm->ntile_y, config->whalo_x, config->whalo_y, &dm->its, &dm->ite, &dm->jts, &dm->jte, &dm->w_neighbor, &dm->e_neighbor, &dm->s_neighbor, &dm->n_neighbor);
        dm->ims = dm->its - config->whalo_x;
        dm->ime = dm->ite + config->whalo_x;
        dm->jms = dm->jts - config->whalo_y;
        dm->jme = dm->jte + config->whalo_y;
        dm->nx = dm->ite - dm->its + 1;
        dm->ny = dm->jte - dm->jts + 1;
        dm->mx = dm->ime - dm->ims + 1;
        dm->my = dm->jme - dm->jms + 1;
    }else if (dir == 'Y') {
        // no partition in y direction
        dm->ntile_y = 1;
        if (dm->ntile_x == 0 || dm->ntile_z == 0) {
            tile_setup_2d(ncpu, config->nx, config->nz, &dm->ntile_x, &dm->ntile_z); // setup the number of tiles in x and z directions
            if (icpu == 0) fprintf(stderr, "xz: nx=%d nz=%d ibnd_x=%d ibnd_z=%d ntile_x=%d ntile_z=%d whalo_x=%d whalo_z=%d\n", config->nx, config->nz, config->ibnd_x, config->ibnd_z, dm->ntile_x, dm->ntile_z, config->whalo_x, config->whalo_z);
        }

        dm->ny = config->ny;
        dm->my = config->ny + 2 * config->whalo_y;
        dm->jts = config->whalo_y;
        dm->jte = config->whalo_y + config->ny - 1;
        dm->jms = dm->jts - config->whalo_y;
        dm->jme = dm->jte + config->whalo_y;
        dm->s_neighbor = MPI_PROC_NULL;
        dm->n_neighbor = MPI_PROC_NULL;

        domain_partition_2d(ncpu, icpu, config->nx, config->nz, config->ibnd_x, config->ibnd_z, dm->ntile_x, dm->ntile_z, config->whalo_x, config->whalo_z, &dm->its, &dm->ite, &dm->kts, &dm->kte, &dm->w_neighbor, &dm->e_neighbor, &dm->d_neighbor, &dm->u_neighbor);
        dm->ims = dm->its - config->whalo_x;
        dm->ime = dm->ite + config->whalo_x;
        dm->kms = dm->kts - config->whalo_z;
        dm->kme = dm->kte + config->whalo_z;
        dm->nx = dm->ite - dm->its + 1;
        dm->nz = dm->kte - dm->kts + 1;
        dm->mx = dm->ime - dm->ims + 1;
        dm->mz = dm->kme - dm->kms + 1;
    }else if (dir == 'X') {
        // no partition in x direction
        dm->ntile_x = 1;
        if (dm->ntile_y == 0 || dm->ntile_z == 0) {
            tile_setup_2d(ncpu, config->ny, config->nz, &dm->ntile_y, &dm->ntile_z); // setup the number of tiles in y and z directions
            if (icpu == 0) fprintf(stderr, "yz: ny=%d nz=%d ibnd_y=%d ibnd_z=%d ntile_y=%d ntile_z=%d whalo_y=%d whalo_z=%d\n", config->ny, config->nz, config->ibnd_y, config->ibnd_z, dm->ntile_y, dm->ntile_z, config->whalo_y, config->whalo_z);
        }

        dm->nx = config->nx;
        dm->mx = config->nx + 2 * config->whalo_x;
        dm->its = config->whalo_x;
        dm->ite = config->whalo_x + config->nx - 1;
        dm->ims = dm->its - config->whalo_x;
        dm->ime = dm->ite + config->whalo_x;
        dm->w_neighbor = MPI_PROC_NULL;
        dm->e_neighbor = MPI_PROC_NULL;

        domain_partition_2d(ncpu, icpu, config->ny, config->nz, config->ibnd_y, config->ibnd_z, dm->ntile_y, dm->ntile_z, config->whalo_y, config->whalo_z, &dm->jts, &dm->jte, &dm->kts, &dm->kte, &dm->s_neighbor, &dm->n_neighbor, &dm->d_neighbor, &dm->u_neighbor);
        dm->jms = dm->jts - config->whalo_y;
        dm->jme = dm->jte + config->whalo_y;
        dm->kms = dm->kts - config->whalo_z;
        dm->kme = dm->kte + config->whalo_z;
        dm->ny = dm->jte - dm->jts + 1;
        dm->nz = dm->kte - dm->kts + 1;
        dm->my = dm->jme - dm->jms + 1;
        dm->mz = dm->kme - dm->kms + 1;
    }
    dm->nslice = dm->nx * dm->ny;
    dm->ncube = dm->nslice * dm->nz;
    dm->mslice = dm->mx * dm->my;
    dm->mcube = dm->mslice * dm->mz;

#if defined(_DEBUG_)
    if (dir == 'Z') {
        fprintf(stderr, "%c%05d: i=(%05d  %05d), j=(%05d  %05d), k=(%03d  %03d), neighbor=(%05d  %05d  %05d  %05d)\n", dir, icpu, dm->its, dm->ite, dm->jts, dm->jte, dm->kts, dm->kte, dm->w_neighbor, dm->e_neighbor, dm->s_neighbor, dm->n_neighbor);
    }else if (dir == 'Y') {
        fprintf(stderr, "%c%05d: i=(%05d  %05d), j=(%05d  %05d), k=(%03d  %03d), neighbor=(%05d  %05d  %05d  %05d)\n", dir, icpu, dm->its, dm->ite, dm->jts, dm->jte, dm->kts, dm->kte, dm->w_neighbor, dm->e_neighbor, dm->d_neighbor, dm->u_neighbor);
    }else if (dir == 'X') {
        fprintf(stderr, "%c%05d: i=(%05d  %05d), j=(%05d  %05d), k=(%03d  %03d), neighbor=(%05d  %05d  %05d  %05d)\n", dir, icpu, dm->its, dm->ite, dm->jts, dm->jte, dm->kts, dm->kte, dm->s_neighbor, dm->n_neighbor, dm->d_neighbor, dm->u_neighbor);
    }
#endif
}

double partition_transpose_index(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max, int radix) {
    struct t_domain *dm_src, *dm_dst;
    int ierr, icpu, ncpu, isrc, idst;
    int i, j, k, m, n;
    int nsend, nrecv;
    int is, ie, js, je, ks, ke;
    double t;

    memset(stat, 0, sizeof(struct t_stat));
    stat->total_time =  MPI_Wtime();

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
#if defined(_COMM_)
    int COMM_OFFSET = 0;
#else
    int COMM_OFFSET = (icpu / da_src->ntile_x) * da_src->ntile_x;
    ncpu = da_src->ntile_x;
    icpu -= COMM_OFFSET;
#endif

    radix = MIN(ncpu, MAX(2, radix));

    double *s_bufr, *r_bufr;
    s_bufr = bufr;
    r_bufr = PTR_ADD(bufr,  ncpu * max);

#if !defined(_NOMEM_)
    /********************************************************************************************************
     * shift forward and pack
     ********************************************************************************************************/
    isrc = icpu;
    dm_src = da_src + isrc;
    for (n = 0; n < ncpu; n++) {
        idst = (icpu + n) % ncpu;
        dm_dst = da_dst + idst;
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
        // pack
        t = MPI_Wtime();
        stat->nbyte_pack += pack(fld_src, s_bufr + max * n, is-dm_src->ims, ie-dm_src->ims, js-dm_src->jms, je-dm_src->jms, ks-dm_src->kms, ke-dm_src->kms, dm_src->mx, dm_src->mslice);
        stat->pack_time += MPI_Wtime() - t;
    }
#endif

    /********************************************************************************************************
     * transpose by modified bruck algorithm, we exchange all the remainder messages to avoid packing
     ********************************************************************************************************/
    for (k = 1; k < ncpu; k *= radix) {
        for (j = 1; j < radix; j++) {
            n = j*k;
            if (n >= ncpu) break;

#if defined(_PACKING_)
            for (nsend = 0, i = n; i < ncpu; i += k*radix) {
                m = MIN(k, ncpu - i) * max;
#if !defined(_NOMEM_)
                memcpy(r_bufr + nsend, s_bufr + i * max, m * sizeof(double));
#endif
                nsend += m;
            }
            nrecv = nsend;
#else
            nsend = nrecv = max*(ncpu-n);
#endif
            ASSERT((m = nsend * sizeof(double)) >= 0, "bytes sent are too large, please increase number of precesses.");
            ASSERT((m = nrecv * sizeof(double)) >= 0, "bytes received are too large, please increase number of precesses.");

            isrc = COMM_OFFSET + (icpu - n + ncpu) % ncpu;
            idst = COMM_OFFSET + (icpu + n) % ncpu;
            t = MPI_Wtime();
#if defined(_PACKING_)
            ierr = MPI_Sendrecv(r_bufr,                 nsend, MPI_DOUBLE, idst, itag,
                                PTR_ADD(r_bufr, nsend), nrecv, MPI_DOUBLE, isrc, itag,
                                comm, MPI_STATUS_IGNORE);
#else
            ierr = MPI_Sendrecv(PTR_ADD(s_bufr, max * n), nsend, MPI_DOUBLE, idst, itag,
                                PTR_ADD(r_bufr, max * n), nrecv, MPI_DOUBLE, isrc, itag,
                                comm, MPI_STATUS_IGNORE);
#endif
            stat->mpi_time += MPI_Wtime() - t;
            stat->nbyte_send += nsend * sizeof(double);
            stat->nbyte_recv += nrecv * sizeof(double);

#if !defined(_NOMEM_)
            for (nrecv = nsend, i = n; i < ncpu; i += k*radix) {
                m = MIN(k, ncpu - i) * max;
#if defined(_PACKING_)
                memcpy(s_bufr + i * max, r_bufr + nrecv,   m * sizeof(double));
                nrecv += m;
#else
                memcpy(s_bufr + i * max, r_bufr + i * max, m * sizeof(double));
#endif
            }
#endif
        }
    }

    /********************************************************************************************************
     * shift backward unpack
     ********************************************************************************************************/
#if !defined(_NOMEM_)
    idst = icpu;
    dm_dst = da_dst + idst;
    for (n = 0; n < ncpu; n++) {
        isrc = ((ncpu - 1 - n) + icpu + 1) % ncpu;
        dm_src = da_src + isrc;
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
        // unpack
        t = MPI_Wtime();
        stat->nbyte_unpack += unpack(s_bufr + max * n, fld_dst, is-dm_dst->ims, ie-dm_dst->ims, js-dm_dst->jms, je-dm_dst->jms, ks-dm_dst->kms, ke-dm_dst->kms, dm_dst->mx, dm_dst->mslice);
        stat->unpack_time += MPI_Wtime() - t;
    }
#endif

    // return total time
    stat->total_time = MPI_Wtime() - stat->total_time;
    return stat->total_time;
}

double partition_transpose_ringall(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max) {
    struct t_domain *dm_src, *dm_dst;
    int ierr, icpu, ncpu, isrc, idst;
    int K, L, m, n;
    int is, ie, js, je, ks, ke;
    double t;

    memset(stat, 0, sizeof(struct t_stat));
    stat->total_time =  MPI_Wtime();

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
#if defined(_COMM_)
    int COMM_OFFSET = 0;
#else
    int COMM_OFFSET = (icpu / da_src->ntile_x) * da_src->ntile_x;
    ncpu = da_src->ntile_x;
    icpu -= COMM_OFFSET;
#endif
    K = ncpu * max;
    ASSERT((m = K * sizeof(double)) >= 0, "bytes received are too large, please increase number of precesses.");

    MPI_Status status;
    double *s_bufr, *r_bufr;
    s_bufr = bufr;
    r_bufr = PTR_ADD(bufr, K);

    /********************************************************************************************************
     * shift forward and pack
     ********************************************************************************************************/
    isrc = icpu;
    dm_src = da_src + isrc;
    for (L = n = 0; n < ncpu; n++) {
        idst = (icpu + n) % ncpu;
        dm_dst = da_dst + idst;
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
#if !defined(_NOMEM_)
        // pack
        t = MPI_Wtime();
        stat->nbyte_pack += pack(fld_src, s_bufr + L, is-dm_src->ims, ie-dm_src->ims, js-dm_src->jms, je-dm_src->jms, ks-dm_src->kms, ke-dm_src->kms, dm_src->mx, dm_src->mslice);
        stat->pack_time += MPI_Wtime() - t;
#endif
        // next
        L += (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
    }

    /********************************************************************************************************
     * unpack itself
     ********************************************************************************************************/
    isrc = icpu;
    dm_src = da_src + isrc;
    idst = icpu;
    dm_dst = da_dst + idst;
    // coordinates of overlaped region
    is = MAX(dm_src->its, dm_dst->its);
    ie = MIN(dm_src->ite, dm_dst->ite);
    js = MAX(dm_src->jts, dm_dst->jts);
    je = MIN(dm_src->jte, dm_dst->jte);
    ks = MAX(dm_src->kts, dm_dst->kts);
    ke = MIN(dm_src->kte, dm_dst->kte);
#if !defined(_NOMEM_)
    // unpack
    t = MPI_Wtime();
    stat->nbyte_unpack += unpack(s_bufr, fld_dst, is-dm_dst->ims, ie-dm_dst->ims, js-dm_dst->jms, je-dm_dst->jms, ks-dm_dst->kms, ke-dm_dst->kms, dm_dst->mx, dm_dst->mslice);
    stat->unpack_time += MPI_Wtime() - t;
#endif
    // point to the buffer for remainder processes
    L -= m = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
    s_bufr = PTR_ADD(s_bufr, m);

    int inext, iprev;
    inext = COMM_OFFSET + (icpu + 1) % ncpu;
    iprev = COMM_OFFSET + (icpu - 1 + ncpu) % ncpu;

    for (n = 1; n < ncpu; n++) {
        /********************************************************************************************************
         * send and receive
         ********************************************************************************************************/
        t = MPI_Wtime();
        ierr = MPI_Sendrecv(s_bufr, L, MPI_DOUBLE, inext, itag,
                            r_bufr, K, MPI_DOUBLE, iprev, itag,
                            comm, &status);
        stat->mpi_time += MPI_Wtime() - t;
        stat->nbyte_send += L * sizeof(double);
        MPI_Get_count(&status, MPI_DOUBLE, &L);
        stat->nbyte_recv += L * sizeof(double);

        /********************************************************************************************************
         * unpack
         ********************************************************************************************************/
        isrc = (icpu - n + ncpu) % ncpu;
        idst = icpu;
        dm_src = da_src + isrc;
        dm_dst = da_dst + idst;
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
#if !defined(_NOMEM_)
        // unpack
        t = MPI_Wtime();
        stat->nbyte_unpack += unpack(r_bufr, fld_dst, is-dm_dst->ims, ie-dm_dst->ims, js-dm_dst->jms, je-dm_dst->jms, ks-dm_dst->kms, ke-dm_dst->kms, dm_dst->mx, dm_dst->mslice);
        stat->unpack_time += MPI_Wtime() - t;
#endif
        // next
        L -= m = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
        // swap the send and receive buffer to avoid memcpy
        bufr = s_bufr;
        s_bufr = PTR_ADD(r_bufr, m); // point to the buffer for remainder processes
        r_bufr = bufr;
    }

    // total time
    stat->total_time = MPI_Wtime() - stat->total_time;
    return stat->total_time;
}

double partition_transpose_ringone(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max, int radix) {
    struct t_domain *dm_src, *dm_dst;
    int ierr, icpu, ncpu, isrc, idst;
    int K, L, m, n;
    int is, ie, js, je, ks, ke;
    double t;

    memset(stat, 0, sizeof(struct t_stat));
    stat->total_time =  MPI_Wtime();

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
#if defined(_COMM_)
    int COMM_OFFSET = 0;
#else
    int COMM_OFFSET = (icpu / da_src->ntile_x) * da_src->ntile_x;
    ncpu = da_src->ntile_x;
    icpu -= COMM_OFFSET;
#endif
    ASSERT((m = max  * sizeof(double)) >= 0, "bytes received are too large, please increase number of processes.");

    radix = MIN(ncpu-1, MAX(1, radix));

    double *s_bufr, *r_bufr;
    s_bufr = bufr;
    r_bufr = PTR_ADD(bufr, radix * max);

#if !defined(_NOMEM_)
    /********************************************************************************************************
     * local copy
     ********************************************************************************************************/
    isrc =  icpu;
    dm_src = da_src + isrc;
    idst =  icpu;
    dm_dst = da_dst + idst;
    // coordinates of overlaped region
    is = MAX(dm_src->its, dm_dst->its);
    ie = MIN(dm_src->ite, dm_dst->ite);
    js = MAX(dm_src->jts, dm_dst->jts);
    je = MIN(dm_src->jte, dm_dst->jte);
    ks = MAX(dm_src->kts, dm_dst->kts);
    ke = MIN(dm_src->kte, dm_dst->kte);
    // pack
    t = MPI_Wtime();
    stat->nbyte_pack += pack(fld_src, s_bufr, is-dm_src->ims, ie-dm_src->ims, js-dm_src->jms, je-dm_src->jms, ks-dm_src->kms, ke-dm_src->kms, dm_src->mx, dm_src->mslice);
    stat->pack_time += MPI_Wtime() - t;
    // unpack
    t = MPI_Wtime();
    stat->nbyte_unpack += unpack(s_bufr, fld_dst, is-dm_dst->ims, ie-dm_dst->ims, js-dm_dst->jms, je-dm_dst->jms, ks-dm_dst->kms, ke-dm_dst->kms, dm_dst->mx, dm_dst->mslice);
    stat->unpack_time += MPI_Wtime() - t;
#endif

    MPI_Request *request = (MPI_Request *)malloc(2*radix*sizeof(MPI_Request));
    MPI_Request *PTR_REQ;
    MPI_Status status;

   for (PTR_REQ = request, L = 0, n = 1; n < ncpu; n++) {
        /********************************************************************************************************
         * pack, send and receive
         ********************************************************************************************************/
        isrc = icpu;
        idst = (icpu + n) % ncpu;
        dm_src = da_src + isrc;
        dm_dst = da_dst + idst;
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
#if !defined(_NOMEM_)
        // pack
        t = MPI_Wtime();
        stat->nbyte_pack += pack(fld_src, PTR_ADD(s_bufr, L*max), is-dm_src->ims, ie-dm_src->ims, js-dm_src->jms, je-dm_src->jms, ks-dm_src->kms, ke-dm_src->kms, dm_src->mx, dm_src->mslice);
        stat->pack_time += MPI_Wtime() - t;
#endif
        // send and receive
        isrc = COMM_OFFSET + (icpu - n + ncpu) % ncpu;
        idst = COMM_OFFSET + (icpu + n) % ncpu;
        m = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
        t = MPI_Wtime();
        ierr = MPI_Irecv(PTR_ADD(r_bufr, L*max), max, MPI_DOUBLE, isrc, itag, comm, PTR_REQ++);
        ierr = MPI_Isend(PTR_ADD(s_bufr, L*max), m,   MPI_DOUBLE, idst, itag, comm, PTR_REQ++);
        stat->mpi_time += MPI_Wtime() - t;
        stat->nbyte_send += m * sizeof(double);

        // continue until we reach the end of current round (MPI_Irecv and MPI_Isend)
        L++;
        if (L < radix && n + 1 < ncpu) continue;


        /********************************************************************************************************
         * wait and unpack
         ********************************************************************************************************/
        do {
            t = MPI_Wtime();
            ierr = MPI_Waitany(PTR_REQ - request, request, &K, &status);
            stat->mpi_time += MPI_Wtime() - t;
            // if finished
            if (K != MPI_UNDEFINED && K % 2 == 0) {
                MPI_Get_count(&status, MPI_DOUBLE, &m);
                stat->nbyte_recv += m * sizeof(double);

#if !defined(_NOMEM_)
                isrc = status.MPI_SOURCE - COMM_OFFSET;
                idst = icpu;
                dm_src = da_src + isrc;
                dm_dst = da_dst + idst;
                // coordinates of overlaped region
                is = MAX(dm_src->its, dm_dst->its);
                ie = MIN(dm_src->ite, dm_dst->ite);
                js = MAX(dm_src->jts, dm_dst->jts);
                je = MIN(dm_src->jte, dm_dst->jte);
                ks = MAX(dm_src->kts, dm_dst->kts);
                ke = MIN(dm_src->kte, dm_dst->kte);
                // unpack
                t = MPI_Wtime();
                stat->nbyte_unpack += unpack(r_bufr + K/2*max, fld_dst, is-dm_dst->ims, ie-dm_dst->ims, js-dm_dst->jms, je-dm_dst->jms, ks-dm_dst->kms, ke-dm_dst->kms, dm_dst->mx, dm_dst->mslice);
                stat->unpack_time += MPI_Wtime() - t;
#endif
            }
        } while (K != MPI_UNDEFINED);

        // next round
        PTR_REQ = request;
        L = 0;
    }

    free(request);

    // return total time
    stat->total_time = MPI_Wtime() - stat->total_time;
    return stat->total_time;
}

double partition_transpose(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max, int nsend, int nrecv) {
    struct t_domain *dm_src, *dm_dst;
    int ierr, icpu, ncpu, isrc, idst;
    int K, L, m, n;
    int is, ie, js, je, ks, ke;
    double t;

    memset(stat, 0, sizeof(struct t_stat));
    stat->total_time = MPI_Wtime();

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
#if defined(_COMM_)
    int COMM_OFFSET = 0;
#else
    int COMM_OFFSET = (icpu / da_src->ntile_x) * da_src->ntile_x;
    ncpu = da_src->ntile_x;
    icpu -= COMM_OFFSET;
#endif
    ASSERT((m = max  * sizeof(double)) >= 0, "bytes received are too large, please increase number of precesses.");

    double *s_bufr, *r_bufr;
    s_bufr = bufr;
    r_bufr = PTR_ADD(bufr, ncpu * max);

#if !defined(_NOMEM_)
    /********************************************************************************************************
     * local copy
     ********************************************************************************************************/
    isrc =  icpu;
    dm_src = da_src + isrc;
    idst =  icpu;
    dm_dst = da_dst + idst;
    // if overlaped
    if ( dm_src->nslice > 0 && dm_dst->nslice > 0 &&
       !(dm_src->ite < dm_dst->its || dm_src->its > dm_dst->ite ||
         dm_src->jte < dm_dst->jts || dm_src->jts > dm_dst->jte ||
         dm_src->kte < dm_dst->kts || dm_src->kts > dm_dst->kte) ) {
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
        // pack
        t = MPI_Wtime();
        stat->nbyte_pack += pack(fld_src, s_bufr, is-dm_src->ims, ie-dm_src->ims, js-dm_src->jms, je-dm_src->jms, ks-dm_src->kms, ke-dm_src->kms, dm_src->mx, dm_src->mslice);
        stat->pack_time += MPI_Wtime() - t;
        // unpack
        t = MPI_Wtime();
        stat->nbyte_unpack += unpack(s_bufr, fld_dst, is-dm_dst->ims, ie-dm_dst->ims, js-dm_dst->jms, je-dm_dst->jms, ks-dm_dst->kms, ke-dm_dst->kms, dm_dst->mx, dm_dst->mslice);
        stat->unpack_time += MPI_Wtime() - t;
    }
#endif

    nrecv--, nsend--;
    MPI_Request *request = (MPI_Request *)malloc((nrecv + nsend) * sizeof(MPI_Request));
    MPI_Status status;


    for (K = L = 0, n = 1; n < ncpu; n++) {
        /********************************************************************************************************
         * receive
         ********************************************************************************************************/
        isrc = (icpu - n + ncpu) % ncpu;
        dm_src = da_src + isrc;
        idst = icpu;
        dm_dst = da_dst + idst;
        // if overlaped
        if ( dm_src->nslice > 0 && dm_dst->nslice > 0 &&
           !(dm_src->ite < dm_dst->its || dm_src->its > dm_dst->ite ||
             dm_src->jte < dm_dst->jts || dm_src->jts > dm_dst->jte ||
             dm_src->kte < dm_dst->kts || dm_src->kts > dm_dst->kte) ) {
            // receive
            t = MPI_Wtime();
            ierr = MPI_Irecv(PTR_ADD(r_bufr, K*max), max, MPI_DOUBLE, COMM_OFFSET + isrc, itag, comm, &request[K]);
            stat->mpi_time += MPI_Wtime() - t;
            K++;
        }

        /********************************************************************************************************
         * send
         ********************************************************************************************************/
        isrc =  icpu;
        dm_src = da_src + isrc;
        idst =  (icpu + n) % ncpu;
        dm_dst = da_dst + idst;
        // if overlaped
        if ( dm_src->nslice > 0 && dm_dst->nslice > 0 &&
           !(dm_src->ite < dm_dst->its || dm_src->its > dm_dst->ite ||
             dm_src->jte < dm_dst->jts || dm_src->jts > dm_dst->jte ||
             dm_src->kte < dm_dst->kts || dm_src->kts > dm_dst->kte) ) {
            // coordinates of overlaped region
            is = MAX(dm_src->its, dm_dst->its);
            ie = MIN(dm_src->ite, dm_dst->ite);
            js = MAX(dm_src->jts, dm_dst->jts);
            je = MIN(dm_src->jte, dm_dst->jte);
            ks = MAX(dm_src->kts, dm_dst->kts);
            ke = MIN(dm_src->kte, dm_dst->kte);
#if !defined(_NOMEM_)
            // pack
            t = MPI_Wtime();
            stat->nbyte_pack += pack(fld_src, s_bufr + L*max, is-dm_src->ims, ie-dm_src->ims, js-dm_src->jms, je-dm_src->jms, ks-dm_src->kms, ke-dm_src->kms, dm_src->mx, dm_src->mslice);
            stat->pack_time += MPI_Wtime() - t;
#endif
            // send
            m = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
            t = MPI_Wtime();
            ierr = MPI_Isend(PTR_ADD(s_bufr, L*max), m, MPI_DOUBLE, COMM_OFFSET + idst, itag, comm, &request[nrecv + L]);
            stat->mpi_time += MPI_Wtime() - t;
            stat->nbyte_send += m * sizeof(double);
            L++;
        }
    }
    ASSERT(K == nrecv && L == nsend, "Something impossible did happened");

    /********************************************************************************************************
     * wait finished
     ********************************************************************************************************/
    do {
        t = MPI_Wtime();
        ierr = MPI_Waitany(nrecv, request, &K, &status);
        stat->mpi_time += MPI_Wtime() - t;
        // if finished
        if (K != MPI_UNDEFINED) {
            MPI_Get_count(&status, MPI_DOUBLE, &m);
            stat->nbyte_recv += m * sizeof(double);

#if !defined(_NOMEM_)
            isrc = status.MPI_SOURCE - COMM_OFFSET;
            dm_src = da_src + isrc;
            dm_dst = da_dst + icpu;
            // coordinates of overlaped region
            is = MAX(dm_src->its, dm_dst->its);
            ie = MIN(dm_src->ite, dm_dst->ite);
            js = MAX(dm_src->jts, dm_dst->jts);
            je = MIN(dm_src->jte, dm_dst->jte);
            ks = MAX(dm_src->kts, dm_dst->kts);
            ke = MIN(dm_src->kte, dm_dst->kte);
            // unpack
            t = MPI_Wtime();
            stat->nbyte_unpack += unpack(r_bufr + K*max, fld_dst, is-dm_dst->ims, ie-dm_dst->ims, js-dm_dst->jms, je-dm_dst->jms, ks-dm_dst->kms, ke-dm_dst->kms, dm_dst->mx, dm_dst->mslice);
            stat->unpack_time += MPI_Wtime() - t;
#endif
        }
    } while (K != MPI_UNDEFINED);

    t = MPI_Wtime();
    ierr = MPI_Waitall(nsend, &request[nrecv], MPI_STATUSES_IGNORE);
    stat->mpi_time += MPI_Wtime() - t;

    free(request);

    // return total time
    stat->total_time = MPI_Wtime() - stat->total_time;
    return stat->total_time;
}

double partition_transpose_alltoallv(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int *count_offset) {
    struct t_domain *dm_src, *dm_dst;
    int ierr, icpu, ncpu, isrc, idst;
    int nsend, nrecv, m;
    int is, ie, js, je, ks, ke;
    double t;

    memset(stat, 0, sizeof(struct t_stat));
    stat->total_time = MPI_Wtime();

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);

    double *s_bufr, *r_bufr;
    s_bufr = bufr;
    r_bufr = PTR_ADD(bufr, da_src[icpu].ncube);

    int *s_count, *s_offset, *r_count, *r_offset;
#if defined(_NOPAD_)
    s_count  = count_offset;
    s_offset = NULL;
    r_count  = count_offset + ncpu;
    r_offset = NULL;
#else
    s_count  = count_offset;
    s_offset = count_offset + ncpu;
    r_count  = count_offset + ncpu * 2;
    r_offset = count_offset + ncpu * 3;
#endif

    /********************************************************************************************************
     * pack and send
     ********************************************************************************************************/
    dm_src = da_src + icpu;
    for (nsend = idst = 0, dm_dst = da_dst; idst < ncpu; idst++, dm_dst++) {
#if !defined(_NOPAD_)
        s_offset[idst] = nsend;
#endif
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
#if !defined(_NOMEM_)
        // pack
        t = MPI_Wtime();
        stat->nbyte_pack += pack(fld_src, s_bufr + nsend, is-dm_src->ims, ie-dm_src->ims, js-dm_src->jms, je-dm_src->jms, ks-dm_src->kms, ke-dm_src->kms, dm_src->mx, dm_src->mslice);
        stat->pack_time += MPI_Wtime() - t;
#endif
        // offset and count
        s_count[idst] = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
        nsend += s_count[idst];

        ASSERT((m = s_count[idst]  * sizeof(double)) >= 0, "bytes sent are too large, please increase number of precesses.");
    }

    /********************************************************************************************************
     * receive
     ********************************************************************************************************/
    dm_dst = da_dst + icpu;
    for (nrecv = isrc = 0, dm_src = da_src; isrc < ncpu; isrc++, dm_src++) {
#if !defined(_NOPAD_)
        r_offset[isrc] = nrecv;
#endif
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
        // offset and count
        r_count[isrc] = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
        nrecv += r_count[isrc];

        ASSERT((m = r_count[isrc]  * sizeof(double)) >= 0, "bytes received are too large, please increase number of precesses.");
    }

    t = MPI_Wtime();
    ierr = MPI_Alltoallv(s_bufr, s_count, s_offset, MPI_DOUBLE, r_bufr, r_count, r_offset, MPI_DOUBLE, comm);
    stat->mpi_time = MPI_Wtime() - t;
    stat->nbyte_send += nsend * sizeof(double);
    stat->nbyte_recv += nrecv * sizeof(double);

    /********************************************************************************************************
     * unpack
     ********************************************************************************************************/
#if !defined(_NOMEM_)
    dm_dst = da_dst + icpu;
    for (nrecv = isrc = 0, dm_src = da_src; isrc < ncpu; isrc++, dm_src++) {
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
        // unpack
        t = MPI_Wtime();
        stat->nbyte_unpack += unpack(r_bufr + nrecv, fld_dst, is-dm_dst->ims, ie-dm_dst->ims, js-dm_dst->jms, je-dm_dst->jms, ks-dm_dst->kms, ke-dm_dst->kms, dm_dst->mx, dm_dst->mslice);
        stat->unpack_time += MPI_Wtime() - t;
        // count
        nrecv += (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
    }
#endif

    // return total time
    stat->total_time = MPI_Wtime() - stat->total_time;
    return stat->total_time;
}

double partition_transpose_alltoall(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max) {
    struct t_domain *dm_src, *dm_dst;
    int ierr, icpu, ncpu, isrc, idst;
    int m;
    int is, ie, js, je, ks, ke;
    double t;

    memset(stat, 0, sizeof(struct t_stat));
    stat->total_time = MPI_Wtime();

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
    ASSERT((m = max  * sizeof(double)) >= 0, "bytes sent or received are too large, please increase number of precesses.");

    double *s_bufr, *r_bufr;
    s_bufr = bufr;
    r_bufr = PTR_ADD(bufr, ncpu * max);

    /********************************************************************************************************
     * pack and send
     ********************************************************************************************************/
#if !defined(_NOMEM_)
    dm_src = da_src + icpu;
    for (idst = 0, dm_dst = da_dst; idst < ncpu; idst++, dm_dst++) {
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
        // pack
        t = MPI_Wtime();
        stat->nbyte_pack += pack(fld_src, s_bufr + idst * max, is-dm_src->ims, ie-dm_src->ims, js-dm_src->jms, je-dm_src->jms, ks-dm_src->kms, ke-dm_src->kms, dm_src->mx, dm_src->mslice);
        stat->pack_time += MPI_Wtime() - t;
    }
#endif

    t = MPI_Wtime();
    ierr = MPI_Alltoall(s_bufr, max, MPI_DOUBLE, r_bufr, max, MPI_DOUBLE, comm);
    stat->mpi_time = MPI_Wtime() - t;
    stat->nbyte_send += ncpu * max * sizeof(double);
    stat->nbyte_recv += ncpu * max * sizeof(double);

    /********************************************************************************************************
     * unpack
     ********************************************************************************************************/
#if !defined(_NOMEM_)
    dm_dst = da_dst + icpu;
    for (isrc = 0, dm_src = da_src; isrc < ncpu; isrc++, dm_src++) {
        // coordinates of overlaped region
        is = MAX(dm_src->its, dm_dst->its);
        ie = MIN(dm_src->ite, dm_dst->ite);
        js = MAX(dm_src->jts, dm_dst->jts);
        je = MIN(dm_src->jte, dm_dst->jte);
        ks = MAX(dm_src->kts, dm_dst->kts);
        ke = MIN(dm_src->kte, dm_dst->kte);
        // unpack
        t = MPI_Wtime();
        stat->nbyte_unpack += unpack(r_bufr + isrc * max, fld_dst, is-dm_dst->ims, ie-dm_dst->ims, js-dm_dst->jms, je-dm_dst->jms, ks-dm_dst->kms, ke-dm_dst->kms, dm_dst->mx, dm_dst->mslice);
        stat->unpack_time += MPI_Wtime() - t;
    }
#endif

    // return total time
    stat->total_time = MPI_Wtime() - stat->total_time;
    return stat->total_time;
}

void partition_transpose_workmem_ring(MPI_Comm comm, int iflg, struct t_domain *da_src, struct t_domain *da_dst, double **bufr, int *max) {
    struct t_domain *dm_src, *dm_dst;
    int icpu, ncpu, isrc, idst;
    int m;
    int is, ie, js, je, ks, ke;

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
#if defined(_COMM_)
    int COMM_OFFSET = 0;
#else
    int COMM_OFFSET = (icpu / da_src->ntile_x) * da_src->ntile_x;
    ncpu = da_src->ntile_x;
    icpu -= COMM_OFFSET;
#endif

    int nx = 0, ny = 0;
    for (isrc = 0; isrc < ncpu; isrc++) {
        dm_src = da_src + isrc;
        if (nx < dm_src->nx) nx = dm_src->nx;
    }
    for (idst = 0; idst < ncpu; idst++) {
        dm_dst = da_dst + idst;
        if (ny < dm_dst->ny) ny = dm_dst->ny;
    }

    *max = nx * ny * dm_src->nz;
    m = (*max) * ncpu * 2;

    if (iflg) {
#if defined(_NOMEM_)
        *bufr = NULL;
#else
        //*bufr = (double *)malloc(m * sizeof(double));
        shmc((void **)bufr, m * sizeof(double), 999);
#endif
#if !defined(_NOMEM_)
        if (*bufr == NULL) {
            fprintf(stderr, "Program terminated: can't allocate %u memory for bufr!\n", m);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
#if defined(_DEBUG_)
        fprintf(stderr, "%d: (%d, %d) allocated %u memory %p\n", icpu, dm_src->ntile_x, dm_src->ntile_y, m, *bufr);
#endif
#endif
    }else{
#if defined(_NOMEM_)
        *bufr = NULL;
#else
        //free(*bufr);
        shmd((void **)bufr, m * sizeof(double), 999);
#endif
#if defined(_DEBUG_)
        fprintf(stderr, "%d: (%d, %d) deallocated memory\n", icpu, dm_src->ntile_x, dm_src->ntile_y);
#endif
    }
}

void partition_transpose_workmem(MPI_Comm comm, int iflg, struct t_domain *da_src, struct t_domain *da_dst, double **bufr, int *max, int *nsend, int *nrecv) {
    struct t_domain *dm_src, *dm_dst;
    int icpu, ncpu, isrc, idst;
    int m;
    int is, ie, js, je, ks, ke;

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
#if defined(_COMM_)
    int COMM_OFFSET = 0;
#else
    int COMM_OFFSET = (icpu / da_src->ntile_x) * da_src->ntile_x;
    ncpu = da_src->ntile_x;
    icpu -= COMM_OFFSET;
#endif

    *max = 0;
    // receive the overlaped region from other processors
    *nrecv = 0;
    dm_dst = da_dst + icpu;
    for (isrc = 0; isrc < ncpu; isrc++) {
        dm_src = da_src + isrc;
        // if overlaped
        if ( dm_src->nslice > 0 && dm_dst->nslice > 0 &&
           !(dm_src->ite < dm_dst->its || dm_src->its > dm_dst->ite ||
             dm_src->jte < dm_dst->jts || dm_src->jts > dm_dst->jte ||
             dm_src->kte < dm_dst->kts || dm_src->kts > dm_dst->kte) ) {
            // coordinates of overlaped region
            is = MAX(dm_src->its, dm_dst->its);
            ie = MIN(dm_src->ite, dm_dst->ite);
            js = MAX(dm_src->jts, dm_dst->jts);
            je = MIN(dm_src->jte, dm_dst->jte);
            ks = MAX(dm_src->kts, dm_dst->kts);
            ke = MIN(dm_src->kte, dm_dst->kte);
            // size of overlaped region
            m = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
            if (*max < m) *max = m;
            // count how many cpus to which I need to receive
            (*nrecv)++;
        }
    }

    // send the overlaped region to other processors
    *nsend = 0;
    dm_src = da_src + icpu;
    for (idst = 0; idst < ncpu; idst++) {
        dm_dst = da_dst + idst;
        // if overlaped
        if ( dm_src->nslice > 0 && dm_dst->nslice > 0 &&
           !(dm_src->ite < dm_dst->its || dm_src->its > dm_dst->ite ||
             dm_src->jte < dm_dst->jts || dm_src->jts > dm_dst->jte ||
             dm_src->kte < dm_dst->kts || dm_src->kts > dm_dst->kte) ) {
            // coordinates of overlaped region
            is = MAX(dm_src->its, dm_dst->its);
            ie = MIN(dm_src->ite, dm_dst->ite);
            js = MAX(dm_src->jts, dm_dst->jts);
            je = MIN(dm_src->jte, dm_dst->jte);
            ks = MAX(dm_src->kts, dm_dst->kts);
            ke = MIN(dm_src->kte, dm_dst->kte);
            // size of overlaped region
            m = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
            if (*max < m) *max = m;
            (*nsend)++;
        }
    }

    m = (*max) * ((*nrecv) + (*nsend));

    if (iflg) {
#if defined(_NOMEM_)
        *bufr = NULL;
#else
        //*bufr = (double *)malloc(m * sizeof(double));
        shmc((void **)bufr, m * sizeof(double), 999);
#endif
#if !defined(_NOMEM_)
        if (*bufr == NULL) {
            fprintf(stderr, "Program terminated: can't allocate %u memory for bufr!\n", m);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
#if defined(_DEBUG_)
        fprintf(stderr, "%d: (%d, %d) allocated %u memory %p\n", icpu, dm_src->ntile_x, dm_src->ntile_y, m, *bufr);
#endif
#endif
    }else{
#if defined(_NOMEM_)
        *bufr = NULL;
#else
        //free(*bufr);
        shmd((void **)bufr, m * sizeof(double), 999);
#endif
#if defined(_DEBUG_)
        fprintf(stderr, "%d: (%d, %d) deallocated memory\n", icpu, dm_src->ntile_x, dm_src->ntile_y);
#endif
    }
}

void partition_transpose_workmem_alltoallv(MPI_Comm comm, int iflg, struct t_domain *dm_src, struct t_domain *dm_dst, double **bufr, int **count_offset) {
    int icpu, ncpu;
    int m, n;

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);

#if defined(_NOPAD_)
    n = ncpu * 2; // send_counts, recv_counts
#else
    n = ncpu * 4; // send_displs, send_counts, recv_displs, recv_counts
#endif
    ASSERT(ncpu * n * sizeof(int) < 135291469824, "MPI_Alltoallv needs more than 128Gb memory for send_counts(ncpu), send_displs(ncpu), recv_counts(ncpu), and recv_displs(ncpu)");
    m = dm_src->ncube + dm_dst->ncube;

    if (iflg) {
#if defined(_NOMEM_)
        *bufr = NULL;
#else
        //*bufr = (double *)malloc(m * sizeof(double));
        shmc((void **)bufr, m * sizeof(double), 999);
#endif
#if !defined(_NOMEM_)
        if (*bufr == NULL) {
            fprintf(stderr, "Program terminated: can't allocate %u memory for bufr!\n", m);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
#if defined(_DEBUG_)
        fprintf(stderr, "%d: (%d, %d) allocated %u memory %p\n", icpu, dm_src->ntile_x, dm_src->ntile_y, m, *bufr);
#endif
#endif
        *count_offset = (int *)malloc(n * sizeof(int));
        if (*count_offset == NULL) {
            fprintf(stderr, "Program terminated: can't allocate %u memory for count_offset!\n", n);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
#if defined(_DEBUG_)
        fprintf(stderr, "%d: (%d, %d) allocated %u memory %p\n", icpu, dm_src->ntile_x, dm_src->ntile_y, n, *count_offset);
#endif
    }else{
#if defined(_NOMEM_)
        *bufr = NULL;
#else
        //free(*bufr);
        shmd((void **)bufr, m * sizeof(double), 999);
#endif
        free(*count_offset);
#if defined(_DEBUG_)
        fprintf(stderr, "%d: (%d, %d) deallocated memory\n", icpu, dm_src->ntile_x, dm_src->ntile_y);
#endif
    }
}
