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

int send_overlap(char dir, char pass, struct t_domain *dm, struct t_domain *dn, int *is, int *ie, int *js, int *je, int *ks, int *ke) {
  int iret = 0;
  if (dir == 'x') {
    *js = dm->jts;
    *je = dm->jte;
    *ks = dm->kts;
    *ke = dm->kte;

    if (pass == '-') { // west neighbor
      if (dn->icpu >= dm->icpu) {
        *is = dn->ite + 1 - dn->nx_g;
        *ie = dn->ite + dn->whalo_x - dn->nx_g;
      }else{
        *is = dn->ite + 1;
        *ie = dn->ite + dn->whalo_x;
      }
      if ( dm->nslice > 0 && dn->nslice > 0 &&
         !(dm->ite < *is || dm->its > *ie) ) {
        // coordinates of overlaped region
        *is = MAX(dm->its, *is);
        *ie = MIN(dm->ite, *ie);

        iret = 1;
      }
    }else if (pass == '+') { // east neighbor
      if (dn->icpu <= dm->icpu) {
        *is = dn->its - dn->whalo_x + dn->nx_g;
        *ie = dn->its - 1 + dn->nx_g;
      }else{
        *is = dn->its - dn->whalo_x;
        *ie = dn->its - 1;
      }
      if ( dm->nslice > 0 && dn->nslice > 0 &&
         !(dm->ite < *is || dm->its > *ie) ) {
        // coordinates of overlaped region
        *is = MAX(dm->its, *is);
        *ie = MIN(dm->ite, *ie);

        iret = 1;
      }
    }
  }else if (dir == 'y') {
    *is = dm->ims;
    *ie = dm->ime;
    *ks = dm->kts;
    *ke = dm->kte;

    if (pass == '-') { // south neighbor
      if (dn->icpu >= dm->icpu) {
        *js = dn->jte + 1 - dn->ny_g;
        *je = dn->jte + dn->whalo_y - dn->ny_g;
      }else{
        *js = dn->jte + 1;
        *je = dn->jte + dn->whalo_y;
      }
      if ( dm->nslice > 0 && dn->nslice > 0 &&
         !(dm->jte < *js || dm->jts > *je) ) {
        // coordinates of overlaped region
        *js = MAX(dm->jts, *js);
        *je = MIN(dm->jte, *je);

        iret = 1;
      }
    }else if (pass == '+') { // north neighbor
      if (dn->icpu <= dm->icpu) {
        *js = dn->jts - dn->whalo_y + dn->ny_g;
        *je = dn->jts - 1 + dn->ny_g;
      }else{
        *js = dn->jts - dn->whalo_y;
        *je = dn->jts - 1;
      }
      if ( dm->nslice > 0 && dn->nslice > 0 &&
         !(dm->jte < *js || dm->jts > *je) ) {
        // coordinates of overlaped region
        *js = MAX(dm->jts, *js);
        *je = MIN(dm->jte, *je);

        iret = 1;
      }
    }
  }
  return iret;
}

inline int recv_overlap(char dir, char pass, struct t_domain *dm, struct t_domain *dn, int *is, int *ie, int *js, int *je, int *ks, int *ke) {
  int iret = 0;

  if (dir == 'x') {
    *js = dm->jts;
    *je = dm->jte;
    *ks = dm->kts;
    *ke = dm->kte;

    if (pass == '-') { // west neighbor
      if (dn->icpu >= dm->icpu) {
        *is = dm->its - dm->whalo_x + dm->nx_g;
        *ie = dm->its - 1 + dm->nx_g;
      }else{
        *is = dm->its - dm->whalo_x;
        *ie = dm->its - 1;
      }
      if ( dm->nslice > 0 && dn->nslice > 0 &&
         !(*ie < dn->its || *is > dn->ite) ) {
        // coordinates of overlaped region
        *is = MAX(*is, dn->its);
        *ie = MIN(*ie, dn->ite);
        if (dn->icpu >= dm->icpu) {
          *is -= dm->nx_g;
          *ie -= dm->nx_g;
        }

        iret = 1;
      }
    }else if (pass == '+') { // east neighbor
      if (dn->icpu <= dm->icpu) {
        *is = dm->ite + 1 - dm->nx_g;
        *ie = dm->ite + dm->whalo_x - dm->nx_g;
      }else{
        *is = dm->ite + 1;
        *ie = dm->ite + dm->whalo_x;
      }
      if ( dm->nslice > 0 && dn->nslice > 0 &&
         !(*ie < dn->its || *is > dn->ite) ) {
        // coordinates of overlaped region
        *is = MAX(*is, dn->its);
        *ie = MIN(*ie, dn->ite);
        if (dn->icpu <= dm->icpu) {
          *is += dm->nx_g;
          *ie += dm->nx_g;
        }

        iret = 1;
      }
    }
  }else if (dir == 'y') {
    *is = dm->ims;
    *ie = dm->ime;
    *ks = dm->kts;
    *ke = dm->kte;

    if (pass == '-') { // south neighbor
      if (dn->icpu >= dm->icpu) {
        *js = dm->jts - dm->whalo_y + dm->ny_g;
        *je = dm->jts - 1 + dm->ny_g;
      }else{
        *js = dm->jts - dm->whalo_y;
        *je = dm->jts - 1;
      }
      if ( dm->nslice > 0 && dn->nslice > 0 &&
         !(*je < dn->jts || *js > dn->jte) ) {
        // coordinates of overlaped region
        *js = MAX(*js, dn->jts);
        *je = MIN(*je, dn->jte);
        if (dn->icpu >= dm->icpu) {
          *js -= dm->ny_g;
          *je -= dm->ny_g;
        }

        iret = 1;
      }
    }else if (pass == '+') { // north neighbor
      if (dn->icpu <= dm->icpu) {
        *js = dm->jte + 1 - dm->ny_g;
        *je = dm->jte + dm->whalo_y - dm->ny_g;
      }else{
        *js = dm->jte + 1;
        *je = dm->jte + dm->whalo_y;
      }
      if ( dm->nslice > 0 && dn->nslice > 0 &&
         !(*je < dn->jts || *js > dn->jte) ) {
        // coordinates of overlaped region
        *js = MAX(*js, dn->jts);
        *je = MIN(*je, dn->jte);
        if (dn->icpu <= dm->icpu) {
          *js += dm->ny_g;
          *je += dm->ny_g;
        }

        iret = 1;
      }
    }
  }

  return iret;
}

double halo_exchange(int icpu, struct t_domain *da, double *fld, struct t_stat *stat, double *bufr) {
  double t;
  struct t_domain *dm, *dn;
  int i_prev, i_next;
  int m_prev, m_next, m_recv;
  int i, j, k, m, n;
  int is, ie, js, je, ks, ke;
  int ierr;

  dm = da + icpu;
  m = MAX(dm->whalo_x * dm->ny * dm->nz, dm->mx * dm->whalo_y * dm->nz);
  if (m <= 0) return 0;
#if defined(_NOMEM_)
  double *r_bufr0 = NULL, *r_bufrm = NULL, *s_bufr0 = NULL, *s_bufrm = NULL;
#else
  double *r_bufr0 = bufr, *r_bufrm = bufr + m, *s_bufr0 = bufr + 2 * m, *s_bufrm = bufr + 3 * m;
#endif

  memset((void *)stat, 0, sizeof(struct t_stat));
  stat->total_time = MPI_Wtime();

  MPI_Status status[4];
  MPI_Request request[4];

for (char dir='x'; dir<='y'; dir++) {
  if (dir == 'x') {
    i_prev = dm->w_neighbor;
    i_next = dm->e_neighbor;
    n = dm->ntile_x;
  } else if (dir == 'y') {
    i_prev = dm->s_neighbor;
    i_next = dm->n_neighbor;
    n = dm->ntile_y;
  }

  while (--n) {
    m_prev = m_next = 0;

    /***************************************************************************
     * co-location and pack
     **************************************************************************/
    // region sent to negative neighbour
    dn = da + i_prev;
    if ( send_overlap(dir, '-', dm, dn, &is, &ie, &js, &je, &ks, &ke) ) {
      // size of overlaped region
      m_prev = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
#if !defined(_NOMEM_)
      // pack
      t = MPI_Wtime();
      stat->nbyte_pack += pack(fld, s_bufr0, is-dm->ims, ie-dm->ims, js-dm->jms, je-dm->jms, ks-dm->kms, ke-dm->kms, dm->mx, dm->mslice);
      stat->pack_time += MPI_Wtime() - t;
#endif
      //fprintf(stderr, "s%c%d#%d<-%d: %3d:%3d, %3d:%3d, %3d:%3d\n", dir, n, i_prev, icpu, is, ie, js, je, ks, ke);
    }

    // region sent to positive neighbour
    dn = da + i_next;
    if ( send_overlap(dir, '+', dm, dn, &is, &ie, &js, &je, &ks, &ke) ) {
      // size of overlaped region
      m_next = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
#if !defined(_NOMEM_)
      // pack
      t = MPI_Wtime();
      stat->nbyte_pack += pack(fld, s_bufrm, is-dm->ims, ie-dm->ims, js-dm->jms, je-dm->jms, ks-dm->kms, ke-dm->kms, dm->mx, dm->mslice);
      stat->pack_time += MPI_Wtime() - t;
#endif
      //fprintf(stderr, "s%c%d#%d->%d: %3d:%3d, %3d:%3d, %3d:%3d\n", dir, n, icpu, i_next, is, ie, js, je, ks, ke);
    }

    // have exchange with all the overlaped processes
    if (m_prev == 0 && m_next == 0) break;

    /***************************************************************************
     * exchange halo region
     **************************************************************************/
#if defined(_SENDRECV_)
    t = MPI_Wtime();

    // positive direction
    ierr = MPI_Sendrecv(s_bufrm, m_next, MPI_DOUBLE, i_next, 1,
                        r_bufr0, m,      MPI_DOUBLE, i_prev, 1,
                        MPI_COMM_WORLD, &status[0]);
    // negative direction
    ierr = MPI_Sendrecv(s_bufr0, m_prev, MPI_DOUBLE, i_prev, 2,
                        r_bufrm, m,      MPI_DOUBLE, i_next, 2,
                        MPI_COMM_WORLD, &status[1]);

    stat->mpi_time += MPI_Wtime() - t;
#elif defined(_ISENDIRECV_)
    t = MPI_Wtime();

    // positive direction
    ierr = MPI_Isend(s_bufrm, m_next, MPI_DOUBLE, i_next, 1, MPI_COMM_WORLD, &request[2]);
    ierr = MPI_Irecv(r_bufr0, m,      MPI_DOUBLE, i_prev, 1, MPI_COMM_WORLD, &request[0]);
    // negative direction
    ierr = MPI_Isend(s_bufr0, m_prev, MPI_DOUBLE, i_prev, 2, MPI_COMM_WORLD, &request[3]);
    ierr = MPI_Irecv(r_bufrm, m,      MPI_DOUBLE, i_next, 2, MPI_COMM_WORLD, &request[1]);

    ierr = MPI_Waitall(4, request, status);

    stat->mpi_time += MPI_Wtime() - t;
#elif defined(_ISENDRECV_)
    t = MPI_Wtime();

    // positive direction
    ierr = MPI_Isend(s_bufrm, m_next, MPI_DOUBLE, i_next, 1, MPI_COMM_WORLD, &request[2]);
    ierr = MPI_Recv( r_bufr0, m,      MPI_DOUBLE, i_prev, 1, MPI_COMM_WORLD, &status[0]);
    // negative direction
    ierr = MPI_Isend(s_bufr0, m_prev, MPI_DOUBLE, i_prev, 2, MPI_COMM_WORLD, &request[3]);
    ierr = MPI_Recv( r_bufrm, m,      MPI_DOUBLE, i_next, 2, MPI_COMM_WORLD, &status[1]);

    ierr = MPI_Waitall(2, request, MPI_STATUSES_IGNORE);

    stat->mpi_time += MPI_Wtime() - t;
#elif defined(_SENDIRECV_)
    t = MPI_Wtime();

    // positive direction
    ierr = MPI_Irecv(r_bufr0, m,      MPI_DOUBLE, i_prev, 1, MPI_COMM_WORLD, &request[0]);
    ierr = MPI_Send( s_bufrm, m_next, MPI_DOUBLE, i_next, 1, MPI_COMM_WORLD);
    // negative direction
    ierr = MPI_Irecv(r_bufrm, m,      MPI_DOUBLE, i_next, 2, MPI_COMM_WORLD, &request[1]);
    ierr = MPI_Send( s_bufr0, m_prev, MPI_DOUBLE, i_prev, 2, MPI_COMM_WORLD);

    ierr = MPI_Waitall(2, request, status);

    stat->mpi_time += MPI_Wtime() - t;
#else
    k = dir == 'x' ? dm->icpu % dm->ntile_x : dm->icpu / dm->ntile_x;

    t = MPI_Wtime();

    // prevent from deadlock
    if (k % 2 == 0 ) {
      // positive direction
      ierr = MPI_Send(s_bufrm, m_next, MPI_DOUBLE, i_next, 1, MPI_COMM_WORLD);
      ierr = MPI_Recv(r_bufr0, m,      MPI_DOUBLE, i_prev, 1, MPI_COMM_WORLD, &status[0]);
      // negative direction
      ierr = MPI_Send(s_bufr0, m_prev, MPI_DOUBLE, i_prev, 2, MPI_COMM_WORLD);
      ierr = MPI_Recv(r_bufrm, m,      MPI_DOUBLE, i_next, 2, MPI_COMM_WORLD, &status[1]);
    }else{
      // positive direction
      ierr = MPI_Recv(r_bufr0, m,      MPI_DOUBLE, i_prev, 1, MPI_COMM_WORLD, &status[0]);
      ierr = MPI_Send(s_bufrm, m_next, MPI_DOUBLE, i_next, 1, MPI_COMM_WORLD);
      // negative direction
      ierr = MPI_Recv(r_bufrm, m,      MPI_DOUBLE, i_next, 2, MPI_COMM_WORLD, &status[1]);
      ierr = MPI_Send(s_bufr0, m_prev, MPI_DOUBLE, i_prev, 2, MPI_COMM_WORLD);
    }

    stat->mpi_time += MPI_Wtime() - t;
#endif

  stat->nbyte_send += (m_prev + m_next) * sizeof(double);
  MPI_Get_count(&status[0], MPI_DOUBLE, &m_recv);
  stat->nbyte_recv += m_recv * sizeof(double);
  MPI_Get_count(&status[1], MPI_DOUBLE, &m_recv);
  stat->nbyte_recv += m_recv * sizeof(double);

    /***************************************************************************
     * co-location and unpack
     **************************************************************************/
    // region receive from negative neighour
    dn = da + i_prev;
    if ( recv_overlap(dir, '-', dm, dn, &is, &ie, &js, &je, &ks, &ke) ) {
#if !defined(_NOMEM_)
      // unpack
      t = MPI_Wtime();
      stat->nbyte_unpack += unpack(r_bufr0, fld, is-dm->ims, ie-dm->ims, js-dm->jms, je-dm->jms, ks-dm->kms, ke-dm->kms, dm->mx, dm->mslice);
      stat->unpack_time += MPI_Wtime() - t;
#endif
      //fprintf(stderr, "r%c%d#%d->%d: %3d:%3d, %3d:%3d, %3d:%3d\n", dir, n, i_prev, icpu, is, ie, js, je, ks, ke);
    }

    // region receive from positive neighbour
    dn = da + i_next;
    if ( recv_overlap(dir, '+', dm, dn, &is, &ie, &js, &je, &ks, &ke) ) {
#if !defined(_NOMEM_)
      // unpack
      t = MPI_Wtime();
      stat->nbyte_unpack += unpack(r_bufrm, fld, is-dm->ims, ie-dm->ims, js-dm->jms, je-dm->jms, ks-dm->kms, ke-dm->kms, dm->mx, dm->mslice);
      stat->unpack_time += MPI_Wtime() - t;
#endif
      //fprintf(stderr, "r%c%d#%d<-%d: %3d:%3d, %3d:%3d, %3d:%3d\n", dir, n, icpu, i_next, is, ie, js, je, ks, ke);
    }

    /***************************************************************************
     * neighbor's neighbor ......
     **************************************************************************/
    if (dir == 'x') {
      dn = da + i_prev;
      i_prev = dn->w_neighbor;
      dn = da + i_next;
      i_next = dn->e_neighbor;
    } else if (dir == 'y') {
      dn = da + i_prev;
      i_prev = dn->s_neighbor;
      dn = da + i_next;
      i_next = dn->n_neighbor;
    }
  }
}
  // return total time
  stat->total_time = MPI_Wtime() - stat->total_time;
  return stat->total_time;
}

void halo_exchange_workmem(int iflg, struct t_domain *dm, double **bufr) {
  int m = MAX(dm->whalo_x * dm->ny * dm->nz, dm->mx * dm->whalo_y * dm->nz);
  if (m <= 0) return;

  if (iflg) {
    // the first 2*m is for recv, the last 2*m is for send, so total 4*m
#if defined(_NOMEM_)
    *bufr = NULL;
#else
    *bufr = (double *)malloc(4 * m * sizeof(double));
#endif
#if !defined(_NOMEM_)
    if (*bufr == NULL) {
      fprintf(stderr, "Program terminated: can't allocate memory for bufr!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
#if defined(_DEBUG_)
    fprintf(stderr, "%d: allocated %d bytes memory %p\n", dm->icpu, m, *bufr);
#endif
#endif
  }else{
#if defined(_NOMEM_)
    *bufr = NULL;
#else
    free(*bufr);
#endif
#if !defined(_NOMEM_) && defined(_DEBUG_)
    fprintf(stderr, "%d: deallocated memory\n", dm->icpu);
#endif
  }
}
