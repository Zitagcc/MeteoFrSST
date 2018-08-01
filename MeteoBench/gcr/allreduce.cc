/********************************************************
 * Created by Yongjun ZHENG on 01 Feb 2017.             *
 * Email: yongjun.zheng@meteo.fr                        *
 * Copyright © 2017 Yongjun ZHENG. All rights reserved. *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <inttypes.h>

#include <mpi.h>

#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define PTR_ADD(a,b) ( (a)==NULL ? NULL : (a)+(b) )
#define PTR_ADDC(a,b) ( (a)==NULL ? NULL : (char *)(a)+(b) )

#define sstmac_app_name allreduce

void mxor(int icpu, int mask, int radix, int *partner) {
    int digit, nextm;
    digit = ((icpu / mask) % radix) * mask;
    nextm = mask *  radix;
    icpu -= digit;
    for (int k = 1; k < radix; k++) {
        digit += mask;
        *partner = icpu + digit % nextm;
        partner++;
    }
}

struct t_stat {
    double mpi_time;
    uint64_t nbyte_send;
    uint64_t nbyte_recv;
};

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

void do_stat(struct t_stat *stat_avg, struct t_stat *stat_max, struct t_stat *stat_min, struct t_stat *stat, int navg) {
    // avg
    stat_avg->mpi_time += stat->mpi_time;
    stat_avg->nbyte_send += stat->nbyte_send;
    stat_avg->nbyte_recv += stat->nbyte_recv;
    if (navg) {
        stat_avg->mpi_time /= navg;
        stat_avg->nbyte_send /= navg;
        stat_avg->nbyte_recv /= navg;
    }
    // max
    if (stat->mpi_time > stat_max->mpi_time) stat_max->mpi_time = stat->mpi_time;
    if (stat->nbyte_send > stat_max->nbyte_send) stat_max->nbyte_send = stat->nbyte_send;
    if (stat->nbyte_recv > stat_max->nbyte_recv) stat_max->nbyte_recv = stat->nbyte_recv;
    // min
    if (stat->mpi_time < stat_min->mpi_time) stat_min->mpi_time = stat->mpi_time;
    if (stat->nbyte_send < stat_min->nbyte_send) stat_min->nbyte_send = stat->nbyte_send;
    if (stat->nbyte_recv < stat_min->nbyte_recv) stat_min->nbyte_recv = stat->nbyte_recv;
}

void output_stat(int icpu, int ncpu, struct t_stat stat, char *fn) {
    int ierr;

    if (icpu == 0) {
        struct t_stat stat_avg, stat_max, stat_min;
        stat_avg = stat;
        stat_max = stat;
        stat_min = stat;

        int k;
        for (k = 1; k < ncpu; k++) {
            ierr = MPI_Recv(&stat, sizeof(struct t_stat), MPI_BYTE, MPI_ANY_SOURCE, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            do_stat(&stat_avg, &stat_max, &stat_min, &stat, k+1==ncpu ? ncpu : 0);
        }

#define PRIu64 "llu"
        FILE *fp = fopen(fn, "wt");
        fprintf(fp, "%15.10lf %10" PRIu64 " %10" PRIu64 " ",  stat_avg.mpi_time, stat_avg.nbyte_send, stat_avg.nbyte_recv);
        fprintf(fp, "%15.10lf %10" PRIu64 " %10" PRIu64 " ",  stat_max.mpi_time, stat_max.nbyte_send, stat_max.nbyte_recv);
        fprintf(fp, "%15.10lf %10" PRIu64 " %10" PRIu64 "\n", stat_min.mpi_time, stat_min.nbyte_send, stat_min.nbyte_recv);
        fclose(fp);
    }else{
        ierr = MPI_Send(&stat, sizeof(struct t_stat), MPI_BYTE, 0, 999, MPI_COMM_WORLD);
    }
}


double gcr_allreduce(MPI_Comm comm, int itag, struct t_stat *stat, int max, double *sbufr, double *rbufr) {
    int ierr, icpu, ncpu, isrc, idst;
    double t;

    memset(stat, 0, sizeof(struct t_stat));

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);

    t = MPI_Wtime();
    ierr = MPI_Allreduce(sbufr, rbufr, max, MPI_DOUBLE, MPI_SUM, comm);
    stat->mpi_time = MPI_Wtime() - t;
    stat->nbyte_send = max * ncpu * sizeof(double);
    stat->nbyte_recv = max * ncpu * sizeof(double);

    return stat->mpi_time;
}

// first reduce then broadcast
double gcr_allreduce_rb(MPI_Comm comm, int itag, struct t_stat *stat, int max, double *sbufr, double *rbufr) {
    int ierr, icpu, ncpu, isrc, idst;
    double t;

    memset(stat, 0, sizeof(struct t_stat));

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);

    t = MPI_Wtime();
    ierr = MPI_Reduce(sbufr,rbufr, max, MPI_DOUBLE, MPI_SUM, 0, comm);
    ierr = MPI_Bcast(rbufr, max, MPI_DOUBLE, 0, comm);
    stat->mpi_time = MPI_Wtime() - t;
    stat->nbyte_send = max * ncpu * sizeof(double);
    stat->nbyte_recv = max * ncpu * sizeof(double);

    return stat->mpi_time;
}

// shift ring-k
double gcr_allreduce_rs(MPI_Comm comm, int itag, struct t_stat *stat, int max, double *sbufr, double *rbufr, int radix) {
    int ierr, icpu, ncpu, isrc, idst;
    int K, L, m, n;
    double t;

    memset(stat, 0, sizeof(struct t_stat));

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);

    radix = MIN(ncpu-1, MAX(1, radix));
#if defined(_NOMEM_)
    double *tmp = NULL;
#else
    double *tmp = (double *)malloc(radix*max*sizeof(double));
    memset(rbufr, 0, max*sizeof(double));
#endif

    MPI_Request *request = (MPI_Request *)malloc(2*radix*sizeof(MPI_Request));
    MPI_Request *PTR_REQ;
    MPI_Status status;

    for (PTR_REQ = request, L = 0, n = 1; n < ncpu; n++) {
        /********************************************************************************************************
         * send and receive
         ********************************************************************************************************/
        isrc = (icpu - n + ncpu) % ncpu;
        idst = (icpu + n) % ncpu;
        t = MPI_Wtime();
        ierr = MPI_Irecv(PTR_ADD(tmp, max*L), max, MPI_DOUBLE, isrc, itag, comm, PTR_REQ++);
        ierr = MPI_Isend(sbufr, max, MPI_DOUBLE, idst, itag, comm, PTR_REQ++);
        stat->mpi_time += MPI_Wtime() - t;
        stat->nbyte_send += max * sizeof(double);

        // continue until we reach the end of current round (MPI_Irecv and MPI_Isend)
        L++;
        if (L < radix && n + 1 < ncpu) continue;


        /********************************************************************************************************
         * wait
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
                int j = max * (K/2);
                for (int i = 0; i < max; i++, j++) {
                    rbufr[i] += tmp[j]; // MPI_SUM
                }
#endif
            }
        } while (K != MPI_UNDEFINED);

        // next round,
        PTR_REQ = request;
        L = 0;
    }

    free(request);
#if !defined(_NOMEM_)
    free(tmp);
#endif

    return stat->mpi_time;
}

// fixed ring-k
double gcr_allreduce_rf(MPI_Comm comm, int itag, struct t_stat *stat, int max, double *sbufr, double *rbufr, int radix) {
    int ierr, icpu, ncpu, isrc, idst;
    int K, L, m, n;
    double t;

    memset(stat, 0, sizeof(struct t_stat));

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);

    radix = MIN(ncpu-1, MAX(1, radix));
#if defined(_NOMEM_)
    double *tmp = NULL;
#else
    double *tmp = (double *)malloc((radix+1)*max*sizeof(double));
    memset(rbufr, 0, max*sizeof(double));
    memcpy(PTR_ADD(tmp, max*radix), sbufr, max*sizeof(double));
#endif

    MPI_Request *request = (MPI_Request *)malloc(2*radix*sizeof(MPI_Request));
    MPI_Request *PTR_REQ;
    MPI_Status status;

    for (PTR_REQ = request, L = 0, n = 1; n < ncpu; n++) {
        /********************************************************************************************************
         * send and receive
         ********************************************************************************************************/
        isrc = (icpu - (L + 1) + ncpu) % ncpu;
        idst = (icpu + (L + 1)) % ncpu;
        t = MPI_Wtime();
        ierr = MPI_Irecv(PTR_ADD(tmp, max*L),     max, MPI_DOUBLE, isrc, itag, comm, PTR_REQ++);
        ierr = MPI_Isend(PTR_ADD(tmp, max*radix), max, MPI_DOUBLE, idst, itag, comm, PTR_REQ++);
        stat->mpi_time += MPI_Wtime() - t;
        stat->nbyte_send += max * sizeof(double);

        // continue until we reach the end of current round (MPI_Irecv and MPI_Isend)
        L++;
        if (L < radix && n + 1 < ncpu) continue;

        /********************************************************************************************************
         * wait
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
                int j = max * (K/2);
                for (int i = 0; i < max; i++, j++) {
                  rbufr[i] += tmp[j]; // MPI_SUM
                }
#endif
            }
        } while (K != MPI_UNDEFINED);

        // next round
        PTR_REQ = request;
        L = 0;
#if !defined(_NOMEM_)
        memcpy(PTR_ADD(tmp, max*radix), PTR_ADD(tmp, max*(radix-1)), max*sizeof(double));
#endif
    }

    free(request);
#if !defined(_NOMEM_)
    free(tmp);
#endif

    return stat->mpi_time;
}

// recursive radix
double gcr_allreduce_rk(MPI_Comm comm, int itag, struct t_stat *stat, int max, double *sbufr, double *rbufr, int radix) {
    int ierr, icpu, ncpu, isrc, idst;
    int i, j, K, m, mask;
    double t;

    memset(stat, 0, sizeof(struct t_stat));

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);

    int por, rem;
    for (i = 0; i < radix; i++) {
        if (radix+i <= 16) {
            por = pow(radix+i, (int)(log(ncpu)/log(radix+i)));
            rem = ncpu - por;
            if (rem <= por) {
                radix += i;
                break;
            }
        }
        
        if (radix-i > 1) {
            por = pow(radix-i, (int)(log(ncpu)/log(radix-i)));
            rem = ncpu - por;
            if (rem <= por) {
                radix -= i;
                break;
            }
        }
    }
    if (icpu == 0) fprintf(stderr, "radix=%d\n", radix);

    int m_radix = radix - 1;
    double *tmp = NULL;
#if !defined(_NOMEM_)
    tmp = (double *)malloc(m_radix*max*sizeof(double));
    memcpy(rbufr, sbufr, max*sizeof(double));
#endif
    int *partner = (int *)malloc(m_radix*sizeof(int));
    MPI_Request *request = (MPI_Request *)malloc(2*m_radix*sizeof(MPI_Request));
    MPI_Request *PTR_REQ;
    MPI_Status status;

    if (icpu < rem) {
        isrc = por + icpu;
        t = MPI_Wtime();
        ierr = MPI_Recv(tmp, max, MPI_DOUBLE, isrc, itag, comm, MPI_STATUS_IGNORE);
        stat->mpi_time += MPI_Wtime() - t;
        stat->nbyte_recv += max * sizeof(double);
#if !defined(_NOMEM_)
        for (i = 0; i < max; i++) {
            rbufr[i] += tmp[i];
        }
#endif
    }else if (icpu >= por) {
        idst = icpu - por;
        t = MPI_Wtime();
        ierr = MPI_Send(sbufr, max, MPI_DOUBLE, idst, itag, comm);
        stat->mpi_time += MPI_Wtime() - t;
        stat->nbyte_send += max * sizeof(double);
    }

    if (icpu < por) {
        for (mask = 1; mask < por; mask *= radix) {
            mxor(icpu, mask, radix, partner);
            
            /********************************************************************************************************
             * pose send and recv
             ********************************************************************************************************/
            for (PTR_REQ = request, i = 0; i < m_radix; i++) {
                isrc = idst = partner[i];
                t = MPI_Wtime();
                ierr = MPI_Irecv(PTR_ADD(tmp, max*i), max, MPI_DOUBLE, isrc, itag, comm, PTR_REQ++);
                ierr = MPI_Isend(rbufr,               max, MPI_DOUBLE, idst, itag, comm, PTR_REQ++);
                stat->mpi_time += MPI_Wtime() - t;
                stat->nbyte_send += max * sizeof(double);
            }

            /********************************************************************************************************
             * wait
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
                    j = (K / 2) * max;
                    for (i = 0; i < max; i++, j++) {
                        rbufr[i] += tmp[j];
                    }
#endif
                }
            } while (K != MPI_UNDEFINED);
        }
    }

    if (icpu < rem) {
        idst = por + icpu;
        t = MPI_Wtime();
        ierr = MPI_Send(rbufr, max, MPI_DOUBLE, idst, itag, comm);
        stat->mpi_time += MPI_Wtime() - t;
        stat->nbyte_send += max * sizeof(double);
    }else if (icpu >= por) {
        isrc = icpu - por;
        t = MPI_Wtime();
        ierr = MPI_Recv(rbufr, max, MPI_DOUBLE, isrc, itag, comm, MPI_STATUS_IGNORE);
        stat->mpi_time += MPI_Wtime() - t;
        stat->nbyte_recv += max * sizeof(double);
    }

    free(partner);
    free(request);
#if !defined(_NOMEM_)
    free(tmp);
#endif

    return stat->mpi_time;
}

int main(int argc, char *argv[]) {
    int icpu, ncpu;
    int mth, max, radix;
    double t;

    MPI_Comm comm = MPI_COMM_WORLD;
    int itag = 111;
    struct t_stat stat;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);

    if (argc < 3) {
        if (icpu==0) fprintf(stderr, "Usage: %s method size radix\n", argv[0]);
        MPI_Finalize();
        return 0;
    }
    mth = atoi(argv[1]);
    max = atoi(argv[2]);
    if (mth == 5) {
        radix = 4;
        if (argc == 4) radix = atoi(argv[3]);
        radix = MAX(2, radix);
        radix = MIN(ncpu, radix);
    }else{
        radix = 0;
    }

#if defined(_NOMEM_)
    double *sbufr = NULL;
    double *rbufr = NULL;
#else
    double *sbufr = (double *)malloc(max*sizeof(double));
    double *rbufr = (double *)malloc(max*sizeof(double));
    for (int i = 0; i < max; i++) {
        sbufr[i] = icpu;
    }
#endif

    switch (mth) {
        case 1:
            t = gcr_allreduce(comm, itag, &stat, max, sbufr, rbufr);
            break;
        case 2:
            t = gcr_allreduce_rb(comm, itag, &stat, max, sbufr, rbufr);
            break;
        case 3:
            t = gcr_allreduce_rs(comm, itag, &stat, max, sbufr, rbufr, radix);
            break;
        case 4:
            t = gcr_allreduce_rf(comm, itag, &stat, max, sbufr, rbufr, radix);
            break;
        case 5:
            t = gcr_allreduce_rk(comm, itag, &stat, max, sbufr, rbufr, radix);
            break;
    }

    char fn[128];
    if (icpu == 0) {
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
        sprintf(fn, "%08d/%d.%08d.%04d", ncpu, mth, max, radix);
    }
    output_stat(icpu, ncpu, stat, fn);

#if !defined(_NOMEM_)
    if (icpu == 0) {
        fprintf(stderr, "%d: m=%d, t=%lf", icpu, mth, t);
        for (int i = 0; i < max; i++) {
            fprintf(stderr, ", %lf", rbufr[i]);
        }
        fprintf(stderr, "\n");
    }

    free(sbufr);
    free(rbufr);
#endif

    MPI_Finalize();
    return 0;
}
