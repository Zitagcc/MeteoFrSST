#ifndef __MPP__
#define __MPP__

#include <mpi.h>
#include <inttypes.h>
#include <iostream>

#define INDEX(i, j, k, ncol, nslice) ( (k)*(nslice) + (j)*(ncol) + (i) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define PTR_ADD(a,b) ( (a)==NULL ? NULL : (a)+(b) )

#define _ASSERT_ 1
#if defined(_ASSERT_)
#   define ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

struct t_config {
    double xs, xe, ys, ye, zs, ze, dx, dy, dz;
    int nx, ny, nz;
    int ibnd_x, ibnd_y, ibnd_z;
    int ntile_x, ntile_y, ntile_z;
    int whalo_x, whalo_y, whalo_z;
    int niter;
    int radix;
    int ntrunc;
};

struct t_domain {
    int ncpu, icpu;
    int nx_g, ny_g, nz_g;
    int ntile_x, ntile_y, ntile_z;
    int whalo_x, whalo_y, whalo_z;
    int nx, ny, nz, mx, my, mz;
    int its, ite, jts, jte, kts, kte;
    int ims, ime, jms, jme, kms, kme;
    int nslice, ncube, mslice, mcube;
    int w_neighbor, e_neighbor, s_neighbor, n_neighbor, d_neighbor, u_neighbor;
    char dir;
};

struct t_stat {
    double cal_time;
    double mpi_time;
    double pack_time;
    double unpack_time;
    double total_time;
    uint64_t nbyte_pack;
    uint64_t nbyte_unpack;
    uint64_t nbyte_send;
    uint64_t nbyte_recv;
};

double mysecond();
int pack(double *src, double *dst, int is, int ie, int js, int je, int ks, int ke, int ncol, int nslice);
int unpack(double *src, double *dst, int is, int ie, int js, int je, int ks, int ke, int ncol, int nslice);
//int prime_factors(const int num, const int m, int *factors);
//void tile_setup(const int ncpu, const int nx, const int ny, int *ntile_x, int *ntile_y);
void domain_partition(char dir, int icpu, int ncpu, struct t_config *config, struct t_domain *dm);
void partition_transpose_workmem_alltoallv(MPI_Comm comm, int iflg, struct t_domain *dm_src, struct t_domain *dm_dst, double **bufr, int **count_offset);
void partition_transpose_workmem(MPI_Comm comm, int iflg, struct t_domain *da_src, struct t_domain *da_dst, double **bufr, int *max, int *nsend, int *nrecv);
void partition_transpose_workmem_ring(MPI_Comm comm, int iflg, struct t_domain *da_src, struct t_domain *da_dst, double **bufr, int *max);
double partition_transpose_alltoallv(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int *count_offset);
double partition_transpose(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max, int nsend, int nrecv);
double partition_transpose_ringall(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max);
double partition_transpose_ringone(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max, int radix);
double partition_transpose_index(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max, int radix);
double partition_transpose_alltoall(MPI_Comm comm, int itag, struct t_domain *da_src, struct t_domain *da_dst, double *fld_src, double *fld_dst, struct t_stat *stat, double *bufr, int max);

#endif
