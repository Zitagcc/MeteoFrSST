#include "mpp.h"
#define PI 3.14159265358979323846264338327950288419716939937510L

void print_time(const char *prefix, int icpu);
void read_config(const char *filename, struct t_config *config);
void fd_derivative(int ntrunc, char dir, double dx, struct t_domain *dm, double *fld, double *deriv);
void output_one(struct t_domain *dm, double *fld);
void output_all(int icpu, struct t_domain *da, double *fld);
//void do_stat(struct t_stat *stat_avg, struct t_stat *stat_max, struct t_stat *stat_min, struct t_stat *stat, int navg);
void output_stat(int icpu, int ncpu, struct t_stat *stat, FILE *fp);
