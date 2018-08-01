#include "mpp.h"
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

void print_time(const char *prefix, int icpu);
void read_config(const char *filename, struct t_config *config);
void fft_derivative(char dir, struct t_domain *dm, double *fld, double *deriv, double L);
void output_one(struct t_domain *dm, double *fld);
void output_all(int icpu, struct t_domain *da, double *fld);
//void do_stat(struct t_stat *stat_avg, struct t_stat *stat_max, struct t_stat *stat_min, struct t_stat *stat, int navg);
void output_stat(int icpu, int ncpu, struct t_stat *stat, FILE *fp);
