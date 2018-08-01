#ifndef __COLL__
#define __COLL__

#include <mpi.h>

int bcast(void *bufr, int count, MPI_Datatype data_type, MPI_Comm comm);
int scatter(void *bufr, int count, MPI_Datatype data_type, MPI_Comm comm);
int gather(void *bufr, int count, MPI_Datatype data_type, MPI_Comm comm);
int allgather(void *send_bufr, int count, MPI_Datatype data_type, void *recv_bufr, MPI_Comm comm);
#endif