/********************************************************
 * Created by Yongjun ZHENG on 05 Sep 2016.             *
 * Email: yongjun.zheng@meteo.fr                        *
 * Copyright © 2016 Yongjun ZHENG. All rights reserved. *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coll.h"

void find_root_leaf(int icpu, int ncpu, int *root, int *left, int *right, int *ncpu_left, int *ncpu_right) {
    *root = -1;
    int curr_cpu = 0;
    while (1) {
        *ncpu_left = ncpu / 2;
        *ncpu_right = ncpu - *ncpu_left - 1; // number of cpu in right branch = ncpu - ncpu_left - 1 (root)
        *left =  curr_cpu + 1;
        *right = *left + *ncpu_left;
        if (icpu == curr_cpu) {
            if (ncpu == 2) *right = -1, *ncpu_right = 0;
            if (ncpu == 1) *left = *right = -1, *ncpu_left = *ncpu_right = 0;
            break;
        }else if (icpu < *right) {
            *root = curr_cpu;
            curr_cpu = *left;
            ncpu = *ncpu_left;
        }else {
            *root = curr_cpu;
            curr_cpu = *right;
            ncpu = *ncpu_right; 
        }
    }
}

int bcast(void *bufr, int count, MPI_Datatype data_type, MPI_Comm comm) {
    int ierr, icpu, ncpu;
    
    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
    
    // width first traversal
    icpu++;
    int root = (icpu >> 1) - 1;
    int left = (icpu << 1) - 1, right = icpu << 1;
    icpu--;
    
    if (root >= 0) ierr = MPI_Recv(bufr, count, data_type, root, 999, comm, MPI_STATUS_IGNORE);
    if (left < ncpu) ierr = MPI_Send(bufr, count, data_type, left, 999, comm);
    if (right < ncpu) ierr = MPI_Send(bufr, count, data_type, right, 999, comm);
    
    return ierr;
}

int scatter(void *bufr, int count, MPI_Datatype data_type, MPI_Comm comm) {
    int ierr, icpu, ncpu;
    
    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
    
    int data_size;
    ierr = MPI_Type_size(data_type, &data_size);
    data_size *= count;
    
    int root, left, right, n_left, n_right, n_subtree;
    find_root_leaf(icpu, ncpu, &root, &left, &right, &n_left, &n_right);
    n_subtree = n_left + 1 + n_right;
    n_left *= data_size;
    n_right *= data_size;
    n_subtree *= data_size;
    
    void *b_left = NULL, *b_right = NULL, *b_subtree = NULL; // temporal buffer
    
    // allocate memory and receive from root
    if (root >= 0) {
        b_subtree = (void *)malloc(n_subtree);
        ierr = MPI_Recv(b_subtree, n_subtree, MPI_CHAR, root, 999, comm, MPI_STATUS_IGNORE);
        memcpy(bufr, b_subtree, data_size);
    }else{
        b_subtree = bufr;
    }

    // send to left 
    if (left >= 0) {
        b_left = (char *)b_subtree + data_size;
        ierr = MPI_Send(b_left, n_left, MPI_CHAR, left, 999, comm);
    }
    
    // send to right
    if (right >= 0) {
        b_right = (char *)b_subtree + data_size + n_left;
        ierr = MPI_Send(b_right, n_right, MPI_CHAR, right, 999, comm);
    }
    
    // deallocate memory
    if (root >= 0) {
        free(b_subtree);
    }

    return ierr;
}

int gather(void *bufr, int count, MPI_Datatype data_type, MPI_Comm comm) {
    int ierr, icpu, ncpu;
    
    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
    
    int data_size;
    ierr = MPI_Type_size(data_type, &data_size);
    data_size *= count;
    
    int root, left, right, n_left, n_right, n_subtree;
    find_root_leaf(icpu, ncpu, &root, &left, &right, &n_left, &n_right);
    n_subtree = n_left + 1 + n_right;
    n_left *= data_size;
    n_right *= data_size;
    n_subtree *= data_size;
    
    void *b_left = NULL, *b_right = NULL, *b_subtree = NULL; // temporal buffer
    
    if (root >= 0) {
        b_subtree = (void *)malloc(n_subtree);
        memcpy(b_subtree, bufr, data_size);
    }else{
        b_subtree = bufr;
    }

    // receive from left 
    if (left >= 0) {
        b_left = (char *)b_subtree + data_size;
        ierr = MPI_Recv(b_left, n_left, MPI_CHAR, left, 999, comm, MPI_STATUS_IGNORE);
    }
    
    // receive from right
    if (right >= 0) {
        b_right = (char *)b_subtree + data_size + n_left;
        ierr = MPI_Recv(b_right, n_right, MPI_CHAR, right, 999, comm, MPI_STATUS_IGNORE);
    }
    
    // send to root
    if (root >= 0) {
        ierr = MPI_Send(b_subtree, n_subtree, MPI_CHAR, root, 999, comm);
    
        // deallocate memory
        free(b_subtree);
    }
    return ierr;
}

int allgather(void *send_bufr, int count, MPI_Datatype data_type, void *recv_bufr, MPI_Comm comm) {
    int ierr, icpu, ncpu;
    
    MPI_Comm_size(comm, &ncpu);
    MPI_Comm_rank(comm, &icpu);
    
    int data_size;
    ierr = MPI_Type_size(data_type, &data_size);
    data_size *= count;
    memcpy((char *)recv_bufr + data_size * icpu, send_bufr, data_size);
    
    int inext, iprev, isrc_next, isrc_prev, n;
    inext = (icpu + 1) % ncpu;
    iprev = (icpu - 1 + ncpu) % ncpu;
    for (n = 1; n < ncpu; n++) {
        // the id of original processor of the message which will send to next processor
        isrc_next = (icpu - n + 1 + ncpu) % ncpu;
        // the id of original processor of the message which will receive from previous processor
        isrc_prev = (icpu - n + ncpu) % ncpu;
        
        // send and receive
        ierr = MPI_Sendrecv((char *)recv_bufr + data_size * isrc_next, count, data_type, inext, 1,
                            (char *)recv_bufr + data_size * isrc_prev, count, data_type, iprev, 1,
                            comm, MPI_STATUS_IGNORE);
    }
    
    return ierr;
}
