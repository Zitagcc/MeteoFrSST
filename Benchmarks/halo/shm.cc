/********************************************************
 * Created by Yongjun ZHENG on 01 Oct 2016.             *
 * Email: yongjun.zheng@meteo.fr                        *
 * Copyright © 2016 Yongjun ZHENG. All rights reserved. *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <fcntl.h>
#include <errno.h>
#include <mpi.h>
#include "shm.h"

#if defined(_SHMX_)
#define SHM_KEY "/SHM"
#endif

// please make sure the mapping per node is not exceeded vm.max_map_count
int shmc(void **mem, int n, int id) {
#if defined(_SHMX_)
    int fd;
    char key[64];
    sprintf(key, "%s.%d", SHM_KEY, id);
    
    // create a new memory object
    if ((fd = shm_open(key, O_RDWR | O_CREAT, 0660)) == -1) {
        fprintf(stderr, "shm_open for malloc failed: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }
    
    // set the memory object's size
    if ((ftruncate(fd, n)) == -1) {
        fprintf(stderr, "ftruncate failed: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }
    
    // map the memory object
    *mem = (void *)mmap(0, n, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (*mem == MAP_FAILED) {
        fprintf(stderr, "mmap failed: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }
    
    // the memory object remains in the system after closing
    close(fd);
#else
    *mem = malloc(n);
    if (*mem == NULL) {
        fprintf(stderr, "Error: can not allocate memory of %d bytes\n", n);
        return EXIT_FAILURE;
    }
#endif
    
    return EXIT_SUCCESS;
}

int shmd(void **mem, int n, int id) {
#if defined(_SHMX_)
    int fd;
    char key[64];
    sprintf(key, "%s.%d", SHM_KEY, id);
    
    // get the memory object
    if ((fd = shm_open(key, O_RDONLY, 0660)) == -1) {
        fprintf(stderr, "shm_open for free failed: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }

    // unmap the memory
    if (munmap(*mem, n) == -1) {
        fprintf(stderr, "munmap failed: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }

    // the memory object remains in the system after closing
    close(fd);

    MPI_Barrier(MPI_COMM_WORLD);

    // remove a memory object
    shm_unlink(key);
#else
    free(*mem);
#endif

    return EXIT_SUCCESS;
}

int shm_malloc(void **mem, int shmid) {
    // attach
    *mem = (void *)shmat(shmid, NULL, 0);
    if (*mem == (void *)-1) {
        perror("Error in shmat: ");
        return -1;
    }
    
    return 0;
}

int shm_free(void **mem) {
    // detach
    if (shmdt(*mem) == -1) {
        perror("Error in shmdt: ");
        return -1;
    }else{
        return 0;
    }
}
