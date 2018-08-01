/********************************************************
 * Created by Yongjun ZHENG on 01 Oct 2016.             *
 * Email: yongjun.zheng@meteo.fr                        *
 * Copyright © 2016 Yongjun ZHENG. All rights reserved. *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <dirent.h>
#include "utils.h"

int main(int argc, char *argv[]) {
    int icpu, ncpu, shmid, mpi_proc_null = 0;
    int iproc, nproc;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    if (argc >= 3) {
        system("for shmid in `ipcs -m | awk 'NR>3{print $2}'`\ndo\nipcrm -m $shmid\ndone\n");

        char path[256];
        strcpy(path, argv[1]);
        ncpu = atoi(argv[2]);
        if (argc == 4) mpi_proc_null = atoi(argv[3]);

        // get key
        key_t key = ftok(path, iproc);
        if (key == (key_t)-1) {
            perror("Error in ftok: ");
            return -1;
        }

        size_t shmsz = sizeof(struct t_config) + sizeof(struct t_domain) * ncpu + sizeof(struct t_stat) * ncpu;

        // create shared memory
        shmid = shmget(key, shmsz, IPC_CREAT | 0600);
        if (shmid == -1) {
            perror("Error in shmget for creating memory: ");
            return -1;
        }

        // attach
        char *GO_MBLOCK = (char *)shmat(shmid, NULL, 0);
        if (GO_MBLOCK == (void *)-1) {
            perror("Error in shmat: ");
            return -1;
        }

        // initialize the shared memory
        struct t_config *config;
        config = (struct t_config *)GO_MBLOCK;
        read_config("./config.xml", config);

        struct t_domain *dm_z = NULL;
        struct t_domain *da_z = NULL;
        da_z = (struct t_domain *)(GO_MBLOCK + sizeof(struct t_config));
        for (icpu = 0; icpu < ncpu; icpu++) {
            dm_z = da_z + icpu;
            // domain partition
            dm_z->ntile_x = config->ntile_x, dm_z->ntile_y = config->ntile_y;
            domain_partition('Z', icpu, ncpu, config, dm_z);

            if (mpi_proc_null) {
                if (dm_z->w_neighbor < 0) dm_z->w_neighbor = mpi_proc_null;
                if (dm_z->e_neighbor < 0) dm_z->e_neighbor = mpi_proc_null;
                if (dm_z->s_neighbor < 0) dm_z->s_neighbor = mpi_proc_null;
                if (dm_z->n_neighbor < 0) dm_z->n_neighbor = mpi_proc_null;
                if (dm_z->d_neighbor < 0) dm_z->d_neighbor = mpi_proc_null;
                if (dm_z->u_neighbor < 0) dm_z->u_neighbor = mpi_proc_null;
            }
        }

        if (0) { // turn on if want to check shared memory
            FILE *fid = fopen("mblock.bin", "wb");
            fwrite(config, sizeof(struct t_config), 1, fid);
            fwrite(da_z, sizeof(struct t_domain), ncpu, fid);
            fclose(fid);
        }

        // detach
        if (shmdt(GO_MBLOCK) == -1) {
            perror("Error in shmdt: ");
            return -1;
        }

        fprintf(stderr, "CREATE(node=%d): shmid=%d, size=%d, MPI_PROC_NULL=%d => %d\n", iproc, shmid, shmsz, MPI_PROC_NULL, mpi_proc_null);

        int *shmid_i = (int *)malloc(sizeof(int) * nproc);
        MPI_Gather(&shmid, 1, MPI_INT, shmid_i, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (iproc == 0) {
            sprintf(path, "%s/GO_MBLOCK_SHMID", argv[1]);
            FILE *fp = fopen(path, "wt");
            for (int k=0; k<nproc; k++) {
                fprintf(fp, "%d:", shmid_i[k]);
            }
            fclose(fp);
        }
        free(shmid_i);
    }else if (argc == 2) {
    	char *token = strtok(argv[1], ":");
    	for (int k=0; k<iproc; k++) {
    		token = strtok(NULL, ":");
    	}
        shmid = atoi(token);

        // delete shared memory
        if (shmctl(shmid, IPC_RMID, NULL) == -1) {
            perror("Error in shmctl: ");
            return -1;
        }

        fprintf(stderr, "DELETE(node=%d): shmid=%d\n", iproc, shmid);
        system("ipcs -m");
    }else{
        fprintf(stderr, "Usage: %s shmid | path ncpu\n", argv[0]);
        return -1;
    }

    MPI_Finalize();
    return 0;
}
