#ifndef __SHM__
#define __SHM__

int shmc(void **mem, int n, int id); // create
int shmd(void **mem, int n, int id); // destroy

int shm_malloc(void **mem, int shmid);
int shm_free(void **mem);

#endif
