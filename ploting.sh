#!/bin/bash
#SBATCH --job-name="halo.exe"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=normal64

#SBATCH --exclusive
#SBATCH --output=00000100/slurm-%j.out

# if [ -f ${HOME}/.bashrc ]; then
#     ${HOME}/.bashrc
# fi

cd Benchmarks/halo 
./auto.sh 1 2 3
smpirun -np 1 -platform ../cluster_crossbar.xml -hostfile ../cluster_hostfile.txt --cfg=smpi/host-speed:100 ./shm_mgr.x 00000100 100 -1
GO_MBLOCK_SHMID=`cat 00000100/GO_MBLOCK_SHMID`
smpirun -np 1  -platform ../cluster_crossbar.xml  -hostfile ../cluster_hostfile.txt --cfg=smpi/host-speed:100 ./halo.exe -s $GO_MBLOCK_SHMID  -f 00000100/parameters.ini
smpirun -np 1  -platform ../cluster_crossbar.xml  -hostfile ../cluster_hostfile.txt --cfg=smpi/host-speed:100 ./shm_mgr.x $GO_MBLOCK_SHMID

