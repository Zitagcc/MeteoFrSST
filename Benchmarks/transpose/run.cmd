#!/bin/bash
smpicxx alltoall.cc -o alltoall.exe
smpirun -np 1000000 -platform 01000000/platform.xml -hostfile 01000000/hostfile ./alltoall.exe --cfg=smpi/simulate-computation:no --cfg=smpi/host-speed:2500000000 --cfg=network/model:SMPI --cfg=contexts/stack-size:64 --cfg=contexts/guard-size:0 1 8000 4
