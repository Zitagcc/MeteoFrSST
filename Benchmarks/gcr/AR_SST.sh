#!/bin/bash

if [ $# -ne 2 ]; then
  echo "Usage: $0 sst topology"
  exit
fi

rm -f allreduce.o allreduce.exe nohup.out
sst++ -fPIC -D_NOMEM_ -c allreduce.cc -o allreduce.o
sst++ allreduce.o -o allreduce.exe

if [ "$1" == "sst" ]; then
  if [ "$2" == "torus" ]; then
    nohup ./AR_SST.py 1:1:8:1000000 100 30 16:0 first_available:block torus:minimal 1:1 10:10:10:10:10 1000:1000:1000:1000:1000 75:25:25:25 &
  elif [ "$2" == "fattree" ]; then
    nohup ./AR_SST.py 1:1:8:1000000 100 30 16:0 first_available:block fattree:fattree 1:1 10:10:10:10:10 1000:1000:1000:1000:1000 4:33:33 &
  elif [ "$2" == "dragonfly" ]; then
    nohup ./AR_SST.py 1:1:8:1000000 100 30 16:0 first_available:block dragonfly:minimal 1:1 10:10:10:10:10 1000:1000:1000:1000:1000 25:25:75:25:74 &
  fi
elif [ "$1" == "simgrid" ]; then
  if [ "$2" == "torus" ]; then
    nohup ./JOB_BAT_SIM.py 28800:14400:200000 100 30 16 N TORUS:390625 1:1 10:10:10000 1:1:1 '25,25,25,25' &
  elif [ "$2" == "fattree" ]; then
    nohup ./JOB_BAT_SIM.py 28800:14400:200000 100 30 16 N FAT_TREE:390625 1:1 10:10:10000 1:1:1 '4;25,25,25,25;1,25,25,25;2,2,2,2' &
  elif [ "$2" == "dragonfly" ]; then
    nohup ./JOB_BAT_SIM.py 28800:14400:200000 100 30 16 N DRAGONFLY:390625 1:1 10:10:10000 1:1:1 '25,1;25,1;25,1;25' &
  fi
else
  nohup ./AR_SST.py 1:2:4:10000 100 30 200 &
fi
