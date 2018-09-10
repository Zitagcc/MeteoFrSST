#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 sst|simgrid|clean [_ALGO_ _STAT_ | topology]"
    exit
fi

if [ $# -eq 1 -a "$1" == "clean" ]; then
  make -f halo.mk clean
  make -f shm_mgr.mk clean

  rm -f *.o *.exe nohup.out
  exit
fi

if [ $# -eq 3 ]; then
  if [ `uname` == "Linux" ]; then
    BK=""
  else
    BK="_bk"
  fi

  if [ "$1" == "sst" ]; then
    sed -i $BK "s/CXX := .*/CXX := sst++/g" halo.mk
    sed -i $BK "s/CC := .*/CC := sstcc/g" halo.mk
  else if [ "$1" == "smpi" ]; then
    sed -i $BK "s/CXX := .*/CXX := smpicxx/g" halo.mk
    sed -i $BK "s/CC := .*/CC := smpicc/g" halo.mk
  else
    sed -i $BK "s/CXX := .*/CXX := mpicxx/g" halo.mk
    sed -i $BK "s/CC := .*/CC := mpicc/g" halo.mk
  fi

  sed -i $BK "s/ALGO := .*/ALGO := $2/g" halo.mk
  sed -i $BK "s/COMM := .*/COMM := $3/g" halo.mk

  make -f halo.mk clean
  make -f halo.mk
  make -f shm_mgr.mk clean
  make -f shm_mgr.mk

  exit
fi


rm -f nohup.out
if [ "$1" == "sst" ]; then
  if [ "$2" == "torus" ]; then
    nohup ./JOB_BAT_SST.py 28800:14400:200000 100 30 16:0 first_available:block torus:minimal 1:1 10:10:10:10:10 1000:1000:1000:1000:1000 25:25:25:25 &
  elif [ "$2" == "fattree" ]; then
    nohup ./JOB_BAT_SST.py 28800:14400:200000 100 30 16:0 first_available:block fattree:fattree 1:1 10:10:10:10:10 1000:1000:1000:1000:1000 4:25:25 &
  elif [ "$2" == "dragonfly" ]; then
    nohup ./JOB_BAT_SST.py 28800:14400:1000000 100 30 16:0 first_available:block dragonfly:minimal 1:1 10:10:10:10:10 1000:1000:1000:1000:1000 25:25:75:25:74 &
  fi
else
  nohup ./JOB_BAT_SST.py 28800:14400:200000 100 30 16:0 first_available:block torus:minimal 1:1 10:10:10:10:10 1000:1000:1000:1000:1000 25:25:25:25 &
  # nohup ./JOB_BAT_SST.py 28800:14400:10000 100 30 200 &
fi
