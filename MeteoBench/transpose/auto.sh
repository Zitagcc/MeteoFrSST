#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 sst|clean [_ALGO_ _COMM_ | topology]"
    exit
fi

if [ $# -eq 1 -a "$1" == "clean" ]; then
  make -f transpose.mk clean
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
    sed -i $BK "s/CXX := .*/CXX := sst++/g" transpose.mk
    sed -i $BK "s/CC := .*/CC := sstcc/g" transpose.mk
  else
    sed -i $BK "s/CXX := .*/CXX := mpicxx/g" transpose.mk
    sed -i $BK "s/CC := .*/CC := mpicc/g" transpose.mk
  fi

  sed -i $BK "s/ALGO := .*/ALGO := $2/g" transpose.mk
  sed -i $BK "s/COMM := .*/COMM := $3/g" transpose.mk

  make -f transpose.mk clean
  make -f transpose.mk
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
    nohup ./JOB_BAT_SST.py 28800:14400:200000 100 30 16:0 first_available:block dragonfly:minimal 1:1 10:10:10:10:10 1000:1000:1000:1000:1000 25:25:25:25:24 &
  fi
else
  nohup ./JOB_BAT_SST.py 28800:14400:10000 100 30 200 &
fi
