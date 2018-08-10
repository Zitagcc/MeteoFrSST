
rm -rf time.dat
for ncp in `seq 5  5 500`;
do

/usr/bin/time -a -f '%e' -o time.tmp smpirun -n $ncp -hostfile ../cluster_hostfile.txt -platform ../cluster_crossbar.xml --cfg=smpi/host-speed:100 ./halo.exe
time=`cat time.tmp`

echo "$ncp $time" >> time.dat
rm time.tmp
done
