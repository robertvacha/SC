#!/bin/bash


for name in `ls -d test_*`; do
    cd $name/new
    rm -f config.last config.last2 energy.dat stat.dat movie wl-new.dat cluster_stat.dat cluster.dat
    cd ../old
    rm -f config.last config.last2 energy.dat stat.dat movie wl-new.dat cluster_stat.dat cluster.dat
    cd ../..
done

cd volumeChange
for name in `ls -d ?h ?l`; do
    cd $name/new
    rm -f config.last config.last2 energy.dat stat.dat movie wl-new.dat cluster_stat.dat cluster.dat
    cd ../old
    rm -f config.last config.last2 energy.dat stat.dat movie wl-new.dat cluster_stat.dat cluster.dat
    cd ../..
done
cd ..
