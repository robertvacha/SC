#!/bin/bash

name=$1
topo_old="top.init_old"
topo_new="top.init_new"
samples=$2  # number of uniformly distributed points on a sphere -> used as dir for spherocylinder
angle_inc=$3  # 3.14rad = 180deg, 0.2rad ~ 11.5deg
coor_inc=$4

var="0"

for x in `seq 1.0 $coor_inc 9.0`; do   #      dir oriented in X
  for y in `seq 1.0 $coor_inc 9.0`; do # patchdir oriented in Y
    for z in `seq 1.0 $coor_inc 9.0`; do

      ### generate config
      cp config.init_prescript config.init
      ./Spc_config_generator $x $y $z $samples $angle_inc >> config.init
      head -n -1 config.init > config.init2
      num=`tail -n 1 config.init` ## number of particles for topo
      mv config.init2 config.init  

      ### generate topo
      cp $topo_new top.init
      printf $num >> top.init

      ### Run new
      ./scOOP | grep -A 99999999 "Production run:" | tail -n +3 > res
   
      ### generate topo
      cp $topo_old top.init
      printf $num >> top.init

      ### run old
      ./sc35 2> /dev/null | grep -A 99999999 "Production run:" | tail -n +3 > res2

      ### Print resutls
      tail -n +3 config.init > temp
      paste res res2 temp -d" " | awk '{if($1 != $2) print}' >> mismatch_$name

      let "var=var+num"

    done  
  done
done

printf "%s sampled configurations " $var

if [ -s "mismatch_$name" ] 
then
	exit 1
else
	exit 0
fi
