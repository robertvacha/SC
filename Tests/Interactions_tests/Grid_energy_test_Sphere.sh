#!/bin/bash

name=$1
topo_old="top.init_old"
topo_new="top.init_new"
coor_inc=$2

#for x in `seq 1.0 $coor_inc 9.0`; do   #      dir oriented in X
#  for y in `seq 1.0 $coor_inc 9.0`; do # patchdir oriented in Y

    x="5.0"
    y="5.0"

    num="0"
    cp config.init_prescript config.init

    for z in `seq 1.0 $coor_inc 5.0`; do
      printf "%s %s %s %s\n" $x $y $z " 1 0 0  0 1 0  0 0" >> config.init
      let "num++"
    done  

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

#  done
#done

printf "%s sampled configurations " $var

if [ -s "mismatch_$name" ] 
then
	exit 1
else
	exit 0
fi
