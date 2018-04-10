#!/bin/bash

function print {
    if [ $? -eq 0 ]
    then
        echo -e $1 " \e[1;32mOK\e[0m."
    else 
        echo -e $1 " \e[1;31mNOK\e[0m."
    fi
}

echo "Using testing versions of SC program"
echo "TEST: generate configuration"
echo "Writeout - Coordinates energy_new_version energy_version_10.1063_1.4933229_fixed_repulse"
echo "Only energy between 1st particle and rest calculated and written out"
echo "Config can have overlaps"
echo "Too high energy truncated to 1000.0"
echo "Cases when energy of both versions match is ommited in output"

g++ main.cpp -std=c++11 -o Spc_config_generator

gcc sc35_fixedRepulse.c -lm -o sc35

rm -f mismatch_*

echo ""
echo "Testing..."


nameA="A 1 "
nameB="B 2 "

sphere=("SPN 1.33 0.95" "SPA 1.33 1.2 1.346954458 1.6")

spherocylinder_n=("PSC 1.33 1.2 1.346954458 0.3 30.0 0.0 3 0.0" "CPSC 1.33 1.2 1.346954458 0.3 30.0 0.0 3 0.0" "CHPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 0.0 5.0" "CHCPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 0.0 5.0" "TPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 0.0 180.0 90 5.0" "TCPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 0.0 180.0 90 5.0" "TCHPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 0.0 180.0 90 5.0 5.0" "TCHCPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 0.0 180.0 90 5.0 5.0")
spherocylinder_o=("PSC 1.33 1.2 1.346954458 0.3 30.0 0.0 3" "CPSC 1.33 1.2 1.346954458 0.3 30.0 0.0 3" "CHPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 5.0" "CHCPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 5.0" "TPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 180.0 90 5.0" "TCPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 180.0 90 5.0" "TCHPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 180.0 90 5.0 5.0" "TCHCPSC 1.333333 1.2 1.346954458 0.3 90 5.0 3 180.0 90 5.0 5.0")


##
## Spherocylinders: PSC, CPSC, CHPSC, CHCPSC, TPSC, TCPSC, TCHPSC, TCHCPSC
##
for ix in ${!spherocylinder_n[*]}; do
  for jx in ${!spherocylinder_n[*]}; do
    if [ $jx -ge $ix ]; then

      cat top.init_prescript > top.init_new
      echo $nameA ${spherocylinder_n[$ix]} >> top.init_new
      echo $nameB ${spherocylinder_n[$jx]} >> top.init_new
      cat top.init_postscript >> top.init_new

      cat top.init_prescript > top.init_old
      echo $nameA ${spherocylinder_o[$ix]} >> top.init_old
      echo $nameB ${spherocylinder_o[$jx]} >> top.init_old
      cat top.init_postscript >> top.init_old

      type_i=`printf ${spherocylinder_n[$ix]} | awk '{print $1}'`
      type_j=`printf ${spherocylinder_n[$jx]} | awk '{print $1}'`
      ./Grid_energy_test_Spherocylinder.sh "TEST_${type_i}_${type_j}" 10 1.0 3.0
      print "TEST_${type_i}_${type_j}"

    fi
  done
done

## Sphere - SPA-SPA, SPN-SPN, SPN-SPA
#echo "Number of items in original array: ${#sphere[*]}"
for ix in ${!sphere[*]}; do
  for jx in ${!sphere[*]}; do
    if [ $jx -ge $ix ]; then

      cat top.init_prescript > top.init_new
      echo $nameA ${sphere[$ix]} >> top.init_new
      echo $nameB ${sphere[$jx]} >> top.init_new
      cat top.init_postscript >> top.init_new

      cat top.init_prescript > top.init_old
      echo $nameA ${sphere[$ix]} >> top.init_old
      echo $nameB ${sphere[$jx]} >> top.init_old
      cat top.init_postscript >> top.init_old

      type_i=`printf ${sphere[$ix]} | awk '{print $1}'`
      type_j=`printf ${sphere[$jx]} | awk '{print $1}'`
      ./Grid_energy_test_Sphere.sh "TEST_${type_i}_${type_j}" 0.0033
      print "TEST_${type_i}_${type_j}"

    fi
  done
done

## Sphere - Spherocylinder
##
for ix in ${!sphere[*]}; do
  for jx in ${!spherocylinder_n[*]}; do
    if [ $jx -ge $ix ]; then

      cat top.init_prescript > top.init_new
      echo $nameA ${sphere[$ix]} >> top.init_new
      echo $nameB ${spherocylinder_n[$jx]} >> top.init_new
      cat top.init_postscript >> top.init_new

      cat top.init_prescript > top.init_old
      echo $nameA ${sphere[$ix]} >> top.init_old
      echo $nameB ${spherocylinder_o[$jx]} >> top.init_old
      cat top.init_postscript >> top.init_old

      type_i=`printf ${sphere[$ix]} | awk '{print $1}'`
      type_j=`printf ${spherocylinder_n[$jx]} | awk '{print $1}'`
      ./Grid_energy_test_Spherocylinder_Sphere.sh "TEST_${type_i}_${type_j}" 1.11
      print "TEST_${type_i}_${type_j}"

    fi
  done
done
exit








