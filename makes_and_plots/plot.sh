#!/bin/bash
if [ $# -ne 1 ]; then
   echo 'incorrect number of commandline arguments'
   exit
fi
n=$1
script=gnuscript$n
#echo ' ' | cat >> script
gnuplot $script
