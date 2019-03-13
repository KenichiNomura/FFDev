#!/bin/sh
#SBATCH -n 12 
#SBATCH -t 14400 
#SBATCH -p priya 
#SBATCH -A lc_pv

source /usr/usc/python/3.6.0/setup.sh 
source /home/rcf-12/knomura/hpc-23/ReaxFFDev/FFDev/Tools/pescan/setup.sh

PESCAN=../../pescan.py

dir=$1
echo "dir : $1"
for f in `find $dir -name "*.xyz"` ; do 
    echo "python3.6 -u ${PESCAN} ${f} > ${f}.log"
    python3.6 -u ${PESCAN} ${f} > ${f}.log
done
