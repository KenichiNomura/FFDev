#!/bin/bash
#SBATCH --ntasks=12 
#SBATCH --account=lc_pv
#SBATCH --partition=priya
#SBATCH --time=23:00:00
#SBATCH --export=none

source /usr/usc/intel/default/setup.sh
source /usr/usc/openmpi/default/setup.sh.intel

export QC=/usr/usc/qchem/5.0
export QCAUX=$QC/qcaux
export QCPLATFORM=LINUX_Ix86_64
export QCRSH=ssh
export PATH=$QC/bin:$PATH
export QCSCRATCH=$TMPDIR

ulimit -s unlimited

echo "starting simulation **************************************"
date
qchem -save -nt 12 Cl.in Cl.out Cl.save
date
echo "simulation finished **************************************"
echo
