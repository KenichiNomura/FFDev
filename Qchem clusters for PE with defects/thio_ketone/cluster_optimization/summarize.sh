#! /bin/bash

## ------------------------------------------
## Run this script inside each cluster folder
## ------------------------------------------

# Create label
echo $PWD | sed 's/.*\///g' | awk '{print "Summary for cluster : " $0}' > summary.txt
echo "" >> summary.txt
echo "" >> summary.txt

# Print cluster energy
grep Convergence CHOSx.out | tail -1 | awk '{printf "The final cluster energy is %19.16f Hartrees = %19.16f eV = %19.16f kCal \n", $2, $2*27.2107, $2*627.503}' >> summary.txt
echo "" >> summary.txt
echo "" >> summary.txt


# Print Mulliken charges
wc -l mol.pdb | awk '{printf "grep -A %d Mulliken CHOSx.out | tail -n %d \n", $1+4, $1+5}' | bash >> summary.txt
echo "" >> summary.txt
