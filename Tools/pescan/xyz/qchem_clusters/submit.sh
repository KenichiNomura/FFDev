#!/bin/sh

for d in halogens/* hydroxy ketone selenium thio_ketone thiol; do
  sbatch run.sh $d
done
