#!/bin/bash

#SBATCH -A A-asoz
#SBATCH -J reformat_nv
#SBATCH -o reformat_out
#SBATCH -n 1
#SBATCH -p serial
#SBATCH -t 10:00:00

#
# This script combines nc_node^3 grid data files to one n_all, ntot_all, v_all.
#
# Inputs
#   1. ./so/nc<nc>/grid_files created by clumping_tree: n, h, vel data
#   2. redshifts_ref.txt file: one redshift per line
#
# Outputs:
#   global/so/<redshift>n_all.dat, ntot_all.dat, v_all.dat
#
#
# If output already exists with correct size, it will be skipped
#


# Configuration
nccpu=12            # nccpu^3 MPI nodes used in cubep3m
ncs="25 50 100"     # ngrid per dim per MPI node

reformat2=reformat2 # full path of the reformat2 executable (if necessary)
dirs="so"           # "so fof" if you also have fof dir




if [ ! -d global ]; then
  mkdir global
fi

for z in `cat redshifts_ref.txt`
do
  echo "redshift $z"
  for dir in $dirs
  do
    if [ ! -d global/$dir ]; then
      mkdir global/$dir
    fi

    cd $dir

    for nc in $ncs
    do
      cd nc$nc
      ncg=`expr $nc \* $nccpu`

      #echo "for $dir $nc"

      if [ ! -d ../../global/$dir/nc$ncg ]; then
        mkdir ../../global/$dir/nc$ncg
      fi

      if [ -f ../../global/$dir/nc$ncg/${z}v_all.dat ]; then
	size=`ls -l ../../global/$dir/nc$ncg/${z}v_all.dat | awk '{print $5}'`
	sizeex=`expr 12 \* $ncg \* $ncg \* $ncg + 12`
	if [ $size -ne $sizeex ]; then
	  echo "SIZE $dir/nc$ncg/${z}v_all.dat is incorrect."
	  echo "$size != $sizeex"
	  echo "reformat2 $nccpu ../../global/$dir/nc$ncg $z"
	  date
	  $reformat2 $nccpu ../../global/$dir/nc$ncg $z .
	  if [ $? -ne 0 ]; then
	      echo "Error reformat2"; exit
	  fi
	fi
      else
	echo "reformat2 $nccpu ../../global/$dir/nc$ncg $z"
	date
	$reformat2 $nccpu ../../global/$dir/nc$ncg $z .
	if [ $? -ne 0 ]; then
	  echo "Error reformat2"; exit
	fi
      fi

      cd ..
    done
    cd ..
  done
done

echo "global/ n_all, ntot_all, v_all all done."
