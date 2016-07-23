#!/bin/bash

#SBATCH -A A-asoz
#SBATCH -J reformat1           # job name
#SBATCH -o reformat1_out  # output and error file name (%j expands to jobID)
#SBATCH -n 1              # total number of mpi tasks requested
#SBATCH -p serial     # queue (partition) -- normal, development, etc.
#SBATCH -t 10:00:00        # run time (hh:mm:ss) - 1.5 hours

#
# This script combines nc_node^3 clumping grid data files to one c_all file.
#
# Inputs
#   1. ./so/nc<nc>/grid_files created by clumping_tree: for clumping data
#   2. redshifts_ref.txt file: one redshift per line
#
# Outputs:
#   global/so/<redshift>c_all.dat
#
#
# If output already exists with correct size, it will be skipped
#

# Configuration

nccpu=12            # nccpu^3 MPI nodes used in cubep3m
ncs="25 50 100"     # ngrid per dim per MPI node

reformat=reformat   # full path of the reformat executable (if necessary)
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

      if [ ! -d ../../global/$dir/nc$ncg ]; then
        mkdir ../../global/$dir/nc$ncg
      fi

      if [ -f ../../global/$dir/nc$ncg/${z}c_all.dat ]; then
	size=`ls -l ../../global/$dir/nc$ncg/${z}c_all.dat | awk '{print $5}'`
	sizeex=`expr 12 \* $ncg \* $ncg \* $ncg + 12`
	if [ $size -ne $sizeex ]; then
	  echo "SIZE $dir/nc$ncg/${z}v_all.dat is incorrect."
	  echo "$size != $sizeex"
	  echo "reformat $nccpu ../../global/$dir/nc$ncg $z"
	  date

	  # c_all file already exist but incorrect file size
	  $reformat $nccpu ../../global/$dir/nc$ncg ${z}c 1
	  if [ $? -ne 0 ]; then
	      echo "Error reformat2"; exit
	  fi
	fi
      else
	# c_all does not exist
	echo "reformat $nccpu ../../global/$dir/nc$ncg $z"
	date
	$reformat $nccpu ../../global/$dir/nc$ncg ${z}c 1
	if [ $? -ne 0 ]; then
	  echo "Error reformat"; exit
	fi
      fi

      cd ..
    done
    cd ..
  done
done

echo "c_all all done."

