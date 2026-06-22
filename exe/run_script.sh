#!/bin/bash -l
#!


export nproc=$1 
export name=$2 
export exec=$3
export pwd=$4


echo "name:" $name
echo "nproc:" $nproc
echo "pwd" $pwd
echo `pwd`

export wdir=$pwd
#!
export OMP_NUM_THREADS=$nproc
export MKL_NUM_THREADS=$nproc
#!
export KMP_LIBRARY=turnaround
export KMP_AFFINITY=disabled

export KMP_STACKSIZE=1gb
export OMP_NESTED=FALSE
#!
hostname
#!
cd   $wdir
echo -e "Changed directory to `pwd`.\n"

export LAUNCH=time  ###"dplace -x2"
export TMPDIR=$wdir
#!
echo "TMPDIR = " $TMPDIR
echo "USER = " $USER
echo "OMP_NUM_THREADS = " $OMP_NUM_THREADS
echo "wdir" $wdir
echo "wdir" $wdir
echo "OMP_NUM_THREADS=" $OMP_NUM_THREADS
cd $wdir
#!
echo $wdir
#!

$LAUNCH $pwd/$exec < $pwd/$name.inp > $pwd/$name.out


#!
echo "DONE"
