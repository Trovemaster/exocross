#!/bin/bash -l
#
# Check if our working directory is on the central file server
#

export exec=j-xcross-1212.x

#exec=j-trove_1105.x


export pwd=`pwd`
echo $pwd

export name=`echo $1 | sed -e 's/\.inp//'`
echo $name

#export TIME=$3
echo $JOB

if [ -e "${name}.o" ]; then
   /bin/rm ${name}.o
fi

if [ -e "${name}.e" ]; then
   /bin/rm ${name}.e
fi

if [ -e "${name}.out" ]; then
  if [ -e "${name}.tmp" ]; then
    /bin/rm ${name}.tmp
  fi
  /bin/mv ${name}.out ${name}.tmp
fi


export TIME="11:59:00"

export MEM=5G
export nproc=$2


qsub -P AllUsers -N $name -e $name.e -o $name.o -l h_rt=$TIME -l mem=$MEM -pe smp $nproc -l tmpfs=100G  -wd $pwd $pwd/run_script.sh $nproc $name $exec $pwd


