#!/bin/bash
. /home/9yelin9/.bash_profile

#$ -q openmp.q@phase09
#$ -pe mpi 16
#$ -j y
#$ -cwd
#$ -o log/$JOB_NAME.log

t0=$(date +%s.%N)
t0_string=$(date)

./mc 30 16 1e5 1e6 1.2 0.05 0.01 0.5

t1=$(date +%s.%N)
t1_string=$(date)

t=$(echo "$t1 - $t0"|bc)
h=$(echo "($t/3600)"|bc)
m=$(echo "($t%3600)/60"|bc)
s=$(echo "($t%3600)%60"|bc)

echo ""
echo "# Job ID       : $JOB_ID"
echo "# Job name     : $JOB_NAME"
echo "# Time start   : $t0_string"
echo "# Time end     : $t1_string"
echo "# Time elapsed : ${h}h ${m}m ${s}s"
echo ""
