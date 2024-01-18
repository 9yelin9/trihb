#!/bin/bash
. /home/9yelin9/.bash_profile

#$ -q openmp.q@phase09
#$ -pe mpi 1
#$ -j y
#$ -cwd
#$ -o log/h1.2000_D0.0000_T-1.0000_vs0.5000_ve0.0200_dv-0.0200_L30_M1e5_tm11011519/$JOB_NAME.log

t0=$(date +%s.%N)
t0_string=$(date)

./mc none eq data/h1.2000_D0.0000_T-1.0000_vs0.5000_ve0.0200_dv-0.0200//L30_M1e5_tm11011519 30 1e5   1.200000   0.000000   0.480000

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
