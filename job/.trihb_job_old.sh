#!/bin/bash
. /home/9yelin9/.bash_profile

#$ -q openmp.q@phase09
#$ -pe mpi 1
#$ -j y
#$ -cwd
#$ -o log/$JOB_NAME.log

t0=$(date +%s.%N)
t0_string=$(date)

dn="h-1.0000_D0.0500_T0.0001_vs0.0000_ve10.0000_dv1.0000"
fn="data/${dn}/script.txt"
L="12"
M="1e4"
v=$(date +"%m%d%H%M")
init="none"
init_mode="eq"

save="data/${dn}/L${L}_M${M}_v${v}"
info="${save}/info.txt"

if [ ! -d "$save" ]; then
	mkdir $save
fi
printf "script ${fn}\nL ${L}\nM ${M}\nv ${v}\ninit ${init}\ninit_mode ${init_mode}\n" > $info

while IFS="" read -r line || [ -n "$line" ]
do
	./mc $init $init_mode $save $L $M $line
done < $fn

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
