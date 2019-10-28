#!/bin/bash

# Lancement: qsub -o .output -j y run_multi_cl.sh

#$ -N morpheus
#$ -wd /workdir2/auder/morpheus/reports
#$ -m abes
#$ -M benjamin@auder.net
#$ -pe make 50
#$ -l h_vmem=1G
rm -f .output

module load R/3.6.1

N=100
n=1e5
nc=50
nstart=5

for d in 2 5 10; do
	for link in "logit" "probit"; do
		R --slave --args N=$N n=$n nc=$nc d=$d link=$link nstart=$nstart <multistart.R >out_${n}_${link}_${d}_${nstart} 2>&1
	done
done
