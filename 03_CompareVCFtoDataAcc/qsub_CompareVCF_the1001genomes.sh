#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N VCFcompare
#PBS -V
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -o log

#Input:  VCF files generated after SNP calling
#Output: Compare accessions

cd $PBS_O_WORKDIR

inFiles=(`ls *.vcf`)
length=${#inFiles[@]}
count=0
for (( i=0; i<$length; i=i+1 ));do
	bash ~/MyScripts/CompareVCFtoDataAcc/03_CompareVCFtoDataAcc/potatoSkin_VCFfile_CompareSNPs.sh ${inFiles[$i]} 100 &
	count=$((count+1))
	if [ $(($count % 16)) = 0 ]
	then
		wait
	fi
done
wait
