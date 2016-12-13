#!/bin/bash
#PBS -S /bin/bash
#PBS -J 0-1535
#PBS -l select=1:ncpus=1:mem=6gb
#PBS -l walltime=24:00:00
#PBS -o logs.osnpcall
#PBS -e logs.osnpcall
 

cd $PBS_O_WORKDIR

refFol="finalLikelihoodScore_swedes243"

posFiles=(`ls 0*/$refFol/*pos.txt`)
num=${#posFiles[@]}
outFol="pairWiseDiff"

mkdir $outFol

for (( i=0; i<$num ;i=i+1 ));do
#	echo "${posFiles[$i]}	${posFiles[$j]}"
	bash /home/GMI/rahul.pisupati/MyScripts/CompareVCFtoDataAcc/11_CompareVCFtoDataAcc/04_CalcHammingDist_twoPOSfiles.sh $PBS_O_WORKDIR/${posFiles[$PBS_ARRAY_INDEX]} $PBS_O_WORKDIR/${posFiles[$i]} >> $PBS_O_WORKDIR/$outFol/plate_${PBS_ARRAY_INDEX}.pairwise.txt &
	if [ $(($i % 16)) = 0 ]
	then
		wait
	fi
done
