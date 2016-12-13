#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -J 1-96
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=12:00:00
#PBS -o logs.ovcfcompare
#PBS -e logs.evcfcompare


#-------------Variables
qthres=100

#Get the SNP positions into a file
cd $PBS_O_WORKDIR
mkdir $outFolder
mkdir $outFolder/logs

#inFile=`ls *.filter.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.vcf//'`
inFile=`ls *.filter.gvcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.gvcf//'`
#inFile=`ls *.pos.txt |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.pos\.txt//'`

rm $PBS_O_WORKDIR/$outFolder/$inFile.ScoreAcc.txt $PBS_O_WORKDIR/$outFolder/$inFile.refScore.txt

function VCFcompare {
#	grep -v "^#" $1.filter.vcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	grep -v "^#" $PBS_O_WORKDIR/$1.filter.gvcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	# Run the python script with matched SNPs
	module load numpy
	module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
	module load pygwas
	module load pandas
	python ~/MyScripts/CompareVCFtoDataAcc/08_CompareVCFtoDataAcc/01_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -p $PBS_O_WORKDIR/$outFolder/$1.pos.txt -o $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt -r $PBS_O_WORKDIR/$outFolder/$1.refScore.txt > $outFolder/logs/log.oVCFc.$inFile 2> $outFolder/logs/log.oVCFc.$inFile
}

VCFcompare $inFile
#module load R
#Rscript 
