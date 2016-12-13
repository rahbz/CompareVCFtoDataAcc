#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -J 0-95
#PBS -l walltime=12:00:00
#PBS -o logs.ovcfcompare
#PBS -e logs.evcfcompare
#Input:  VCF files generated after SNP calling
#Output: Compare accessions


cd $PBS_O_WORKDIR

qthres=10000

hdf5="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary.hdf5"
hdf5_acc="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary_acc.hdf5"
# Corresponding total number of SNPs generated from the CalculateSNPsearchAcc.py
totSNPs="/lustre/scratch/users/rahul.pisupati/wholeImputed_totalSNPs_peracc.txt"

#inFiles=(`ls *.filter.gvcf | sed 's/\.filter\.gvcf//'`)
inFiles=(`ls *.pos.txt | sed 's/\.pos\.txt//'`)

function VCFcompare {
	inFile=$1
#	grep -v "^#" $1.filter.gvcf | awk -v qth=$2 '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0/ {print $1 "\t" $2 "\t" $6 "\t" $(NF-1) "\t" $NF}' | sed 's/^Chr//' > $1.pos.txt
	awk -v qth=$2 '$3 > qth {print $0}' $inFile.pos.txt  > $inFile.pos1.txt
	numSNP=`wc -l $inFile.pos1.txt`
	# Run the Rscript
	module load R
	Rscript ~/MyScripts/04_CompareVCFtoDataAcc/01_01_FindMatchingSNPpositions.R $hdf5 $inFile.pos1.txt $inFile.matchedSNPs 
	# Run the python script with matched SNPs
	module load numpy
	module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
	module load pygwas
	python ~/MyScripts/04_CompareVCFtoDataAcc/01_02_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -t $totSNPs -m $inFile.matchedSNPs -n $numSNP -o $inFile.ScoreAcc.txt -r $inFile.refScore.txt 
}

VCFcompare ${inFiles[$PBS_ARRAY_INDEX]} $qthres
