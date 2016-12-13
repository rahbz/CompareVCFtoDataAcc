#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -J 1-96
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=12:00:00
#PBS -o logs.ovcfcompare
#PBS -e logs.evcfcompare

#hdf5="/lustre/scratch/users/rahul.pisupati/simulateSNPs/wholeImputed_all_chromosomes_filtered.hdf5"
#hdf5_acc="/lustre/scratch/users/rahul.pisupati/simulateSNPs/wholeImputed_all_chromosomes_filtered_acc.hdf5"
#totSNPs="/lustre/scratch/users/rahul.pisupati/simulateSNPs/wholeImputed_totalSNPs_peracc.txt"
#kinFile="/lustre/scratch/users/rahul.pisupati/simulateSNPs/wholeImputed_kinship_filtered.hdf5"

#-------------Variables
qthres=100
#Path to hdf5 file
hdf5="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary.hdf5"
hdf5_acc="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary_acc.hdf5"

# Corresponding total number of SNPs generated from the CalculateSNPsearchAcc.py
totSNPs="/lustre/scratch/users/rahul.pisupati/wholeImputed_totalSNPs_peracc.txt"
#----------------------------

#Get the SNP positions into a file
cd $PBS_O_WORKDIR
mkdir logs
inFile=`ls *.filter.gvcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.gvcf//'`
#inFile=`ls *.pos.txt |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.pos\.txt//'`

function VCFcompare {
	grep -v "^#" $1.filter.gvcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0/ {print $1 "\t" $2 "\t" $(NF-1) "\t" $NF}' | sed 's/^Chr//' > $1.pos.txt
	numSNP=`wc -l $1.pos.txt`
	# Run the Rscript
	module load R
	Rscript ~/MyScripts/04_CompareVCFtoDataAcc/01_01_FindMatchingSNPpositions.R $hdf5 $1.pos.txt $1.matchedSNPs > logs/log.oVCFc.${1} 2> logs/log.oVCFc.${1}
	# Run the python script with matched SNPs
	module load numpy
	module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
	module load pygwas
	python ~/MyScripts/04_CompareVCFtoDataAcc/01_02_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -t $totSNPs -m $1.matchedSNPs -n $numSNP -o $1.ScoreAcc.txt -r $1.refScore.txt > logs/log.oVCFc.$i 2> logs/log.oVCFc.$i
}

VCFcompare $inFile
