#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -J 1-96
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=12:00:00
#PBS -o logs.ovcfcompare
#PBS -e logs.evcfcompare

#hdf5="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
#totSNPs="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.numInformative.txt"
#kinFile="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.kinship.ibs.hdf5"
#totSNPs="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.totSNP.count.txt"

hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.hdf5"
hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.acc.hdf5"
kinFile="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.kinship.ibs.hdf5"
totSNPs="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.numInformative.txt"


#kinFile="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.998/the1001genomes_filtered_all_chroms.kinship.ibs.hdf5"
#totSNPs="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.998/the1001genomes_filtered_all_chroms.totalSNPcount.txt"

#-------------Variables
qthres=100
#Path to hdf5 file
#hdf5="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary.hdf5"
#hdf5_acc="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary_acc.hdf5"
#kinFile="/lustre/scratch/users/rahul.pisupati/wholeImputed_kinship_ibs_binary_mac5.h5py"
#totSNPs="/lustre/scratch/users/rahul.pisupati/wholeImputed_totalSNPs_peracc.txt"
#----------------------------

#Get the SNP positions into a file
cd $PBS_O_WORKDIR
outFolder="newCollapsedScores"
mkdir $outFolder
mkdir $outFolder/logs
#inFile=`ls *.filter.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.vcf//'`
inFile=`ls *.filter.gvcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.gvcf//'`
#inFile=`ls *.pos.txt |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.pos\.txt//'`
rm $PBS_O_WORKDIR/$outFolder/$inFile.ScoreAcc.txt $PBS_O_WORKDIR/$outFolder/$inFile.refScore.txt

function VCFcompare {
#	grep -v "^#" $1.filter.vcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	grep -v "^#" $1.filter.gvcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	# Run the python script with matched SNPs
	module load numpy
	module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
	module load pygwas
	module load pandas
	python ~/MyScripts/CompareVCFtoDataAcc/07_CompareVCFtoDataAcc/01_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -t $totSNPs -p $PBS_O_WORKDIR/$outFolder/$1.pos.txt -o $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt -r $PBS_O_WORKDIR/$outFolder/$1.refScore.txt > $outFolder/logs/log.oVCFc.$inFile 2> $outFolder/logs/log.oVCFc.$inFile
}

VCFcompare $inFile
#module load R
#Rscript 
