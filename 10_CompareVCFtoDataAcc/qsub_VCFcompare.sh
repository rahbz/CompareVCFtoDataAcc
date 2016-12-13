#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -J 1-96
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=12:00:00
#PBS -N GEAlgo
#PBS -o logs.ovcfcompare
#PBS -e logs.ovcfcompare

hdf5="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
kinFile="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.kinship.ibs.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_30k/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_30k/the1001genomes_filtered_all_chroms.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_30k/the1001genomes_filtered_all_chroms.kinship.ibs.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.kinship.ibs.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/VCF_1163g/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.hetfiltered.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/VCF_1163g/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.hetfiltered.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/VCF_1163g/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.hetfiltered.csv"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/VQSR_mapping_2015_Fernando/1135g.181k.prior15.ts99.5.BIALLELIC.hetfiltered.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/VQSR_mapping_2015_Fernando/1135g.181k.prior15.ts99.5.BIALLELIC.hetfiltered.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/rahul/VQSR_mapping_2015_Fernando/1135g.181k.prior15.ts99.5.BIALLELIC.hetfiltered.kinship.ibs.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.kinship.ibs.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.995_newFracScore/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.995_newFracScore/the1001genomes_filtered_all_chroms.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.995_newFracScore/the1001genomes_filtered_all_chroms.kinship.ibs.hdf5"

#-------------Variables
qthres=50
#Path to hdf5 file
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary_acc.hdf5"
#outFolder="newLikelihoodScore_wholeImputed"
#kinFile="/lustre/scratch/users/rahul.pisupati/wholeImputed_kinship_ibs_binary_mac5.h5py"

#hdf5="/lustre/scratch/users/rahul.pisupati/250kdata_all_chromosomes_binary.hdf5"
#hdf5_acc="/lustre/scratch/users/rahul.pisupati/250kdata_all_chromosomes_binary.acc.hdf5"
#kinFile=""

#----------------------------

#Get the SNP positions into a file
cd $PBS_O_WORKDIR
#outFolder="newRefineFracScore_newPval"
outFolder="newLikelihoodScore_1135SNP_1"
#outFolder="newLikelihoodScore_collapsed5k"
mkdir $outFolder
mkdir $outFolder/logs

#inFile=`ls *.genotyped.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.genotyped\.vcf//'`
inFile=`ls *.filter.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.vcf//'`
#inFile=`ls *.filter.gvcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.gvcf//'`
#inFile=`ls *.pos.txt |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.pos\.txt//'`
rm $PBS_O_WORKDIR/$outFolder/$inFile.ScoreAcc.txt $PBS_O_WORKDIR/$outFolder/$inFile.refScore.txt

function VCFcompare {
#	grep -v "^#" $1.genotyped.vcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	grep -v "^#" $1.filter.vcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
#	grep -v "^#" $1.filter.gvcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	module load numpy/1.9.1-goolf-1.4.10-Python-2.7.3
	module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
	module load pygwas
	module load pandas/0.16.2-goolf-1.4.10-Python-2.7.3
	module load ipython
	python ~/MyScripts/CompareVCFtoDataAcc/10_CompareVCFtoDataAcc/01_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -p $PBS_O_WORKDIR/$outFolder/$1.pos.txt -o $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt -r $PBS_O_WORKDIR/$outFolder/$1.refScore.txt > $outFolder/logs/log.oVCFc.$inFile 2> $outFolder/logs/log.oVCFc.$inFile
	#python ~/MyScripts/CompareVCFtoDataAcc/10_CompareVCFtoDataAcc/01_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -p $PBS_O_WORKDIR/$1.pos.txt -o $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt -r $PBS_O_WORKDIR/$outFolder/$1.refScore.txt > $outFolder/logs/log.oVCFc.$inFile 2> $outFolder/logs/log.oVCFc.$inFile
	# Get the matched and mismatch positions for the TopHit
	TopAcc=`sort -k5,5g $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt | head -n 1 | cut -f1`
	python ~/MyScripts/CompareVCFtoDataAcc/10_CompareVCFtoDataAcc/03_getMatchdPos.py -d $hdf5 -e $hdf5_acc -p $PBS_O_WORKDIR/$outFolder/$1.pos.txt -m $PBS_O_WORKDIR/$outFolder/$1.matchedSNPs -n $PBS_O_WORKDIR/$outFolder/$1.nmatchedSNPs -a $TopAcc > $outFolder/logs/log.matSNP.$inFile 2> $outFolder/logs/log.matSNP.$inFile
}

VCFcompare $inFile
#module load R
#Rscript 
