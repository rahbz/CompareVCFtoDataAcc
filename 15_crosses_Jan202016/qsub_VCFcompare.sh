#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -J 1-48
#PBS -N GENCross
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=12:00:00
#PBS -o logs.ovcfcompare
#PBS -e logs.ovcfcompare

hdf5="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
outFolder="newLikelihoodMatrix_1135SNP_300Kb"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_30k/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_30k/the1001genomes_filtered_all_chroms.acc.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.acc.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/VCF_1163g_Jan20/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.BINARY.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/VCF_1163g_Jan20/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.BINARY.acc.hdf5"
#outFolder="newLikelihoodMatrix_1163g_300kb"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/VQSR_mapping_2015_Fernando/1135g.181k.prior15.ts99.5.BIALLELIC.hetfiltered.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/VQSR_mapping_2015_Fernando/1135g.181k.prior15.ts99.5.BIALLELIC.hetfiltered.acc.hdf5"

#hdf5="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/998.Swedes.243.SNPmatrix/Swedes243.177k.prior15.gauss6.ts99.5.BIALLELIC.hdf5"
#hdf5_acc="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/998.Swedes.243.SNPmatrix/Swedes243.177k.prior15.gauss6.ts99.5.BIALLELIC.acc.hdf5"
#outFolder="newLikelihoodMatrix_swedes243_300kb"


#-------------Variables
qthres=50
#Path to hdf5 file
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary_acc.hdf5"
#outFolder="newLikelihoodMatrix_wholeImputed_300kb"

#hdf5="/lustre/scratch/users/rahul.pisupati/250kdata_all_chromosomes_binary.hdf5"
#hdf5_acc="/lustre/scratch/users/rahul.pisupati/250kdata_all_chromosomes_binary.acc.hdf5"
#kinFile=""

#----------------------------

#Get the SNP positions into a file
cd $PBS_O_WORKDIR
bin=300000
mkdir $outFolder
mkdir $outFolder/logs
verFol="~/MyScripts/CompareVCF_Crosses/02_Jan202016"

#inFile=`ls *.genotyped.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.genotyped\.vcf//'`
inFile=`ls *.filter.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.vcf//'`
#inFile=`ls *.filter.gvcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.gvcf//'`
#inFile=`ls *.pos.txt |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.pos\.txt//'`
rm $PBS_O_WORKDIR/$outFolder/$inFile.matLikeliAcc.txt

function VCFcompare {
#	grep -v "^#" $1.filter.vcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	grep -v "^#" $1.filter.vcf | awk -v qth=$qthres '$7 ~ /PASS/ && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//i' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	module load numpy/1.9.1-goolf-1.4.10-Python-2.7.3
	module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
	module load pygwas
	module load pandas/0.16.2-goolf-1.4.10-Python-2.7.3
	module load ipython
	python ~/MyScripts/CompareVCF_Crosses/02_Jan202016/01_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -p $PBS_O_WORKDIR/$outFolder/$1.pos.txt -o $PBS_O_WORKDIR/$outFolder/$1.matLikeliAcc.txt -b $bin -s $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt > $outFolder/logs/log.oVCFc.$inFile 2> $outFolder/logs/log.oVCFc.$inFile
}

VCFcompare $inFile
