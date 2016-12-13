#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -J 1-98
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -l walltime=24:00:00
#PBS -N snpmatch
#PBS -o logs.ovcfcompare
#PBS -e logs.ovcfcompare

hdf5="/lustre/scratch/projects/the1001genomes/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
outFolder="finalLikelihoodScore_1135SNP"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_30k/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_30k/the1001genomes_filtered_all_chroms.acc.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.acc.hdf5"

#f5="/lustre/scratch/projects/the1001genomes/VCF_1163g_Jan20/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.hetfiltered.hdf5"
#df5_acc="/lustre/scratch/projects/the1001genomes/VCF_1163g_Jan20/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.hetfiltered.acc.hdf5"
#utFolder="finalLikelihoodScore_1163SNP"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/VQSR_mapping_2015_Fernando/1135g.181k.prior15.ts99.5.BIALLELIC.hetfiltered.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/VQSR_mapping_2015_Fernando/1135g.181k.prior15.ts99.5.BIALLELIC.hetfiltered.acc.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.acc.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.995_newFracScore/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.995_newFracScore/the1001genomes_filtered_all_chroms.acc.hdf5"

#hdf5="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/999.Swedes.198.SNPmatrix/Swede.198.hetfiltered.bialleleic.22Jun15.hdf5"
#hdf5_acc="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/999.Swedes.198.SNPmatrix/Swede.198.hetfiltered.bialleleic.22Jun15.acc.hdf5"
#outFolder="finalLikelihoodScore_swedes198"

#hdf5="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/998.Swedes.243.SNPmatrix/Swedes243.177k.prior15.gauss6.ts99.5.BIALLELIC.hdf5"
#hdf5_acc="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/998.Swedes.243.SNPmatrix/Swedes243.177k.prior15.gauss6.ts99.5.BIALLELIC.acc.hdf5"
#outFolder="finalLikelihoodScore_swedes243"

#hdf5="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/997.Swedes.220.10May2016/01_2.7M_Swedes220.175k.prior15.gauss4.ts99.8.BIALLELIC.hdf5"
#hdf5_acc="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/997.Swedes.220.10May2016/01_2.7M_Swedes220.175k.prior15.gauss4.ts99.8.BIALLELIC.acc.hdf5"
#outFolder="finalLikelihoodScore_01.swedes220"

#hdf5="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/997.Swedes.220.10May2016/02_2.3M_Swedes220.175k.prior15.gauss4.ts99.5.BIALLELIC.hdf5"
#hdf5_acc="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/997.Swedes.220.10May2016/02_2.3M_Swedes220.175k.prior15.gauss4.ts99.5.BIALLELIC.acc.hdf5"
#outFolder="finalLikelihoodScore_02.swedes220"


#-------------Variables
#Path to hdf5 file
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary_acc.hdf5"
#outFolder="finalLikelihoodScore_wholeImputed"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/VCF_250k_nonImputed_callmethod74_1Dec2015/call_method_74.binary.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/VCF_250k_nonImputed_callmethod74_1Dec2015/call_method_74.binary.acc.hdf5"
#outFolder="finalLikelihoodScore_250k_unImputed"
#----------------------------

## Aquilegia, gokce files
#hdf5="/lustre/scratch/projects/aquilegia/998.rahul_snpmatch/all.49inds.plus.semiaq.filter.multi.indel.repeat.BIALLELIC.hdf5"
#hdf5_acc="/lustre/scratch/projects/aquilegia/998.rahul_snpmatch/all.49inds.plus.semiaq.filter.multi.indel.repeat.BIALLELIC.acc.hdf5"
#outFolder="finalLikelihoodScore_50Acc"


#Get the SNP positions into a file
cd $PBS_O_WORKDIR
mkdir $outFolder
mkdir $outFolder/logs

qthres=50

#inFile=`ls *.genotyped.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.genotyped\.vcf//'`
inFile=`ls *.filter.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.vcf//'`
#inFile=`ls *.filter.gvcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.gvcf//'`
#inFile=`ls *.pos.txt |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.pos\.txt//'`
rm $PBS_O_WORKDIR/$outFolder/$inFile.ScoreAcc.txt $PBS_O_WORKDIR/$outFolder/$inFile.refScore.txt

function VCFcompare {
#	grep -v "^#" $1.genotyped.vcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//i' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	grep -v "^#" $1.filter.vcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//i' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
#	grep -v "^#" $1.filter.gvcf | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^Chr//i' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
#	grep -v "^#" $1.filter.vcf | awk '$7 == "PASS" && $NF ~ /^1\/1|^0\/1|^1\/0|^0\/0/ {print $1 "\t" $2 "\t" substr($NF,1,3)}' | sed 's/^chr//i' > $PBS_O_WORKDIR/$outFolder/$1.pos.txt
	module load numpy/1.9.1-goolf-1.4.10-Python-2.7.3
	module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
	module load pygwas
	module load pandas/0.16.2-goolf-1.4.10-Python-2.7.3
	module load ipython
	python ~/MyScripts/CompareVCFtoDataAcc/11_CompareVCFtoDataAcc/01_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -p $PBS_O_WORKDIR/$outFolder/$1.pos.txt -o $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt > $outFolder/logs/log.oVCFc.$inFile 2> $outFolder/logs/log.oVCFc.$inFile
#	python ~/MyScripts/CompareVCFtoDataAcc/11_CompareVCFtoDataAcc/01_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -p $PBS_O_WORKDIR/$1.pos.txt -o $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt > $outFolder/logs/log.oVCFc.$inFile 2> $outFolder/logs/log.oVCFc.$inFile
	# Get the matched and mismatch positions for the TopHit
#	TopAcc=`sort -k5,5g $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt | head -n 1 | cut -f1`
#	python ~/MyScripts/CompareVCFtoDataAcc/11_CompareVCFtoDataAcc/03_getMatchdPos.py -d $hdf5 -e $hdf5_acc -p $PBS_O_WORKDIR/$outFolder/$1.pos.txt -m $PBS_O_WORKDIR/$outFolder/$1.matchedSNPs -n $PBS_O_WORKDIR/$outFolder/$1.nmatchedSNPs -a $TopAcc > $outFolder/logs/log.matSNP.$inFile 2> $outFolder/logs/log.matSNP.$inFile
}

VCFcompare $inFile
#module load R
#Rscript 
