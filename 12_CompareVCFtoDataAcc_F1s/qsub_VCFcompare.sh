#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -J 1-5
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=24:00:00
#PBS -N GENAlgo
#PBS -o logs.ovcfcompare
#PBS -e logs.ovcfcompare

hdf5="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
outFolder="newFractionScore_F1s_1135SNP"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_30k/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_30k/the1001genomes_filtered_all_chroms.acc.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.acc.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/VCF_1163g/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.hetfiltered.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/VCF_1163g/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.hetfiltered.acc.hdf5"

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
#outFolder="newFractionScore_F1s_swedes243"

#-------------Variables
#Path to hdf5 file
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary_acc.hdf5"
#outFolder="newFractionScore_F1s_wholeImputed"


#hdf5="/lustre/scratch/projects/the1001genomes/rahul/VCF_250k_nonImputed_callmethod74_1Dec2015/call_method_74.binary.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/VCF_250k_nonImputed_callmethod74_1Dec2015/call_method_74.binary.acc.hdf5"
#outFolder="newFractionScore_F1s_250k_unImputed"

#----------------------------

#Get the SNP positions into a file
cd $PBS_O_WORKDIR
mkdir $outFolder
mkdir $outFolder/logs

#inFile=`ls *.genotyped.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.genotyped\.vcf//'`
inFile=`ls *.filter.vcf |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.filter\.vcf//'`
rm $PBS_O_WORKDIR/$outFolder/$inFile.ScoreAcc.txt $PBS_O_WORKDIR/$outFolder/$inFile.refScore.txt

function VCFcompare {
	module load numpy/1.9.1-goolf-1.4.10-Python-2.7.3
	module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
	module load pygwas
	module load vcfnp
	module load pandas/0.16.2-goolf-1.4.10-Python-2.7.3
	module load ipython
	python ~/MyScripts/CompareVCFtoDataAcc/12_CompareVCFtoDataAcc_F1s/01_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -i $PBS_O_WORKDIR/$1.filter.vcf -o $PBS_O_WORKDIR/$outFolder/$1.ScoreAcc.txt -r $PBS_O_WORKDIR/$outFolder/$1.refScore.txt > $outFolder/logs/log.oVCFc.$inFile 2> $outFolder/logs/log.oVCFc.$inFile
}

VCFcompare $inFile
