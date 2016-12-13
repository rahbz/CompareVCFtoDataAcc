#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N matchedSNPs
#PBS -J 1-22
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=6:00:00
#PBS -o logs.ogetSNP
#PBS -e logs.egetSNP

#hdf5="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.kinship.ibs.hdf5"


hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.hdf5"
hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.acc.hdf5"
kinFile="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.kinship.ibs.hdf5"

cd $PBS_O_WORKDIR
inFile=`ls *.pos.txt |head -n $PBS_ARRAY_INDEX | tail -n 1 | sed 's/\.pos\.txt//'`
#acc=6974
#acc=9970
acc=6909
#acc=5253
mkdir logs

module load statsmodels
module load numpy/1.9.1-goolf-1.4.10-Python-2.7.3
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas
module load pandas/0.16.2-goolf-1.4.10-Python-2.7.3
module load ipython

python ~/MyScripts/CompareVCFtoDataAcc/09_CompareVCFtoDataAcc/03_getMatchdPos.py -d $hdf5 -e $hdf5_acc -a $acc -p $inFile.pos.txt -m $inFile.matchedSNPs.$acc -n $inFile.nmatchedSNPs.$acc > logs/log.$inFile 2> logs/log.$inFile


