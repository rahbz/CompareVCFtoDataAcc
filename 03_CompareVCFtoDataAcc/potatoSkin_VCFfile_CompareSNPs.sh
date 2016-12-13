#!/bin/bash

inFile=`echo $1 | sed 's/\.vcf$//'`

#-------------Variables
qthres=$2
#Path to hdf5 file
hdf5=/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary.hdf5
hdf5_acc=/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary_acc.hdf5
# Corresponding total number of SNPs generated from the CalculateSNPsearchAcc.py
totSNPs=/lustre/scratch/users/rahul.pisupati/wholeImputed_totalSNPs_peracc.txt

#----------------------------

#Get the SNP positions into a file
paste <(grep -v '^#' $1|awk -v qth=$qthres 'length($4) == 1 && length($5) == 1 && $6 > qth {print $1 "\t" $2}') <(grep -v '^#' $1|awk -v qth=$qthres 'length($4) == 1 && length($5) == 1 && $6 > qth {print $8}'|cut -f2 -d";")|sed 's/^Chr//' > $inFile.pos.txt
#grep -v '#' $1 | awk -v qth=$qthres 'length($4) == 1 && length($5) == 1 && $6 > qth {print $1 "\t" $2}'|sed 's/^Chr//' > $inFile.pos.txt
numSNP=`wc -l $inFile.pos.txt`


# Run the Rscript
module load R
Rscript ~/MyScripts/CompareVCFtoDataAcc/03_CompareVCFtoDataAcc/01_01_FindMatchingSNPpositions.R $hdf5 $inFile.pos.txt $inFile.matchedSNPs

# Run the python script with matched SNPs
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas
python ~/MyScripts/CompareVCFtoDataAcc/03_CompareVCFtoDataAcc/01_02_CalcHammingDist.py -t $totSNPs -m $inFile.matchedSNPs -n $numSNP -o $inFile.ScoreAcc.txt -d $hdf5

# Make refining, generates a .ReqAcc.txt file
Rscript ~/MyScripts/CompareVCFtoDataAcc/03_CompareVCFtoDataAcc/02_01_getMessyAccessionsList.R ./ $inFile
# below python script use the SNP matrix chunked in the columns
if [ -f "$inFile.ReqAcc.txt" ]
then
	python ~/MyScripts/CompareVCFtoDataAcc/03_CompareVCFtoDataAcc/02_02_CalcHammingDist_perAcc.py -p $inFile.pos.txt -r $inFile.ReqAcc.txt -o $inFile.refScoreAcc.txt -d $hdf5_acc
fi




