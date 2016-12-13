#!/bin/bash

#inFile=`echo $1 | sed 's/\.vcf$//'`
inFile=$1

#-------------Variables
qthres=$2
#Path to hdf5 file
hdf5=$3
hdf5_acc=$4
# Corresponding total number of SNPs generated from the CalculateSNPsearchAcc.py
totSNPs=$5

#----------------------------

#Get the SNP positions into a file
#grep -v "^#" $1 | awk -v qth=$qthres '$6 > qth && $NF ~ /^1\/1|^0\/1|^1\/0/ {print $1 "\t" $2 "\t" $6 "\t" $(NF-1) "\t" $NF}' | sed 's/^Chr//' > $1.pos.txt

numSNP=`wc -l $inFile.pos.txt`


# Run the Rscript
module load R
Rscript ~/MyScripts/04_CompareVCFtoDataAcc/01_01_FindMatchingSNPpositions.R $hdf5 $inFile.pos.txt $inFile.matchedSNPs 

# Run the python script with matched SNPs
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas
python ~/MyScripts/04_CompareVCFtoDataAcc/01_02_CalcHammingDist.py -d $hdf5 -e $hdf5_acc -t $totSNPs -m $1.matchedSNPs -n $numSNP -o $1.ScoreAcc.txt -r $1.refScore.txt 


