#!/bin/bash

#Not a qsub script, have to call from interactive shell

#Input the VCF file from the command line
inFile=`echo $1 | sed 's/\.vcf$//'`

#Quality threshold that has to be used
qthres=$2

#Get the SNP positions into a file
grep -v '#' $1 | awk -v qth=$qthres 'length($4) == 1 && length($5) == 1 && $6 > qth {print $1 "\t" $2}'|sed 's/^Chr//' > $inFile.pos.txt
numSNP=`wc -l $inFile.pos.txt`

# Run the python script with matched SNPs
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas

#HDF5 file for the binaries
hdf5=/lustre/scratch/users/rahul.pisupati/all_chromosomes_binary.hdf5

python ~/MyScripts/04_CompareVCFtoDataAcc/00_CalcHammingDist.py -t /lustre/scratch/users/rahul.pisupati/totalSNPsNum_1001genomes.txt -p $inFile.pos.txt -o $inFile.ScoreAcc.txt -d $hdf5

