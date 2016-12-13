#!/bin/bash

#Not a qsub script, have to call from interactive shell

csvFile=$1

# Run Rscript to get all the AccList

# Run the python script with matched SNPs
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas

python ~/MyScripts/03_CompareVCFtoDataAcc/CalcHammingDist.py -t /lustre/scratch/users/rahul.pisupati/totalSNPsNum_1001genomes.txt -p $inFile.pos.txt -o $inFile.ScoreAcc.txt

