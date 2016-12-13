
workingDir=``
cd $workingDir

#-------------Variables
qthres=100
#Path to hdf5 file
hdf5="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary.hdf5"
hdf5_acc="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary_acc.hdf5"
# Corresponding total number of SNPs generated from the CalculateSNPsearchAcc.py
totSNPs="/lustre/scratch/users/rahul.pisupati/wholeImputed_totalSNPs_peracc.txt"
pathmerged="/lustre/scratch/projects/the1001genomes/rahul/ecotypesids_merged_filtered.csv"
outCSV="intermediate_modified.csv"
#----------------------------

qsub ~/MyScripts/04_CompareVCFtoDataAcc/qsub_VCFcompare.sh

# Make CSV file
#module load R
#id=`ls *.filter.gvcf |head -n 1 | sed 's/\#.*$/\#/'`
#acc=`ls acc_*`
#Rscript ~/MyScripts/04_CompareVCFtoDataAcc/02_makeCSVTable_CompareAccessions.R $hdf5 $id $acc $pathmerged $outCSV
