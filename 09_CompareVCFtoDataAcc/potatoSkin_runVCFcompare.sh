
hdf5="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
hdf5_acc="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
kinFile="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.kinship.ibs.hdf5"

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/dth_0.9998_newModScore/the1001genomes_filtered_all_chroms.kinship.ibs.hdf5"


#hdf5="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary.hdf5"
#hdf5_acc="/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary_acc.hdf5"
#kinFile="/lustre/scratch/users/rahul.pisupati/wholeImputed_kinship_ibs_binary_mac5.h5py"


refFol="newFracScores"

declare -a folders=("/lustre/scratch/projects/the1001genomes/rahul/01_C7190ANXX_20150623B_demux_4_rd2_1/cohort/"
"/lustre/scratch/projects/the1001genomes/rahul/02_C7190ANXX_20150623B_demux_5_rd2_2/cohort/"
"/lustre/scratch/projects/the1001genomes/rahul/03_C71EWANXX_20150609B_demux_3_rd3_1/cohort/"
"/lustre/scratch/projects/the1001genomes/rahul/04_C6KD7ANXX_20150327B_demux_2_rd1_2/cohort/"
"/lustre/scratch/projects/the1001genomes/rahul/05_C70WNANXX_20150520B_demux_8_rd4_3/cohort/"
"/lustre/scratch/projects/the1001genomes/rahul/06_C71EWANXX_20150609B_demux_2_rd4_1/cohort/"
"/lustre/scratch/projects/the1001genomes/rahul/07_C719TANXX_20150710B_demux_3_rd1_12/cohort/"
"/lustre/scratch/projects/the1001genomes/rahul/09_C719TANXX_20150710B_demux_5/")

length=${#folders[@]}


for (( i=0; i<$length ;i=i+1 ));do
	cd ${folders[$i]}
	qsub -voutFolder=$refFol,hdf5=$hdf5,hdf5_acc=$hdf5_acc,kinFile=$kinFile ~/MyScripts/CompareVCFtoDataAcc/09_CompareVCFtoDataAcc/qsub_VCFcompare.1.sh
done
