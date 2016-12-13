
## Load modules R
module load R


#refFol="newLikelihoodScore_1135SNP"
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
#ecos="/lustre/scratch/projects/the1001genomes/rahul/ecotypesids_merged_added.csv"

#refFol="newLikelihoodScore_collapsed5k"
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.hdf5"
#ecos="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.ecotypesids_merged.csv"

refFol="newLikelihoodScore_wholeImputed"
hdf5="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary.hdf5"
ecos="/lustre/scratch/projects/the1001genomes/rahul/ecotypesids_merged_added.csv"


folders=(`ls -d /lustre/scratch/projects/the1001genomes/rahul/[01]*`)
ids=(`echo "rd2_1,rd2_2,rd3_1,rd1_2,rd4_3,rd4_1,rd1_12,rd2_8,rd2_9,rd2_11,rd2_12,rd2_6" | sed 's/,/\n/g'`)


script="~/MyScripts/CompareVCFtoDataAcc/10_CompareVCFtoDataAcc/02_makeCSVTable_CompareAccessions.R"
output="intermediate_modified.csv"

length=${#folders[@]}
for (( i=0; i<$length ;i=i+1 ));do
	acc=`ls ${folders[$i]}/acc*`
	cd ${folders[$i]}/cohort/$refFol
#	Rscript ~/MyScripts/CompareVCFtoDataAcc/10_CompareVCFtoDataAcc/02_makeCSVTable_CompareAccessions.R $ecos $output $acc ${ids[$i]}
done

IntFiles=( "${folders[@]/%//cohort/${refFol}/${output}}" )
(head -n 1 ${IntFiles[0]} && tail -n +2 ${IntFiles[@]} ) | grep -v "^=" | sed '/^\s*$/d' 
