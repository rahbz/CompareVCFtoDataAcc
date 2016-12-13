
## Load modules R
module load R

#refFol="finalLikelihoodScore_1135SNP"
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
#ecos="/lustre/scratch/projects/the1001genomes/rahul/ecotypesids_merged_added.csv"

#refFol="newLikelihoodScore_1135SNP"
#refFol="finalLikelihoodScore_1135SNP"
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
#ecos="/lustre/scratch/projects/the1001genomes/rahul/ecotypesids_merged_added.csv"

#refFol="newLikelihoodScore_collapsed5k"
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.hdf5"
#ecos="/lustre/scratch/projects/the1001genomes/rahul/simulateSNPs_the1001genomes_unimputed/num_5k/the1001genomes_filtered_all_chroms.ecotypesids_merged.csv"

#refFol="finalLikelihoodScore_wholeImputed"
#hdf5="/lustre/scratch/projects/the1001genomes/rahul/wholeImputed_all_chromosomes_binary.hdf5"
#ecos="/lustre/scratch/projects/the1001genomes/rahul/ecotypesids_merged_added.csv"

#fol0=`pwd`
#folders=(`ls -d /lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/001.plate.raw.vcfs/0*`)
#ids=(`echo "C6D6,C6KD,C6K7,C5E7,C6KY,C5M7,C6K6,C6K8,C6D5,C5M6,C5M8,C6KF,C6D4,C7V1,C7G5,C7G6" | sed 's/,/\n/g'`)
#refFol="newLikelihoodMatrix_swedes243_300kb"


#fol0=`pwd`
#folders=(`ls -d /lustre/scratch/projects/the1001genomes/rahul/fromEnvel/[01]*`)
#refFol="newLikelihoodMatrix_wholeImputed_300kb"
#refFol="newLikelihoodMatrix_1135SNP_300Kb"


fol0=`pwd`
folders=(`ls -d /lustre/scratch/projects/the1001genomes/rahul/BiSulphide_Manus_Seqs/0*`)
ids=(`echo "Lib14a,Lib14b,Lib14c,Lib15a,Lib15b,6,7,8,9,10,11,12" | sed 's/,/\n/g'`)
refFol="newLikelihoodMatrix_1135SNP_300Kb"


script="~/MyScripts/CompareVCF_Crosses/01_Jan092016/02_makeCSVTable_Fernando_readingF2s.R"
output="intermediate_modified.csv"

length=${#folders[@]}
for (( i=0; i<$length ;i=i+1 ));do
#	acc=`ls ${folders[$i]}/acc*`
#	cd ${folders[$i]}/cohort/$refFol
	cd ${folders[$i]}/bsmap/$refFol
#	cd ${folders[$i]}/$refFol
#	Rscript ~/MyScripts/CompareVCF_Crosses/02_Jan202016/02_makeCSVTable_Fernando_readingF2s.R $output $(($i + 1))
#	Rscript ~/MyScripts/CompareVCF_Crosses/02_Jan202016/02_makeCSVTable_Fernando_readingF2s.R $output ${ids[$i]}
done


#IntFiles=( "${folders[@]/%//cohort/${refFol}/${output}}" )
IntFiles=( "${folders[@]/%//bsmap/${refFol}/${output}}" )
(head -n 1 ${IntFiles[0]} && tail -n +2 ${IntFiles[@]} ) | grep -v "^=" | sed '/^\s*$/d' 

#IntFiles=( "${folders[@]/%//${refFol}/${output}}" )
#(head -n 1 ${IntFiles[0]} && tail -n +2 ${IntFiles[@]} ) | grep -v "^=" | sed '/^\s*$/d' 
