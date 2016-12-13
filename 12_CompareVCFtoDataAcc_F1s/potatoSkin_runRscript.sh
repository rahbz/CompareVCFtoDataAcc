
## Load modules R
module load R/3.1.0-goolf-1.4.10

#refFol="finalLikelihoodScore_1135SNP"
hdf5="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
ecos="/lustre/scratch/projects/the1001genomes/rahul/ecotypesids_merged_added.csv"

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


folders=(`ls -d /lustre/scratch/projects/the1001genomes/rahul/[01]*`)
ids=(`echo "rd2_1,rd2_2,rd3_1,rd1_2,rd4_3,rd4_1,rd1_12,rd2_8,rd2_9,rd2_11,rd2_12,rd2_6" | sed 's/,/\n/g'`)

fol0=`pwd`
#folders=(`ls -d /lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/003.plate.genotyping/0*`)
#hdf5="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/998.Swedes.243.SNPmatrix/Swedes243.177k.prior15.gauss6.ts99.5.BIALLELIC.hdf5"
#ids=(`echo "C6D6,C6KD,C6K7,C5E7,C6KY,C5M7,C6K6,C6K8,C6D5,C5M6,C5M8,C6KF,C6D4,C7V1,C7G5,C7G6" | sed 's/,/\n/g'`)
#refFol="newFractionScore_F1s_swedes243"
refFol="newFractionScore_F1s_1135SNP"

#fol0=`pwd`
#folders=(`ls -d /lustre/scratch/projects/the1001genomes/rahul/BiSulphide_Manus_Seqs/0*`)
#ids=(`echo "Lib14a,Lib14b,Lib14c,Lib15a,Lib15b,6,7,8,9,10,11,12" | sed 's/,/\n/g'`)
#refFol="finalLikelihoodScore_1135SNP"

#fol0=`pwd`
#folders=(`ls -d /lustre/scratch/projects/the1001genomes/rahul/Data_250k/0*`)
#ids=(`echo "250k_n3,250k_n4,250k_n5,250k_n6,250k_n2,250k_n9" | sed 's/,/\n/g'`)
#refFol="finalLikelihoodScore_250k_unImputed"


script="~/MyScripts/CompareVCFtoDataAcc/11_CompareVCFtoDataAcc/02_makeCSVTable_CompareAccessions.R"
output="intermediate_modified.csv"

length=${#folders[@]}
for (( i=0; i<$length ;i=i+1 ));do
#	acc=`ls ${folders[$i]}/acc*`
	cd ${folders[$i]}/cohort/$refFol
#	cd ${folders[$i]}/bsmap/$refFol
#	cd ${folders[$i]}/$refFol
#	Rscript ~/MyScripts/CompareVCFtoDataAcc/12_CompareVCFtoDataAcc_F1s/05_makeCSVTable_fieldExp.R $hdf5 $output ${ids[$i]}
#	Rscript ~/MyScripts/CompareVCFtoDataAcc/12_CompareVCFtoDataAcc_F1s/02_makeCSVTable_CompareAccessions.R $ecos $output $acc ${ids[$i]}
done


#IntFiles=( "${folders[@]/%//bsmap/${refFol}/${output}}" )
IntFiles=( "${folders[@]/%//cohort/${refFol}/${output}}" )
(head -n 1 ${IntFiles[0]} && tail -n +2 ${IntFiles[@]} ) | grep -v "^=" | sed '/^\s*$/d' 

#IntFiles=( "${folders[@]/%//${refFol}/${output}}" )
#(head -n 1 ${IntFiles[0]} && tail -n +2 ${IntFiles[@]} ) | grep -v "^=" | sed '/^\s*$/d' 
