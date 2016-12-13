

folders=(`ls -d /lustre/scratch/projects/the1001genomes/rahul/[01]*`)


length=${#folders[@]}

for (( i=0; i<$length ;i=i+1 ));do
	cd ${folders[$i]}/cohort
	qsub ~/MyScripts/CompareVCFtoDataAcc/10_CompareVCFtoDataAcc/qsub_VCFcompare.sh -N GEAlgo_plate_${i}
done
