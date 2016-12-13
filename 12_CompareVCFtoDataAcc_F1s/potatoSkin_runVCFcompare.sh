
fol0=`pwd`

folders=(`ls -d [01]*`)


length=${#folders[@]}

for (( i=0; i<$length ;i=i+1 ));do
#	cd $fol0/${folders[$i]}/cohort
#	qsub -N GEAlgo_plate_${i} -P the1001genomes ~/MyScripts/CompareVCFtoDataAcc/12_CompareVCFtoDataAcc_F1s/qsub_VCFcompare.sh
	cd $fol0/${folders[$i]}
	qsub -N GEAlgo_plate_${i} -P field_experiments ~/MyScripts/CompareVCFtoDataAcc/12_CompareVCFtoDataAcc_F1s/qsub_VCFcompare.sh
done
