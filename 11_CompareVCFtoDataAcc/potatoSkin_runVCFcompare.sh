
fol0=`pwd`

folders=(`ls -d 0*`)


length=${#folders[@]}

for (( i=0; i<$length ;i=i+1 ));do
#	cd $fol0/${folders[$i]}/cohort
	cd $fol0/${folders[$i]}/bsmap/
	qsub -N GEAlgo_plate_${i} -P the1001genomes ~/MyScripts/CompareVCFtoDataAcc/11_CompareVCFtoDataAcc/qsub_VCFcompare.sh
#	cd $fol0/${folders[$i]}
#	qsub -N GEAlgo_plate_${i} -P field_experiments ~/MyScripts/CompareVCFtoDataAcc/11_CompareVCFtoDataAcc/qsub_VCFcompare.sh
done
