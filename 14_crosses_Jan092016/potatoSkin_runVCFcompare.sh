
fol0=`pwd`

folders=(`ls -d 0*`)


length=${#folders[@]}

for (( i=0; i<$length ;i=i+1 ));do
#	cd ${folders[$i]}/cohort
	cd $fol0/${folders[$i]}/cohort/
#	echo "$fol0/${folders[$i]}"
#	qsub -N GENcross_p.${i} -P field_experiments ~/MyScripts/CompareVCF_Crosses/01_Jan092016/qsub_VCFcompare.sh
	qsub -N GENcross_p.${i} -P the1001genomes ~/MyScripts/CompareVCF_Crosses/01_Jan092016/qsub_VCFcompare.sh
done
