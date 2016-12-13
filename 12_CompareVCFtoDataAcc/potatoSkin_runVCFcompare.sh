
fol0=`pwd`

folders=(`ls -d [01]*`)


length=${#folders[@]}

for (( i=0; i<$length ;i=i+1 ));do
	cd $fol0/${folders[$i]}/cohort
#	cd $fol0/${folders[$i]}/bsmap/
	numVCF=`ls *.filter.vcf |wc -l`
	echo $numVCF
	qsub -N snpmatch_$((i + 1)) -J 1-$numVCF -P the1001genomes ~/MyScripts/CompareVCFtoDataAcc/12_CompareVCFtoDataAcc/qsub_VCFcompare.sh
#	cd $fol0/${folders[$i]}
#        numVCF=`ls *filter.vcf | wc -l`
#	echo $numVCF
#	qsub -N snpmatch_02_${i} -P field_experiments ~/MyScripts/CompareVCFtoDataAcc/12_CompareVCFtoDataAcc/qsub_VCFcompare.sh
done
