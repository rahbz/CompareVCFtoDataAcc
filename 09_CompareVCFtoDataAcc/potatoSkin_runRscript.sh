

refFol="newFracScores"

declare -a folders=("/lustre/scratch/projects/the1001genomes/rahul/01_C7190ANXX_20150623B_demux_4_rd2_1/cohort/$refFol/"
"/lustre/scratch/projects/the1001genomes/rahul/02_C7190ANXX_20150623B_demux_5_rd2_2/cohort/$refFol/"
"/lustre/scratch/projects/the1001genomes/rahul/03_C71EWANXX_20150609B_demux_3_rd3_1/cohort/$refFol/"
"/lustre/scratch/projects/the1001genomes/rahul/04_C6KD7ANXX_20150327B_demux_2_rd1_2/cohort/$refFol/"
"/lustre/scratch/projects/the1001genomes/rahul/05_C70WNANXX_20150520B_demux_8_rd4_3/cohort/$refFol/"
"/lustre/scratch/projects/the1001genomes/rahul/06_C71EWANXX_20150609B_demux_2_rd4_1/cohort/$refFol/"
"/lustre/scratch/projects/the1001genomes/rahul/07_C719TANXX_20150710B_demux_3_rd1_12/cohort/$refFol/")

length=${#folders[@]}


for (( i=0; i<$length ;i=i+1 ));do
	cd ${folders[$i]}
	acc=`ls ../../`
	qsub ~/MyScripts/CompareVCFtoDataAcc/08_CompareVCFtoDataAcc/qsub_VCFcompare.sh
done
