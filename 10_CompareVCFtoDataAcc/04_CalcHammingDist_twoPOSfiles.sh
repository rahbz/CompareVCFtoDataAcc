#!/usr/bin


# Input two pos files


tempFol="/lustre/scratch/projects/the1001genomes/rahul/TempFiles/"

module load Perl

info=`(awk '{print $1 "," $2}' $1 && awk '{print $1 "," $2}' $2)|sort -T $tempFol | uniq -d |wc -l`
mat=`(awk '{print $1 "," $2 "," $3}' $1 && awk '{print $1 "," $2 "," $3}' $2)|sort -T $tempFol | uniq -d |wc -l`
sc=`perl -e 'print $ARGV[0]/$ARGV[1]' $mat $info`



echo "$1	$2	$mat	$info	$sc"

