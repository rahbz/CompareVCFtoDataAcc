## BASH script to get the het calls from the POS files


posFiles=(`ls *pos.txt | sed 's/\.pos\.txt$//'`)
num=${#posFiles[@]}

module load Perl

for (( i=0; i<$num ;i=i+1 ));do
	hets=(`cut -f3 ${posFiles[$i]}.pos.txt | sort | uniq -c | awk '{print $1}'`)
	echo "${posFiles[$i]},${hets[0]},${hets[1]},${hets[2]}"
done
