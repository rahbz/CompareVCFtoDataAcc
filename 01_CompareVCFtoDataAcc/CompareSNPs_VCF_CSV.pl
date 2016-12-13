#!/usr/bin/perl
use warnings;
use strict;

# Check the given SNPs in the CSV table

my $inVCFfile= $ARGV[1]; chomp $inVCFfile;
my $dataSNPs = $ARGV[0]; chomp $dataSNPs;
# Modules needed to run the python script
my $headDataSNPs = `head -n 1 $dataSNPs`;chomp $headDataSNPs;
my @allPopIDs = split(',', $headDataSNPs);
# Create a bash generated search pattern file
# Create a hash for the all the populations in the DATA SNPs
my %ScoreList ;
for (my $i = 2; $i < @allPopIDs ; $i++){
	$ScoreList{$allPopIDs[$i]} = 0;
}

open (inputSNPs, "$inVCFfile") or die $!;
my $totalvcfSNPs = `grep -v "^#" $inVCFfile |wc -l`; chomp $totalvcfSNPs;
print $totalvcfSNPs,"\n";
foreach my $eachVCFline (<inputSNPs>){
	my @eachVCFtotal = split("\t", $eachVCFline);
	my $chrNo = $eachVCFtotal[0]; $chrNo =~ s/Chr//i;
	my $chrPos= $eachVCFtotal[1]; 
#	print $tailNo,"\n";
	my $dataSNPlist = `grep -m 1 "^$chrNo\t$chrPos" $dataSNPs `;chomp $dataSNPlist;
	if ($dataSNPlist){
#		print $tempLineNo, "\n";
		open(alleleInfo, "fileTemp");
		my @dataSNParray = <alleleInfo>;
		for(my $i = 2; $i < @allPopIDs;$i ++){
			if($dataSNParray[$i] == 1){
#				print $dataSNParray[$i], "\n";
				$ScoreList{$allPopIDs[$i]} ++;
#				print $ScoreList{$allPopIDs[$i]}, "\n"; 
			}
		}
	}
}

open (outScore, ">ScoreOut.txt");
for(my $i = 2; $i < @allPopIDs; $i ++){
	my $totaldataSNPs = `cut -f $i -d "," $dataSNPs|tail -n +2 |awk '{sum += \$0} END {print sum}'`;
	chomp $totaldataSNPs;
	my $probScore = 100*($ScoreList{$allPopIDs[$i]}*$ScoreList{$allPopIDs[$i]})/($totalvcfSNPs*$totaldataSNPs);
	print outScore $allPopIDs[$i],"\t", $probScore,"\n";
}

