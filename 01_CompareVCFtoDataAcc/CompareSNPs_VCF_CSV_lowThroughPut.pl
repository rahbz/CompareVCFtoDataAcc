#!/usr/bin/perl
use warnings;
use strict;

# Written by Rahul Bharadwaj
# Check the given SNPs in the CSV table ???


my $dataSNPs = $ARGV[0]; chomp $dataSNPs;
my $inVCFfile= $ARGV[1]; chomp $inVCFfile;

# Make sure the database SNP table is always sorted.....
#system("(head -n 1 $dataSNPs && tail -n +2 $dataSNPs |sort -k1,2n) > temp_dataSNPs.sorted");

my $headDataSNPs = `head -n 1 $dataSNPs`;chomp $headDataSNPs;
my @allPopIDs = split(',', $headDataSNPs);
# Create a bash generated search pattern file
# Create a hash for the all the populations in the DATA SNPs
my %ScoreList ;
for (my $i = 2; $i < @allPopIDs ; $i++){
	$ScoreList{$allPopIDs[$i]} = 0;
}

open (inputSNPs, "$inVCFfile") or die $!;
my $tempLineNo = 0;
system("cp $dataSNPs temp_dataSNPs1");
foreach my $eachVCFline (<inputSNPs>){
	my $lineNo = $tempLineNo + 1;
	my @eachVCFtotal = split("\t", $eachVCFline);
	my $chrNo = $eachVCFtotal[0]; $chrNo =~ s/Chr//i;
	my $chrPos= $eachVCFtotal[1]; 
	my $tailNo = "+"."$lineNo";
	system("tail -n $tailNo temp_dataSNPs1 > temp_dataSNPs2");
#	print $tailNo,"\n";
	my $dataSNPlist = `grep -n -m 1 "^$chrNo,$chrPos," temp_dataSNPs2`;chomp $dataSNPlist;
	if ($dataSNPlist){
		my @splitSNPlist = split(":", $dataSNPlist);
		$tempLineNo = $splitSNPlist[0];
#		print $tempLineNo, "\n";
		my @dataSNParray = split(",", $splitSNPlist[1]);
		for(my $i = 2; $i < @allPopIDs;$i ++){
			if($dataSNParray[$i] == 1){
#				print $dataSNParray[$i], "\n";
				$ScoreList{$allPopIDs[$i]} ++;
#				print $ScoreList{$allPopIDs[$i]}, "\n"; 
			}
		}
	}
	system("cp temp_dataSNPs2 temp_dataSNPs1");
}

open (outScore, ">ScoreOut.txt");
foreach my $k (keys %ScoreList){
	print outScore $k,"\t", $ScoreList{$k},"\n";
}
