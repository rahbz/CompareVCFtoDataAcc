#!/bin/python


import numpy
import vcf
from pygwas.core import genotype
from optparse import OptionParser

#Takes in the matched SNPs file and the VCF file
#  and calculate the allele frequency at each position
# may be use vcf tools 

#module load PyVCF

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-m", "--matched_SNPs_list", dest="matchedSNPs", help="Input the matched SNPs *.matchedSNPs", type="string")
inOptions.add_option("-p", "--pos_file", dest="posFile", help="Input position file from the VCF file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")

(options, args) = inOptions.parse_args()

GenotypeData = genotype.load_hdf5_genotype_data('/lustre/scratch/users/rahul.pisupati/all_chromosomes_binary.hdf5')

matchedSNPs = open(options.matchedSNPs, 'r')
posfile = open(options.posFile, 'r')
TargetSNPs = tuple(posfile.read().rstrip().split("\n"))

for matsnp in matchedSNP.readlines()[0:]:
  modLine = tuple(map(float, matsnp.rstrip().split(",")))
  for eachsnp in modLine:
    for i in range(0,5):
      if (eachsnp > GenotypeData.chr_regions[i][0] and eachsnp < GenotypeData.chr_regions[i][1]):
        chrNo = i
        break
    chrPos = GenotypeData.positions[eachsnp]
    

