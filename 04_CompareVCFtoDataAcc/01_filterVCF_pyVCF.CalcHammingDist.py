#!/usr/bin/python
import math
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
import vcf
import re
from pygwas.core import genotype

def filterVCF(vcffile,qthres):
  vcf_f = vcf.Reader(open(vcffile,'r'))
  tarchr = numpy.empty(0)
  tarpos = numpy.empty(0)
  targt = numpy.empty(0)
  tarqual = numpy.empty(0)
  for eachPos in vcf_f:
#    if eachPos.samples[0]['GT'] != 'None':
#      if eachPos.QUAL > qthres:
    chr = int(re.sub('Chr','',eachPos.CHROM)) - 1
    tarchr = numpy.append(tarchr, chr)
    tarpos = numpy.append(tarpos, eachPos.POS)
    targt = numpy.append(targt, eachPos.samples[0]['GT'])
    tarqual = numpy.append(tarqual, eachPos.QUAL)
  return (tarchr, tarpos, targt, tarqual)

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-i", "--vcf_file", dest="inFile", help="Input VCF file", type="string")
inOptions.add_option("-t", "--qual_thres", dest="qualThres", help="Filter VCF file for quality threshold", type="float")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-n", "--file_num_snps", dest="file_num_snps", help="Output from the CalculateSNPseachAcc.py script", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")

(options, args) = inOptions.parse_args()


GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)

(TargetChr, TargetPos, TargetGT, TargetQUAL) = filterVCF(options.inFile, options.qualThres)

NumSNPs = len(TargetPos)
print "Total number of SNPs scanned", NumSNPs

MatchedSNPs = 0
ScoreList = numpy.zeros(len(GenotypeData.accessions))
# We know there is only one sample in the gVCF file
for chr in range(0,len(GenotypeData.chr_regions)):
  tarchrind = numpy.where(TargetChr == chr)[0]
  tarchrpos = TargetPos[tarchrind]
  start = GenotypeData.chr_regions[chr][0]
  end = GenotypeData.chr_regions[chr][1]
  positions = GenotypeData.positions[start:end]
  matsnpind = numpy.where(numpy.in1d(positions, tarchrpos))[0] + start
  MatchedSNPs = MatchedSNPs + len(matsnpind)
  tupMatSNPind = tuple(matsnpind)
  ScoreList = ScoreList + numpy.sum(GenotypeData.snps[tupMatSNPind,], axis = 0)

print "Number of matched SNPs", MatchedSNPs

# Calculate the number of SNPs in all the accessions
filetotal = open(options.file_num_snps, 'r')
numberSNPs = filetotal.read().split("\n")
filetotal.close()
# Calculate the final probability based on the score count and total count
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  numsnp = numberSNPs[i].split("\t")[1]
  outScore = (float(ScoreList[i])*float(ScoreList[i]))/(int(numsnp)*int(NumSNPs))
  outfile.write(GenotypeData.accessions[i])
  outfile.write("\t")
  outfile.write("%s" % int(ScoreList[i]))
  outfile.write("\t")
  outfile.write("%s" % int(numsnp))
  outfile.write("\t")
  outfile.write("%s\t" % outScore)
  outfile.write("%s\n" % NumSNPs)
outfile.close()
