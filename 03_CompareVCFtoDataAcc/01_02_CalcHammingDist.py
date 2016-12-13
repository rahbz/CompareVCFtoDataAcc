#!/usr/bin/python
import math
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
from pygwas.core import genotype


#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-m", "--matched_pos", dest="inFile", help="Input Matched Pos file from Rscript", type="string")
inOptions.add_option("-n", "--num_snps_vcf", dest="numsnps", help="Number of SNPs taken for analysis from VCF file", type="float")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-t", "--file_num_snps", dest="file_num_snps", help="Output from the CalculateSNPseachAcc.py script", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")

(options, args) = inOptions.parse_args()


GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)

# Create a numpy array containing all the positions

NumSNP = int(options.numsnps)
ScoreList = numpy.zeros(len(GenotypeData.accessions))

matchedSNP = open(options.inFile, 'r')
for matsnp in matchedSNP.readlines()[0:]:
  modLine = tuple(map(float, matsnp.rstrip().split(",")))
  ScoreList = ScoreList + numpy.sum(GenotypeData.snps[modLine,], axis = 0)
    
print "Total number of SNPs scanned", NumSNP
#print ScoreList
# Calculate the number of SNPs in all the accessions
filetotal = open(options.file_num_snps, 'r')
numberSNPs = filetotal.read().split("\n")
filetotal.close()
# Calculate the final probability based on the score count and total count
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
#  print getsnp
  numsnp = numberSNPs[i].split("\t")[1]
#  outScore = float(ScoreList[i])/int(numsnp)
  outScore = (float(ScoreList[i])*float(ScoreList[i]))/(int(numsnp)*int(NumSNP))
  outfile.write(GenotypeData.accessions[i])
  outfile.write("\t")
  outfile.write("%s" % int(ScoreList[i]))
  outfile.write("\t")
  outfile.write("%s" % int(numsnp))
  outfile.write("\t")
  outfile.write("%s\n" % outScore)
outfile.close()
