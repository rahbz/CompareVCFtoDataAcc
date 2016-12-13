#!/usr/bin/python
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas
# module load h5py

import numpy
import vcfnp
import logging
from pygwas.core import genotype





#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-i", "--vcf_file", dest="vcfFile", help="VCF file for the sample", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-a", "--p1", dest="p1", help="Expected parent 1", type="string")
inOptions.add_option("-b", "--p2", dest="p2", help="Expected parent 2", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")

(options, args) = inOptions.parse_args()

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
num_lines = len(GenotypeData.accessions)

logging.info("Reading the VCF file")
vcf = vcfnp.variants(options.vcfFile, cache=True).view(numpy.recarray)
vcfD = vcfnp.calldata_2d(options.vcfFile, cache=True).view(numpy.recarray)


p1i = numpy.where(GenotypeData.accessions == options.p1)[0]
p2i = numpy.where(GenotypeData.accessions == options.p2)[0]

p1 = GenotypeData_acc.snps[:,p1i]
p2 = GenotypeData_acc.snps[:,p2i]

DPthres = numpy.mean(vcf.DP[numpy.where(vcf.DP > 0)[0]]) * 2
print "Threshold for depth is set at: ", DPthres

snpsREQ = numpy.where((vcfD.is_called[:,0]) & (vcf.QUAL > 30) & (vcf.DP > 0) & (vcf.DP < DPthres))[0]
snpCHR = numpy.array(numpy.chararray.replace(vcf.CHROM[snpsREQ], "Chr", "")).astype("int8")
snpPOS = numpy.array(vcf.POS[snpsREQ])
snpGT = vcfD.GT[snpsREQ, 0]   ## since one sample 
snpPL = vcfD.PL[snpsREQ, 0]
snpGQ = vcfD.GQ[snpsREQ, 0]
snpDP = vcf.DP[snpsREQ]

NumMatSNPs = 0
chunk_size = 1000
Likelihood = 0

TotalmatTarSNPind = numpy.zeros(0, dtype="int")
TotalmatAccInd = numpy.zeros(0, dtype="int")

for i in range(1,6):
  perchrTarPos = numpy.where(snpCHR == i)[0]
  perchrtarSNPpos = snpPOS[perchrTarPos]
  logging.info("Loaded %s chromosome positions from the position file", i)
  start = GenotypeData.chr_regions[i-1][0]
  end = GenotypeData.chr_regions[i-1][1]
  chrpositions = GenotypeData.positions[start:end]
  matchedAccInd = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
  matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, chrpositions))[0]
  matchedTarPL = snpPL[perchrTarPos[matchedTarInd],]
  matchedTarDP = snpDP[perchrTarPos[matchedTarInd]]
 
  TotalmatTarSNPind = numpy.append(TotalmatTarSNPind, perchrTarPos[matchedTarInd]) 
  TotalmatAccInd = numpy.append(TotalmatAccInd, matchedAccInd)

  gtp1 = p1[matchedAccInd]
  gtp2 = p2[matchedAccInd]

  homalt = numpy.where((gtp1 == 1) & (gtp2 == 1))[0]
  homref = numpy.where((gtp1 == 0) & (gtp2 == 0))[0]
  het = numpy.where((gtp1 != -1) & (gtp2 != -1) & (gtp1 != gtp2))[0]
  
  Likelihood = Likelihood + numpy.sum(matchedTarPL[homalt, 2]) + numpy.sum(matchedTarPL[homref, 0]) + numpy.sum(matchedTarPL[het, 1])
  NumMatSNPs = NumMatSNPs + len(matchedAccInd)


print float(Likelihood)/NumMatSNPs, "\t", NumMatSNPs
  
