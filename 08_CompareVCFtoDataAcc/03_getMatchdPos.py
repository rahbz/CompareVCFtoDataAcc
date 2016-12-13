#!/usr/bin/python
import math
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas
import logging
import numpy
import numpy.ma
import pandas
from pygwas.core import genotype
import scipy

inOptions = OptionParser()
inOptions.add_option("-p", "--pos_file", dest="posFile", help="Position file removing the header from VCF", type="string")
inOptions.add_option("-a", "--acc_to_check", dest="accToCheck", help = "Top hit accession from the ScoreAcc file", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-m", "--outMatFile", dest="outNonMat", help="Output positions for non-matched positions")
inOptions.add_option("-n", "--outNonMatFile", dest="outMat", help="Output positions for matched positions")
(options, args) = inOptions.parse_args()

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
num_lines = len(GenotypeData.accessions)
Acc = numpy.where(GenotypeData.accessions == options.accToCheck)[0]
AccSNPs = GenotypeData_acc.snps[:, Acc]


logging.debug("Reading the position file")
targetSNPs = pandas.read_table(options.posFile, header=None, usecols=[0,1,2])

TotNonMatPos = numpy.zeros(0, dtype="uint32")
TotNonMatChr = numpy.zeros(0, dtype="int8")
TotMatPos = numpy.zeros(0, dtype="uint32")
TotMatChr = numpy.zeros(0, dtype="uint32")
MatHetCount = 0
NmatHetCount = 0

for i in range(1,6):
  perchrTarInd = numpy.where(targetSNPs[0] == i)[0]
  perchrtarSNPpos = numpy.array(targetSNPs[1][perchrTarInd])
  logging.debug("Loaded %s chromosome positions from the position file", i)
  start = GenotypeData.chr_regions[i-1][0]
  end = GenotypeData.chr_regions[i-1][1]
  chrpositions = GenotypeData.positions[start:end]

  matchedAccInd = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
  matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, chrpositions))[0]
  matchedTarGTs = targetSNPs[2][perchrTarInd[matchedTarInd]].values
  TarGTs = numpy.zeros(len(matchedTarGTs), dtype="int8")
  TarGTs[numpy.where(matchedTarGTs != "0/0")[0]] = 1
  AccTarSNPs = AccSNPs[matchedAccInd]
  
  ImpPos = numpy.where(AccTarSNPs < 0)[0]
  nmat_all = numpy.where(AccTarSNPs != TarGTs)[0]
  matInd = numpy.where(AccTarSNPs == TarGTs)[0]
  nmatInd = numpy.setdiff1d(nmat_all, ImpPos)
  
  MatHetCount = MatHetCount + numpy.unique(numpy.array(matchedTarGTs[matInd]), return_counts=True)[1][1]
  NmatHetCount = NmatHetCount + numpy.unique(numpy.array(matchedTarGTs[nmatInd]), return_counts=True)[1][1]
  
  nmat = GenotypeData.positions[matchedAccInd[nmatInd]]
  mat = GenotypeData.positions[matchedAccInd[matInd]]

  logging.debug("%d NonMat, %d Mat, %s TotInfo", len(nmat), len(matInd), len(AccTarSNPs)-len(ImpPos))
  TotNonMatPos = numpy.append(TotNonMatPos, nmat)
  TotNonMatChr = numpy.append(TotNonMatChr, numpy.repeat(i, len(nmat)))
  TotMatPos = numpy.append(TotMatPos, mat)
  TotMatChr = numpy.append(TotMatChr, numpy.repeat(i, len(mat)))

print MatHetCount, "\t", NmatHetCount, "\t", len(TotMatPos), "\t", len(TotNonMatPos)
numpy.savetxt(options.outNonMat, numpy.column_stack((TotNonMatChr, TotNonMatPos)), fmt=['%d','%d'] )
numpy.savetxt(options.outMat, numpy.column_stack((TotMatChr, TotMatPos)), fmt=['%d','%d'] )


