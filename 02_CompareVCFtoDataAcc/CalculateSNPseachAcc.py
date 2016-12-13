#!/usr/bin/python
import sys
import math
from subprocess import call
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
from pygwas.core import genotype

GenotypeData = genotype.load_hdf5_genotype_data('/lustre/scratch/users/rahul.pisupati/all_chromosomes_binary.hdf5')

# Calculate the number of SNPs in all the accessions
# Takes a really long time

# Calculate the final probability based on the score count and total count
#outfile = open("totalSNPsNum_1001genomes.txt", 'w')
for i in range(0, len(GenotypeData.accessions)):
  outScore = sum(GenotypeData.snps[:, i])
  print GenotypeData.accessions[i], "\t", outScore
#  print "Written count for", i+1, "accessions", "Accession:", GenotypeData.accessions[i], "Count:", outScore

#outfile.close()
