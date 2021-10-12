from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
import rpy2
import numpy as np

#utils = rpackages.importr('utils')
#utils.chooseCRANmirror(ind=1)
#utils.install_packages('BiocManager')
#install = BiocManager.install
#install('SeqArray')
#BiocManager = importr('BiocManager')

#Imports
SeqArray = importr("SeqArray")
base = importr('base')

#Create and open output text file
base.file_create("OSSE\sample.txt")
output = base.file("OSSE\sample.txt")
base.open(output, "w")

#Open gds file as a gds object
file = SeqArray.seqOpen('OSSE\hapmap_r23a.gds')

#Load samples into memory
sampleID = SeqArray.seqGetData(file, "sample.id")

#Insert function here
def pyReturnData (data):
    pass
returnData = robjects.r('function(data) {\n'\
  'chrom = data[[1]]\n'\
  'pos = data[[2]]\n'\
  'alleles <- c(strsplit(data[[3]][1], ",")[[1]], ".")\n'\
  'genotype = data[[4]]\n'\
  '\n'\
  'internalWrite <- function(alt, hom_het, id) {\n'\
    'writeLines(c(chrom, pos, alleles[1], alt, hom_het, id), output, sep = " ")\n'\
    'writeLines("\n", output)\n'\
  '}\n'\
  '\n'\
  'genotype[is.na(genotype)] = length(alleles) - 1\n'\
  '\n'\
  'for (col in 1:ncol(genotype)) {\n'\
    'num1 = genotype[1, col] + 1\n'\
    'num2 = genotype[2, col] + 1\n'\
    'id = 0\n'\
    '\n'\
    'if (num1 == 0 & num2 == 0) {\n'\
    '}\n'\
    'else if (num1 != 0 & num2 == 0) {\n'\
      'internalWrite(alleles[num1], "het", id)\n'\
    '}\n'\
    'else if (num1 == 0 & num2 != 0) {\n'\
      'internalWrite(alleles[num2], "het", id)\n'\
    '}\n'\
    'else if (num1 == num2) {\n'\
      'internalWrite(alleles[num1], "hom", id)\n'\
    '}\n'\
    'else {\n'\
      'internalWrite(alleles[num1], "het", id)\n'\
      'internalWrite(alleles[num2], "het", id)\n'\
    '}\n'\
  '}\n'\
'}')

#Apply function variant by variant
SeqArray.seqApply(file, StrVector(["position", "position", "allele", "genotype"]), FUN= returnData, margin="by.variant")
#This breaks
SeqArray.seqApply(file, StrVector(["chromosome", "position", "allele", "genotype"]), FUN= base.print, margin="by.variant")

SeqArray.seqApply(file, "genotype", FUN= base.print, margin="by.variant")

#Close any opened files
SeqArray.seqClose(file)
base.close(output)


