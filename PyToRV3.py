from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
import rpy2
import rpy2.rinterface as ri
import numpy as np
import time
import threading
from rpy2.robjects import numpy2ri

#Imports
SeqArray = importr("SeqArray")
base = importr('base')

#Open gds file as a gds object
file = SeqArray.seqOpen('OSSE\hapmap_r23a.gds')

#Load samples into memory
sam = SeqArray.seqGetData(file, "sample.id")
robjects.r('sampleID = {sample}'.format(sample = sam.r_repr()))

#Insert function here
robjects.r('''
output <- function(data) {
  lines <- matrix(, nrow = 0, ncol = 6)
  chrom = data[[1]]
  pos = data[[2]]
  alleles <- c(strsplit(data[[3]][1], ",")[[1]], ".")
  genotype = data[[4]]
  
  internalWrite <- function(alt, hom_het, id) {
    lines <<- rbind(lines, c(chrom, pos, alleles[1], alt, hom_het, id))
  }
  
  genotype[is.na(genotype)] = length(alleles) - 1
  
  for (col in 1:ncol(genotype)) {
    num1 = genotype[1, col] + 1
    num2 = genotype[2, col] + 1
    id = sampleID[col]
    
    if (num1 == 1 & num2 == 1) {
      #nothing to print
    }
    else if (num1 != 1 & num2 == 1) {
      internalWrite(alleles[num1], "het", id)
    }
    else if (num1 == 1 & num2 != 1) {
      internalWrite(alleles[num2], "het", id)
    }
    else if (num1 == num2) {
      internalWrite(alleles[num1], "hom", id)
    }
    else {
      internalWrite(alleles[num1], "het", id)
      internalWrite(alleles[num2], "het", id)
    }
  }
  return(lines)
}
''')
output = robjects.globalenv['output']

@ri.rternalize
def printWithStops(data):
  lines = output(data)
  rows = len(lines) // 6
  for line in range(rows):
      input("...")
      print(lines[line] + "     " + lines[line + rows] + "     " + lines[line + (2* rows)] + "     " + lines[line + (3* rows)] + "     " + lines[line + (4* rows)] + "     " + lines[line + (5* rows)])
  return 


#Apply function by variant 
SeqArray.seqApply(file, StrVector(["chromosome", "position", "allele", "genotype"]), FUN= printWithStops, margin="by.variant")

#Close any opened files
SeqArray.seqClose(file)