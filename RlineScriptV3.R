library(gdsfmt)
library(SeqArray)

file.create("sample.txt")
output <- file("sample.txt")
open(output, "w")

f <- seqOpen("hapmap_r23a.gds")

sampleID = seqGetData(f, "sample.id")

returnData <- function(data) {
  chrom = data[[1]]
  pos = data[[2]]
  alleles <- c(strsplit(data[[3]][1], ",")[[1]], ".")
  genotype = data[[4]]
  
  internalWrite <- function(alt, hom_het, id) {
    writeLines(c(chrom, pos, alleles[1], alt, hom_het, id), output, sep = " ")
    #formatting but may slow down process
    writeLines("\n", output)
  }
  
  genotype[is.na(genotype)] = length(alleles) - 1
  
  for (col in 1:ncol(genotype)) {
    num1 = genotype[1, col] + 1
    num2 = genotype[2, col] + 1
    id = sampleID[col]
    
    if (num1 == 0 & num2 == 0) {
      #nothing to print
    }
    else if (num1 != 0 & num2 == 0) {
      internalWrite(alleles[num1], "het", id)
    }
    else if (num1 == 0 & num2 != 0) {
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
}

seqApply(f, c("chromosome", "position", "allele", "genotype"), FUN= returnData, margin="by.variant")

seqClose(f)
close(output)