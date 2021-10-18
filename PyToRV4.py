from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects.vectors import StrVector
import rpy2.rinterface as ri

class CravatConverter:
    def __init__(self):
        self.format_name = 'gds'

    def check_format(self, f):
        return f[len(f) - 4:] == '.gds'

    def setup(self, f):
        #R imports
        SeqArray = importr("SeqArray")
        base = importr('base')

        #Open file f as a seq file
        file = SeqArray.seqOpen(f)

        #Pull all the sample IDs into memory
        sampleID = SeqArray.seqGetData(file, "sample.id")

        #Throw all the sample IDs into the R space
        robjects.r('sampleID = {sample}'.format(sample = sampleID.r_repr()))

        #create an output function in the R space (this will be responsible for returning all the information of a given variant as a matrix)
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

        #get the output function from the R space and put it in the python space
        output = robjects.globalenv['output']

        #This is a function who's purpose is to wrap the output function. Using this we can output variants line by line 
        @ri.rternalize
        def printWithStops(data):
            lines = output(data)
            rows = len(lines) // 6
            for line in range(rows):
                input("...")
                print(lines[line] + "     " + lines[line + rows] + "     " + lines[line + (2* rows)] + "     " + lines[line + (3* rows)] + "     " + lines[line + (4* rows)] + "     " + lines[line + (5* rows)])


        #Apply function by variant 
        SeqArray.seqApply(file, StrVector(["position", "position", "allele", "genotype"]), FUN= printWithStops, margin="by.variant")

        #Close any opened files
        SeqArray.seqClose(file)

    def convert_line(self, l):
        pass
