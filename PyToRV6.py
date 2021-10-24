from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects.vectors import StrVector
import rpy2.rinterface as ri

def convert_file(*args, **kargs):
    #Import into python space
    SeqArray = importr("SeqArray")
    base = importr('base')

    #Open gds file as a gds object
    file = SeqArray.seqOpen(args[0])

    #Load sample names into memory
    sample = SeqArray.seqGetData(file, "sample.id")

    #moves the samples from the python space to the R space. This two step process is done
    #so we don't need to load the libraries into both spaces
    robjects.r('sampleID = {sample}'.format(sample = sample.r_repr()))

    #finds the total number of variants
    numOfVariants = SeqArray.seqSummary(file, varname = "variant.id")[0]

    #make the output function in R. This will take in a variant and parse the relevant 
    #information. It will output a list of strings where every string has the following
    #information: chrom pos ref alt hom/het
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
    #we now move the prior function from the R space to the python space
    output = robjects.globalenv['output']

    #im not sure how this works but this little ri.rternalize allows me to pass
    #this python function as an R function
    @ri.rternalize
    #Note: data is one variant node in the gds file. This will serve as the callback
    #function for seqApply 
    def printWithStops(data):
        result = []
        #calls the python wrapped, R function 'output'. Hence 'lines' is a nicely parsed
        #list of strings
        lines = output(data)
        #for reasons that I should figure out but have yet to do so, lines appears to be some
        #matrix data structure that I'm not familiar with in python and hence I do this weird
        #extra set of parsing (although it may not be required) to generate a new list of strings
        rows = len(lines) // 6
        for line in range(rows):
            result.append(" ".join([lines[line], lines[line+ rows], lines[line+ (2* rows)], lines[line+ (3* rows)], lines[line+ (4* rows)], lines[line+ (5* rows)]]))
        #converting result to a strvector is crucial because rpy2 will automatically convert a strvector
        #into a R compatiable data type (probably a vector of strings) but won't do this
        #for a general array. This is relevant because if seqApply (an R function) outputs an array of strings
        #(a non R type) then we get an error
        return StrVector(result)

    #edge case, if we have sufficiently small amount of variants just compute everything
    if numOfVariants <= 10:
        result = []
        SeqArray.seqApply(file, StrVector(["position", "position", "allele", "genotype"]), FUN= printWithStops, margin="by.variant", as_is="list")
    #else, we buffer and call seqApply multiple times
    else:
        prev = 1
        #this is a bit lazy of me but the R seq() function is similar to the python
        #range() but has some properties that are nice for this application and hence
        #I use it over range() 
        for i in  robjects.r.seq(10, numOfVariants, 10):
            #here we are telling the seqArray library to only consider the next 10 variants
            SeqArray.seqResetFilter(file, verbose = False)
            robjects.r('array = c({previous}:{next})'.format(previous = prev, next = i))
            SeqArray.seqSetFilter(file, variant_sel = robjects.globalenv['array'], verbose = False)

            #all the information for the next 10 variants have been extracted and parsed as an array
            #of vectors of strings (or maybe a vector of vectors of strings but the difference doesn't matter
            # now because vectors and arrays are both iterable)
            result = SeqArray.seqApply(file, StrVector(["position", "position", "allele", "genotype"]), FUN= printWithStops, margin="by.variant", as_is="list")
            prev = i
            #go through the variants and then return the lines one by one
            for variant in result:
                for line in variant:
                    yield line

for line in convert_file('OSSE\hapmap_r23a.gds'):
    input("...")
    print(line)

#TODO: check all my loops and ranges. I may (or may not) be repeatedly counting some variants 
#c(prev, prev + 10) => (next cycle) c(prev + 10, prev+20)
#since the first value is inclusive then I may be overcounting every tenth variant 

#TODO: currently I buffer every 10 variants. Make this variable

#TODO: In reference to the last todo, make the buffer size automatic based on total size of
# file and maybe also in refernece to pc specs

#TODO: on line 103 I call reset filter everytime, I think this is necessary but it may not be
# more testing/research is required

#TODO: on line 84 I parse my data again, this may not be needed