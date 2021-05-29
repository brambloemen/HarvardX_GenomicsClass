

# coordinating information: GSE5859Subset example -------------------------


load("../10_batch_effects/GSE5859Subset.rda")

dim(geneExpression)

dim(geneAnnotation)

dim(sampleInfo)

all.equal(sampleInfo$filename, colnames(geneExpression))

# geneExpression rows correspond to annotated genes in geneAnnotation
all.equal(rownames(geneExpression), geneAnnotation$PROBEID)

# sample info rows correspond to geneExpression columns
options(digits=2)
cbind(sampleInfo[1:3,], colnames(geneExpression)[1:3], 
      t(geneExpression)[1:3,1:4])


rownames(sampleInfo) = sampleInfo$filename
rownames(geneAnnotation) = geneAnnotation$PROBEID

# construct the ExpressionSet
library(Biobase)
es5859 = ExpressionSet(assayData=geneExpression)
pData(es5859) = sampleInfo
fData(es5859) = geneAnnotation
es5859

# example: use the fData (=feature data) to find 
es5859[which(fData(es5859)$CHR=="chrY"),]

# all methods that can act on expressionSet data
methods(class="ExpressionSet")


# endomorphism concept of expressionSet data: if G= vector of features, S=vector of samples, then we can do:
G <- which(fData(es5859)$CHR=="chrY")
S <- which(pData(es5859)$ethnicity=="ASN")
es5859[G,S] #all chrY genes of ASN ethnicity samples
# similar to subsetting matrices


# add annotation
annotation(es5859) = "hgfocus.db"

library(annotate)
mi = pmid2MIAME("17206142") #pubmed id
mi
abstract(mi)
experimentData(es5859) = mi
es5859

nchar(abstract(es5859))
substr(abstract(es5859),1,50)


# GEO (NIH), GEOquery  ---------------

library(GEOquery)
glioMA = getGEO("GSE78703")[[1]]

glioMA


# arrayexpress (EMBL) ------------------------------------------------------------
# note: arrayexpress is now moved into the Biostudies database

# query to ask for all glioblastoma 
library(ArrayExpress)
sets = queryAE(keywords = "glioblastoma", species = "homo+sapiens")
dim(sets)

sets[5:7,-c(7,8)]

initdir = dir()
if (!file.exists("E-MTAB-5797.sdrf.txt")) nano = getAE("E-MTAB-5797")

afterget = dir()
setdiff(afterget, initdir)






