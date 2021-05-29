.libPaths('C:/Users/BrBl1834/R/win-library')


# General considerations --------------------------------------------------

library(limma)
library(Biobase)
library(data.table)
library(GEOquery)
library(erma)
library(RNAseqData.HNRNPC.bam.chr14)
library(airway)
library(annotate)
library(minfi)
library(locfit)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(BiocStyle)
library(GenomicFiles)
library(GenomicAlignments)
library(MultiAssayExperiment)
library(RaggedExperiment)
library(VariantAnnotation)
library(VariantTools)
library(bigrquery)
library(dplyr)
library(magrittr)
library(curatedTCGAData)
library(ArrayExpress)

methods(class="SummarizedExperiment")
browseVignettes("SummarizedExperiment")

library(airway)
data(airway)
airway

metadata(airway)

colData(airway)
rowRanges(airway)

rowRanges(airway)$ENSG00000172057

# see to how many unique exons/ranges the diff rows can be reduced
reduce(rowRanges(airway)$ENSG00000172057)


names(colData(airway))

table(airway$dex) # main treatment factor



# Handling the ArrayExpress deposit of Illumina 450k Methylation arrays --------


# run previous script from course to load the ArrayExpress for next part of exercise
# source("./management_genome_scale_exp.R")
# library(data.table)
# sd5797 = fread("E-MTAB-5797.sdrf.txt")
# head(sd5797[,c(3,16,18)])



# External storage of large assay data - HDF5Array, saveHDF5Summar --------
## Measuring memory consumption
gc()

## Demonstrating HDF5 for external storage
library(HDF5Array)

data(airway)
airass = assay(airway)  # obtain numerical data, then save as HDF5
href = writeHDF5Array(airass, "airass.h5", "airway")

## HDF5-backed SummarizedExperiment
saveHDF5SummarizedExperiment(airway, "externalAirway", replace=TRUE)
newse = loadHDF5SummarizedExperiment("externalAirway")
newse

assay(newse[c("ENSG00000000005", "LRG_99"), 
            which(newse$dex == "trt")]) # use familiar subsetting

gc()



# GenomicFiles: families of files of a given type -------------------------

library(RNAseqData.HNRNPC.bam.chr14)
library(GenomicFiles)
gf = GenomicFiles(files=RNAseqData.HNRNPC.bam.chr14_BAMFILES)
gf

# define GRanges object that contains the HNRNPC region of interest
hn = GRanges("chr14", IRanges(21677296, 21737638), strand="-")
rowRanges(gf) = hn

# 
library(GenomicAlignments)
# define MAP function: extract alignments that overlaps the region of interest
MAP = function(r, f) 
  # read genomic alignments from BAM file
  readGAlignmentPairs(f, param=ScanBamParam(which=r))
ali = reduceByRange(gf, MAP=MAP)

# length=number of reads aligned to the specified ROI (defined by hn or rowranges(gf)=GRanges object)
sapply(ali[[1]], length)



# BED collections ---------------------------------------------------------


