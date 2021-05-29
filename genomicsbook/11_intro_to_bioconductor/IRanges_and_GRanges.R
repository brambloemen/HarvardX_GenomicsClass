 
### IRanges and GRanges
library(IRanges)
library(GenomicRanges)


# IRanges -----------------------------------------------------------------


ir <- IRanges(5,10)
start(ir)
end(ir)
width(ir)

IRanges(start=c(3,5,17), end=c(10,8,20))

# shift range
ir
shift(ir, -2)

# other operations
narrow(ir, start=2) # start at 2nd number in range

narrow(ir, end=5) #end at 5th number in range

ir
flank(ir, width=3, start=TRUE, both=FALSE) #take a flank (outside the range)
flank(ir, width=-3, start=TRUE, both=FALSE) #take a flank (inside the range

flank(ir, width=3, start=FALSE, both=FALSE)
flank(ir, width=3, start=TRUE, both=TRUE) #both: take flank on both sides of start/end

ir
ir * 2
ir * -2
ir +2
ir -2


# set up a plotting window so we can look at range operations
plot(0,0,xlim=c(0,23),ylim=c(0,13),type="n",xlab="",ylab="",xaxt="n")
axis(1,0:15)
abline(v=0:14 + .5,col=rgb(0,0,0,.5))
# plot the original IRange
plotir <- function(ir,i) { arrows(start(ir)-.5,i,end(ir)+.5,i,code=3,angle=90,lwd=3) }
plotir(ir,1)
# draw a red shadow for the original IRange
polygon(c(start(ir)-.5,start(ir)-.5,end(ir)+.5,end(ir)+.5),c(-1,15,15,-1),col=rgb(1,0,0,.2),border=NA)
# draw the different ranges
plotir(shift(ir,-2), 2)
plotir(narrow(ir, start=2), 3)
plotir(narrow(ir, end=5), 4)
plotir(flank(ir, width=3, start=TRUE, both=FALSE), 5)
plotir(flank(ir, width=3, start=FALSE, both=FALSE), 6)
plotir(flank(ir, width=3, start=TRUE, both=TRUE), 7)
plotir(ir * 2, 8)
plotir(ir * -2, 9)
plotir(ir + 2, 10)
plotir(ir - 2, 11)
plotir(resize(ir, 1), 12)
text(rep(15,12), 1:12, c("ir","shift(ir,-2)","narrow(ir,start=2)",
                         "narrow(ir,end=5)",
                         "flank(ir, start=T, both=F)",
                         "flank(ir, start=F, both=F)",
                         "flank(ir, start=T, both=T)",
                         "ir * 2","ir * -2","ir + 2","ir - 2",
                         "resize(ir, 1)"), pos=4)


# inter range operations
(ir <- IRanges(start=c(3,5,17), end=c(10,8,20)))
range(ir) #total range, gaps included
reduce(ir) #remove gaps, overlaps
gaps(ir) #ranges of gaps
disjoin(ir) #make non-overlapping ranges, but including all start/endpoints


# GRanges -----------------------------------------------------------------

gr <- GRanges("chrZ", IRanges(start=c(5,10),end=c(35,45)),
             strand="+", seqlengths=c(chrZ=100L))
gr


genome(gr) <- "hg19"
gr

seqnames(gr)
seqlengths(gr)

shift(gr, 10)

shift(gr, 80) #gives warning: range is shifted beyond sequence boundaries
trim(shift(gr, 80)) #trim so that range stays within sequence

mcols(gr)
mcols(gr)$value <- c(-1,4)
gr
mcols(gr)$value <- NULL


gr2 <- GRanges("chrZ",IRanges(11:13,51:53))
grl <- GRangesList(gr, gr2)
grl
length(grl)
elementNROWS(grl)
width(grl)
sum(width(grl))

# we can add metadata: each row provides metadata for the corresponding GRrange object
mcols(grl)$value <- c(5,7)
mcols(grl)
grl

# overlaps
(gr1 <- GRanges("chrZ",IRanges(c(1,11,21,31,41),width=5),strand="*"))
(gr2 <- GRanges("chrZ",IRanges(c(19,33),c(38,35)),strand="*"))

fo <- findOverlaps(gr1,gr2) #strand specific, but there is an option to ignore strand
fo #interpretation: 
# 3th element of query (gr1, range 21-25) overlapped with first element of subject (gr2, range 19-38)
# 4th element of query with both the first and second element of subject, ...
queryHits(fo)
subjectHits(fo)

# similarly: %over% function
gr1 %over% gr2 #for each range in gr1, indicate whether an overlap has been found


# RLE: run length encoding
(r <- Rle(c(1,1,1,0,0,-2,-2,-2,rep(-1,20))))
str(r)
as.numeric(r)

# views object: window to look into rle object (e.g. large sequences)
(v <- Views(r, start=c(4,2), end=c(7,6)))



# Applications with genomic elements: strand-aware operations -------------

ir <- IRanges(c(3, 8, 14, 15, 19, 34, 40),
              width = c(12, 6, 6, 15, 6, 2, 7))
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...)
{
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}
plotGRanges = function (x, xlim = x, col = "black", sep = 0.5, xlimits = c(0, 
                                                                           60), ...) 
{
  main = deparse(substitute(x))
  ch = as.character(seqnames(x)[1])
  x = ranges(x)
  height <- 1
  if (is(xlim, "Ranges")) 
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim = xlimits, c(0, max(bins) * (height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom + height, 
       col = col, ...)
  title(main, xlab = ch)
  axis(1)
}
par(mfrow=c(4,1), mar=c(4,2,2,2))
plotRanges(ir, xlim=c(0,60))
plotRanges(reduce(ir), xlim=c(0,60))
plotRanges(disjoin(ir), xlim=c(0,60))
plotRanges(gaps(ir), xlim=c(0,60))

# extension to GRanges
gir = GRanges(seqnames="chr1", ir, strand=c(rep("+", 4), rep("-",3)))

par(mfrow=c(4,1), mar=c(4,2,2,2))
plotGRanges(gir, xlim=c(0,60))
plotGRanges(resize(gir,1), xlim=c(0,60),col="green")
plotGRanges(flank(gir,3), xlim=c(0,60), col="purple")
plotGRanges(flank(gir,2,start=FALSE), xlim=c(0,60), col="brown")

# application to visualize methylation array data
# obtain the data from a methylation array experiment
library(ArrayExpress)
if (!file.exists("E-MTAB-5797.sdrf.txt")) nano = getAE("E-MTAB-5797")
library(minfi)
pref = unique(substr(dir(patt="idat"),1,17)) # find the prefix strings: dir identifies all files in directory which match the patt="idat" regex
raw = read.metharray(pref)
glioMeth = preprocessQuantile(raw) # generate SummarizedExperiment

# visualize
MbyGene = function(mset, symbol="TP53", rad=5000) {
  # phase 1: annotated GRanges for the gene
  require(erma)
  require(Gviz)
  gmod = suppressMessages(genemodel(symbol))     # erma utility
  gseq = as.character(seqnames(gmod)[1])
  gmod$transcript = symbol
  # phase 2: filter down to the region of interest
  mlim = mset[which(seqnames(mset)==gseq),] # restrict to chromosome
  # now focus the methylation data to vicinity of gene
  d1 = subsetByOverlaps(GRanges(rowRanges(mlim),,, getM(mlim)), 
                        range(gmod)+rad)
  # phase 3: use Gviz
  plotTracks(list(DataTrack(d1), 
                  GeneRegionTrack(gmod, 
                                  transcriptAnnotation="transcript", name=gseq), 
                  GenomeAxisTrack(name=gseq, showTitle=TRUE)))
}

MbyGene(glioMeth, symbol="TERT")
