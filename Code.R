#https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/methylKit/inst/doc/methylKit.html
#install "methylKit" and "genomation"



if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("methylKit")



if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("genomation")


library(methylKit)
library(genomation) 
library(GenomicRanges)

#Reading the methylation call files

#This code uses defult data you can copy and paste your own data in C:\Users\00090473\AppData\Local\R\win-library\4.2\methylKit\extdata 

library(methylKit)
file.list=list( system.file("extdata", 
                            "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata",
                            "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", 
                            "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", 
                            "control2.myCpG.txt", package = "methylKit") )


# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=list("test1","test2","ctrl1","ctrl2"),
           assembly="hg18",
           treatment=c(1,1,0,0),
           context="CpG"
           )



#Reading the methylation call files and store them as flat file database

library(methylKit)
file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )


# read the files to a methylRawListDB object: myobjDB 
# and save in databases in folder methylDB
myobjDB=methRead(file.list,
           sample.id=list("test1","test2","ctrl1","ctrl2"),
           assembly="hg18",
           treatment=c(1,1,0,0),
           context="CpG",
           dbtype = "tabix",
           dbdir = "methylDB"
           )

print(myobjDB[[1]]@dbpath)


#Reading the methylation calls from sorted Bismark alignments

my.methRaw=processBismarkAln( location = 
                                system.file("extdata",
                                                "test.fastq_bismark.sorted.min.sam", 
                                                  package = "methylKit"),
                         sample.id="test1", assembly="hg18", 
                         read.context="CpG", save.folder=getwd())

#Descriptive statistics on samples

getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)


getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)


getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)


#Filtering samples based on read coverage

filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)


#Comparative analysis
#Merging samples


meth=unite(myobj, destrand=FALSE)

head(meth)

meth.min=unite(myobj,min.per.group=1L)

getCorrelation(meth,plot=TRUE)


#Clustering samples

clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)


PCASamples(meth, screeplot=TRUE)


PCASamples(meth)


#Batch effects

# make some batch data frame
# this is a bogus data frame
# we don't have batch information
# for the example data
sampleAnnotation=data.frame(batch_id=c("a","a","b","b"),
                            age=c(19,34,23,40))

as=assocComp(mBase=meth,sampleAnnotation)
as

# construct a new object by removing the first pricipal component
# from percent methylation value matrix
newObj=removeComp(meth,comp=1)

mat=percMethylation(meth)

# do some changes in the matrix
# this is just a toy example
# ideally you want to correct the matrix
# for batch effects
mat[mat==100]=80
 
# reconstruct the methylBase from the corrected matrix
newobj=reconstruct(mat,meth)


#Tiling windows analysis

tiles=tileMethylCounts(myobj,win.size=1000,step.size=1000)
head(tiles[[1]],3)

#Finding differentially methylated bases or regions

myDiff=calculateDiffMeth(meth)


# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
#
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
#
#
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)


#Correcting for overdispersion

sim.methylBase1<-dataSim(replicates=6,sites=1000,
                         treatment=c(rep(1,3),rep(0,3)),
                        sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
                        )

my.diffMeth<-calculateDiffMeth(sim.methylBase1[1:,],
                                overdispersion="MN",test="Chisq",mc.cores=1)


#Accounting for covariates

covariates=data.frame(age=c(30,80,34,30,80,40))
sim.methylBase<-dataSim(replicates=6,sites=1000,
                        treatment=c(rep(1,3),rep(0,3)),
                        covariates=covariates,
                        sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
                        )
my.diffMeth3<-calculateDiffMeth(sim.methylBase,
                                covariates=covariates,
                                overdispersion="MN",test="Chisq",mc.cores=1)


#Finding differentially methylated bases using multiple-cores

myDiff=calculateDiffMeth(meth,mc.cores=2)


#Annotating differentially methylated bases or regions

library(genomation)

# read the gene BED file
gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                           package = "methylKit"))


#
# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data
#
annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                        package = "methylKit"),
                           feature.flank.name=c("CpGi","shores"))
#
# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")


#Regional analysis

promoters=regionCounts(myobj,gene.obj$promoters)

head(promoters[[1]])

#Convenience functions for annotation objects
diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

# target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn))

getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)


plotTargetAnnotation(diffAnn,precedence=TRUE,
    main="differential methylation annotation")


plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
       main="differential methylation annotation")


getFeatsWithTargetsStats(diffAnn,percentage=TRUE)


#methylKit convenience functions
#Coercing methylKit objects to GRanges

class(meth)

as(meth,"GRanges")

class(myDiff)

as(myDiff,"GRanges")


#Converting methylKit objects to methylDB objects and vice versa

class(myobjDB[[1]])


as(myobjDB[[1]],"methylRaw")


data(methylKit)
 
objDB=makeMethylDB(methylBase.obj,"exMethylDB")


#Selection: subsetting methylKit Objects

select(meth,1:5) # get first 10 rows of a methylBase object

myDiff[21:25,] # get 5 rows of a methylDiff object

#selectByOverlap

library(GenomicRanges)
my.win=GRanges(seqnames="chr21",
               ranges=IRanges(start=seq(from=9764513,by=10000,length.out=20),width=5000) )
 
# selects the records that lie inside the regions
selectByOverlap(myobj[[1]],my.win)


#reorganize(): reorganizing samples and treatment vector within methylKit objects

# creates a new methylRawList object
myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
# creates a new methylBase object
meth2 =reorganize(meth,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )


#percMethylation(): Getting percent methylation matrix from methylBase objects

# creates a matrix containing percent methylation values
perc.meth=percMethylation(meth)


#methSeg(): segmentation of methylation or differential methylation profiles

mbw=readRDS("H1.chr21.chr22.rds")  # put this file in your WD

 # it finds the optimal number of componets as 6
 res=methSeg(mbw,diagnostic.plot=TRUE,maxInt=100,minSeg=10)

 # however the BIC stabilizes after 4, we can also try 4 componets
 res=methSeg(mbw,diagnostic.plot=TRUE,maxInt=100,minSeg=10,G=1:4)

 # get segments to BED file
 methSeg2bed(res,filename="H1.chr21.chr22.trial.seg.bed")

