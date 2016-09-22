source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("piano", dependencies=TRUE)
source("http://bioconductor.org/biocLite.R")
biocLite("GO.db")
source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")
 
library("piano")
library("Biobase")
 
 
## Load gene set file (source: http://www.broadinstitute.org/gsea/msigdb/collections.jsp)
## Follow the link above, look up: "BP: GO biological process" and thexn
## "Download GMT Files": entrez genes ids
myGsc <- loadGSC("C:/Users/adilm/Documents/R_documents/MyAnalysis/c5.bp.v4.0.entrez.gmt") #Your own path
head(myGsc)
 
## Load P-values and fold changes
DESeq_HCCvsMTCH = as.matrix(read.table("E:/CollaborativeProjects/Elias/HCC_Data_July2015/HCC-DE.txt")) # The DE file for HCC
DESeq_HCCvsMTCH = as.matrix(read.table("E:/CollaborativeProjects/Elias/HCC_Data_July2015/HCC-DE_ACSS1.txt")) # The DE file for ACSS1
DESeq_HCCvsMTCH = as.matrix(read.table("E:/CollaborativeProjects/Elias/HCC_Data_July2015/HCC-DE_ACSS2.txt")) # The DE file for ACSS1
 
DESeq_HCCvsMTCH = as.matrix(read.table("E:/CollaborativeProjects/Stefano/HCC-TM6SF2.txt")) # The DE file for TM6SF2
DESeq_HCCvsMTCH = as.matrix(read.table("E:/CollaborativeProjects/Stefano/HCC-PNPLA3.txt")) # The DE file for TM6SF2
 
head(DESeq_HCCvsMTCH)
pval= DESeq_HCCvsMTCH[ ,3] #extract pvalues
fc= DESeq_HCCvsMTCH[ ,2]  #extract fold changes
pval = as.matrix(pval)
fc = as.matrix(fc)
#put entrez ids as rownames
rownamn = DESeq_HCCvsMTCH[ ,1]
rownamn = as.matrix(rownamn)
rownames(pval) = paste(rownamn)
rownames(fc) = paste(rownamn)
#check everything OK:
head(pval)
head(fc)
 
## Replace NA with 0, complains otherwise
pval[is.na(pval)] <- 1
fc[is.na(fc)] <- 0
 
## Run GSEA
gsaRes <- runGSA(pval, fc, gsc=myGsc, geneSetStat="reporter", signifMethod="nullDist", gsSizeLim=c(10,10000), ncpus=4)
 
 
setwd('E:/CollaborativeProjects/Elias/HCC_Data_July2015')
## Write output to excel file
GSAsummaryTable(gsaRes, save=TRUE, file="E:/CollaborativeProjects/Elias/HCC_Data_July2015/PIANO_output_HCC-DE.xls")
GSAsummaryTable(gsaRes, save=TRUE, file="E:/CollaborativeProjects/Elias/HCC_Data_July2015/PIANO_output_ACSS1.xls")
GSAsummaryTable(gsaRes, save=TRUE, file="E:/CollaborativeProjects/Elias/HCC_Data_July2015/PIANO_output_ACSS2.xls")
 
setwd('E:/CollaborativeProjects/Stefano')
GSAsummaryTable(gsaRes, save=TRUE, file="E:/CollaborativeProjects/Stefano/PIANO_output_TM6SF2.xls")
GSAsummaryTable(gsaRes, save=TRUE, file="E:/CollaborativeProjects/Stefano/PIANO_output_PNPLA3.xls")
 
## network plot
par(mar = rep(2, 4))
networkPlot(gsaRes, class='non', significance=0.0000005,lay=5, ncharLabel=75)
 
## heatmap
dev.new(width=20,height=20)
GSAheatmap(gsaRes, cutoff=25, adjusted=FALSE, ncharLabel=75, cellnote="none", columnnames="full", colorkey=TRUE, colorgrad=NULL, cex=0.5)
 
