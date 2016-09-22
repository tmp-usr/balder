#######################################
# DIFFERENTIAL EXPRESSION with DESeq
#######################################
 
#Install and load
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
 
#import data, every other row from the txt file containing un-normalized and normalized values.
 
#First import whole data matrix
#CountTable_O_4_4 = as.matrix(read.table("E:/CollaborativeProjects/BertNederland/TissueRNAseq/O_5-7-Human_HPA.txt", sep="\t", header=TRUE, row.names=1))
#CountTable_O_4_4 = as.matrix(read.table("E:/CollaborativeProjects/BertNederland/TissueRNAseq/L_5-7-Human.txt", sep="\t", header=TRUE, row.names=1))
#CountTable_O_4_4 = as.matrix(read.table("E:/CollaborativeProjects/BertNederland/TissueRNAseq/M_5-7-Human.txt", sep="\t", header=TRUE, row.names=1))
#CountTable_O_4_4 = as.matrix(read.table("E:/CollaborativeProjects/BertNederland/TissueRNAseq/O_5-7-Human.txt", sep="\t", header=TRUE, row.names=1))
#CountTable_O_4_4 = as.matrix(read.table("E:/CollaborativeProjects/BertNederland/TissueRNAseq/S_5-7-Human.txt", sep="\t", header=TRUE, row.names=1))
 
#CountTable_AB = as.matrix(read.table("E:/CollaborativeProjects/PfizerProject/DataSets/Count_AB.txt", sep="\t", header=TRUE, row.names=1))

#CountTable_ABC = as.matrix(read.table("E:/CollaborativeProjects/PfizerProject/DataSets/Count_ABC_DeSeq-Alb-C4.txt", sep="\t", header=TRUE, row.names=1))
 
#Check
#CountTable_ABC[1:4,]
#dim(CountTable_ABC)
#CountTable = CountTable_ABC[ ,c(1:6,7:12)]    #A versus B
#CountTable = CountTable_ABC[ ,c(1:6,13:17)]   #A versus C
#CountTable = CountTable_ABC[ ,c(7:12,13:17)]  #B versus C

hek_df= as.matrix(read.table('cell_line_counts.tsv',sep="\t", header= TRUE, row.names= 1 ))
#CountTable [1:4,]
dim(hek_df)
 
#Rownames
#rownames_S_5_7 = as.matrix(read.table("C:/Allt/TCGA-Assembler/HCC/Data from Adil/rownamn_S_5_7.txt"))
# And put rownames on
#rownames(CountTable_S_5_7) = paste(rownames_S_5_7)
 
 
#Tell DEseq which samples are IR and IS
condition = factor(c('X293_E', 'X293_E', 'X293_E', 'X293_F', 'X293_F', 'X293_F', 'X293_H', 'X293_H', 'X293_H', 'X293_T', 'X293_T', 'X293_T', 'Freestyle293', 'Freestyle293', 'Freestyle293', 'HEK293', 'HEK293', 'HEK293'))
#condition = factor(c('293_E', '293_E', '293_E', '293_F', '293_F', '293_F', '293_H', '293_H', '293_H', '293_T', '293_T', '293_T', 'Freestyle293', 'Freestyle293', 'Freestyle293', 'HEK293', 'HEK293', 'HEK293'))


# Create the "count data set" which contains the raw reads and the meta data (conditions)

dds2 <-DESeqDataSetFromMatrix(countData=hek_df,colData=as.data.frame(condition),design=~condition)

#sizeFactors( cds )
#head( counts( cds, normalized=TRUE ) )
 
# NORMALISATION
cds = estimateSizeFactors(dds2 )
 
 
# VARIANCE ESTIMATION
#cds = estimateDispersions( cds ) # Gives warning message: Dispersion fit did not converge.
#cds = estimateDispersions( cds, sharingMode =  "fit-only") # Gives warning message: Dispersion fit did not converge.
cds = estimateDispersions( cds, fitType="local") #This option does not give the warning message: Dispersion fit did not converge.
plotDispEsts( cds )


res= nbinomWaldTest(cds)
####summary(cds)
####cds2 =as.matrix(cds)
####write.table(cds$CountDataSet, file="E:/CollaborativeProjects/PfizerProject/DataSets/DESeq_AB_normalized.txt")
 
# AND FINALY: DIFFERENTIAL EXPRESSION
res_Freestyle293_293_H= nbinomWaldTest(cds, "Freestyle293","293_H")
res_Freestyle293_293_F= nbinomWaldTest(cds, "Freestyle293","293_F")
res_Freestyle293_293_E= nbinomWaldTest(cds, "Freestyle293","293_E")
res_Freestyle293_HEK293= nbinomWaldTest(cds, "Freestyle293","HEK293")
res_Freestyle293_293_T= nbinomWaldTest(cds, "Freestyle293","293_T")
res_293_H_293_F= nbinomWaldTest(cds, "293_H","293_F")
res_293_H_293_E= nbinomWaldTest(cds, "293_H","293_E")
res_293_H_HEK293= nbinomWaldTest(cds, "293_H","HEK293")
res_293_H_293_T= nbinomWaldTest(cds, "293_H","293_T")
res_293_F_293_E= nbinomWaldTest(cds, "293_F","293_E")
res_293_F_HEK293= nbinomWaldTest(cds, "293_F","HEK293")
res_293_F_293_T= nbinomWaldTest(cds, "293_F","293_T")
res_293_E_HEK293= nbinomWaldTest(cds, "293_E","HEK293")
res_293_E_293_T= nbinomWaldTest(cds, "293_E","293_T")
res_HEK293_293_T= nbinomWaldTest(cds, "HEK293","293_T")


#head(res)

write.table(res_Freestyle293_293_H, file="res_Freestyle293_vs_293_H.tsv")
write.table(res_Freestyle293_293_F, file="res_Freestyle293_vs_293_F.tsv")
write.table(res_Freestyle293_293_E, file="res_Freestyle293_vs_293_E.tsv")
write.table(res_Freestyle293_HEK293, file="res_Freestyle293_vs_HEK293.tsv")
write.table(res_Freestyle293_293_T, file="res_Freestyle293_vs_293_T.tsv")
write.table(res_293_H_293_F, file="res_293_H_vs_293_F.tsv")
write.table(res_293_H_293_E, file="res_293_H_vs_293_E.tsv")
write.table(res_293_H_HEK293, file="res_293_H_vs_HEK293.tsv")
write.table(res_293_H_293_T, file="res_293_H_vs_293_T.tsv")
write.table(res_293_F_293_E, file="res_293_F_vs_293_E.tsv")
write.table(res_293_F_HEK293, file="res_293_F_vs_HEK293.tsv")
write.table(res_293_F_293_T, file="res_293_F_vs_293_T.tsv")
write.table(res_293_E_HEK293, file="res_293_E_vs_HEK293.tsv")
write.table(res_293_E_293_T, file="res_293_E_vs_293_T.tsv")
write.table(res_HEK293_293_T, file="res_HEK293_vs_293_T.tsv")

 
## Quality control
plotMA(res_HEK293_293_T)
hist(res_HEK293_293_T$pval, breaks=100, col="skyblue", border="slateblue", main="")

