r_imports= """
#library(vegan)
library(DESeq2)
#library(edgeR)
#library(reshape2)
#require (labdsv)

"""
runMetaMDS= """runMetaMDS = function(data, analysis_type){
    
    #pdf(paste("r_NMDS.pdf",sep=""))
    
    #title= paste(analysis_type, "NMDS Ordination")
    #env= metadata[rownames(data),]
    #names= sapply(env$Triclosan, toString )
    #print(rownames(data)) 
    #print(colnames(data)) 
    if (analysis_type == "16S" | analysis_type == "18S" | analysis_type == "phylum")
    {
        dist= "bray"
    }
    else {
        dist= "horn"
    }

    nmds = metaMDS(data, type="n", distance= dist, trymax=500, stress= .2)
    #sites <- scores(nmds, display = "sites")
    #print(sites)

    #plot(nmds, main= "18S data")
    #points(sites, col = "blue", pch = 16, cex = 1.2)
    #ordipointlabel(nmds, display='sites')
    #text(nmds, labels= names)
    #legend("bottomright", paste("Stress:", sprintf("%.2f",nmds$stress)),pch=1 ,  inset = .02)
    #dev.off()
    return(nmds) 
    #plot(nmds, display = "sites")
    #n99 = scores(nmds)
    #nmds$stress
    #output = as.data.frame(merge(n99,mdata, by.x = "row.names", by.y = "row.names"))
}
"""

mdsScores= """mdsScores= function(nmds)
{
return(scores(nmds))
}
"""

#!!! merge data with metadata for some arbitrary reason


buildDiversityIndices= """buildDiversityIndices <- function(otu_data)
{

#adonis(md[val.cols] ~ md$Triclosan, permutations=1000, method="bray")
chao1 = as.data.frame(t(estimateR(otu_data)))$S.chao1  #chao1 and ACEs
shannon =  diversity(otu_data)                   # shannon
richness  =  specnumber(otu_data)
pielou  =  shannon/log(richness)                         #Pielou's
#a5 = cbind(pielou,shannon, merged_data)
return(list(chao1, shannon, richness, pielou))
}
"""

drawLogAxis= """drawlogaxis <- function(side,range)
{
	par(tck=0.02)
#	d <- log(range,10)
	d <- range
	mlog <- floor(min(d))
	Mlog <- ceiling(max(d))
	#SeqLog <- c(mlog:Mlog)
	SeqLog = c(-1.0,-0.5, 0,0.5,1,1.5,2,2.5,3)
    #expression(paste("10"^"th"))
    #labels= parse(text = sprintf("10^%.1f", SeqLog)) 
    #Nlog <- (Mlog-mlog)+1
	labels= c("0","0.316","1","3.16","10","31.6","100","316","1000")
    axis(side, at=SeqLog, labels=labels)
	#ats <- log(seq(from=2,to=9,by=1),10)
	#mod <- NULL
	#for(i in SeqLog)
	#{
	#	mod <- c(mod,rep(i,length(ats)))
	#}
	#ats <- rep(ats,Nlog)
	#ats <- ats+mod
	#par(tck=0.02/3)
	#axis(side,at=ats,labels=NA)
}
"""



calCommunityResponse= """

calCommunityResponse = function(otu_data)

{

va = as.matrix(metaMDSdist(otu_data, distance="bray", trymax=2))
va1= melt(va)
va2 = va1[order(va1$Var1),]
return(va2)

}
"""


plotDiversityIndex= """plotDiversityIndex <- function (type_short, type_long, merged_data, index_data, analysis_type)
{
title= paste(analysis_type, "diversity")
outF= paste('triclosan',type_short,'log10.pdf', sep="_")
pdf(outF, width= 7, height= 7)
plot(log10(merged_data$Triclosan+.1),index_data, xlab = "Triclosan concentration (log10) [nM]", ylab= type_long, main= title, xaxt='n', pch=20)
#text(log10(md$Triclosan+1), index_data, labels= merged_data$Triclosan,  adj = c( 1, 1 ), pos= rep(1, length(merged_data$Triclosan)))
drawlogaxis(1, c(0,3))
abline(lm(index_data ~ log10(merged_data$Triclosan +0.1)))
dev.off()
}
"""



plotRarefaction= """plotRarefaction = function(data, analysis_type, analysis_dir)
{
sac= specaccum(data, method= "rarefaction", gamma="jack1")
outF= paste(analysis_dir, "rarefaction.pdf", sep="/")
pdf(outF, width= 7, height= 7)
plot(sac, ci.type="polygon", ci.col="grey",xlab="Samples", ylab="OTU Richness")
dev.off()
}
"""



getDoseResponseModel= """getModel <- function(df, x,y,model, log=TRUE)
{
    
    if(log == TRUE){
        y= log10(y+.1)
    } 

    if(model == "ll"){
        
        df_m=  drm(x ~ y, data = df, fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    }
    
    else if (model == "weibull1"){
        df_m=  drm(x ~ y, data = df, fct = W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    
    }

    else if (model == "weibull2"){
        df_m=  drm(x ~ y, data = df, fct = W2.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    
    }

    df_m
}
"""



getFit= """getFit <- function(df, phylum){

models=c('weibull1',"weibull2", "ll")
attach(df)

for(i in 1:length(models))
{
    model_name= models[i]
    print(model_name)
    result= tryCatch({
    model= getModel(df, phylum, triclosan, model_name)
}, warning = function(war) {
 
  # warning handler picks up where error was generated
  print(paste("MY_WARNING:  ",war))

}, error = function(err) {
    #print(paste("MY_ERROR:  ",err))
    model= getModel(df, phylum, triclosan, model_name, log= FALSE)
})

    p_value= modelFit(model)$p[2]
    p_value    
    #if(p_value > 0.05){
    #    print(paste(model_name, p_value))
    #}
    #}
}
}
"""


testTaxa= """testTaxa <- function(df_otu, metadata){
cols= colnames(df_otu)
for(i in 1:length(cols)){
    phylum= cols[i]
    print(phylum)
    df_phylum= data.frame(df_otu[,phylum],  metadata[rownames(df_otu),]["Triclosan"])
    colnames(df_phylum)= c(phylum,"triclosan")
    #getFit(df_phylum, phylum)

}
}
"""


testDeseq2 = """testDeseq2 <- function(Y, group, filter=FALSE, returnEffect=FALSE){



  dds2 <-DESeqDataSetFromMatrix(countData=Y,colData=as.data.frame(group),design=~group)

  dds2 <- estimateSizeFactors(dds2)
  dds2 <- estimateDispersions(dds2)
  dds2 <- nbinomWaldTest(dds2)
  if (filter) {
    res <- results(dds2)
  } else {
    res <- results(dds2,cooksCutoff=FALSE, independentFiltering=FALSE)
  }

  if ( returnEffect) {
    return(list(res[,5],2^res[,2]))
  } else {
    return(res[,5])
  }
}
"""



normalizeDeseq2="""normalizeDeseq2 <- function(Y, groups)
{
  
    
    dds2 <-DESeqDataSetFromMatrix(countData=Y, colData=as.data.frame(groups), design=~1)
    sizeFactors <- estimateSizeFactors(dds2)
    
    return(counts(sizeFactors, normalized=TRUE))

}

"""


correlateDeseq2="""correlateDeseq2 <- function(otu_data, metadata)
{
    ### index column is messing with the function
    dds= DESeqDataSetFromMatrix(otu_data, metadata, ~Triclosan)
    dds= DESeq(dds)
    corr= results(dds)
    return(list(corr$stat, corr$pvalue, corr$padj))

}
"""



testEdgeR= """testEdgeR <- function(Y, group,returnEffect=FALSE){

  dge <- DGEList(counts=Y,group=group)

  mod <- model.matrix(~group)
  dge <- calcNormFactors(dge,method="TMM")
  dge <- estimateGLMCommonDisp(dge,mod)
  dge <- estimateGLMTrendedDisp(dge,mod)
  dge <- estimateGLMTagwiseDisp(dge,mod)
  fitdge <- glmFit(dge,mod)
  lrt <- glmLRT(fitdge,coef=2)
  if ( returnEffect) {
    return(list(lrt$table[,4],exp(fitdge$coefficients[,2])))
  } else {
    return(lrt$table[,4])
  }
}
"""









