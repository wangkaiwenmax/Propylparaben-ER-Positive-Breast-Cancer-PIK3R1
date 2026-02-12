####01差异分析####

library(edgeR)
library(pheatmap)

inputFile="input.txt"       
pFilter=0.05                     
logFCfilter=0.585                  
conFile="sample1.txt"            
treatFile="sample2.txt"          

rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

sample1=read.table(conFile,sep="\t",header=F,check.names=F)
sample2=read.table(treatFile,sep="\t",header=F,check.names=F)
conData=data[,as.vector(sample1[,1])]
treatData=data[,as.vector(sample2[,1])]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

group=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("con","treat"))
ordered_tags <- topTags(et, n=100000)
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$PValue)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

diffSig = diff[(diff$PValue < pFilter & (diff$logFC>logFCfilter | diff$logFC<(-logFCfilter))),]
normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)  
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffGeneExp.txt",sep="\t",quote=F,col.names=F)        
write.table(diff,file="all.xls",sep="\t",quote=F)

diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut, file="diff.txt",sep="\t",quote=F,col.names=F)

#####02  GO分析####

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

pvalueFilter <- 0.05        
p.adjustFilter <- 0.05     

rt <- read.table("input.txt", header=F, sep="\t", check.names=F)  
genes <- unique(as.vector(rt[,1]))  
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
gene <- entrezIDs[entrezIDs != "NA"]  
kk <- enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO <- as.data.frame(kk)
GO <- GO[(GO$pvalue < pvalueFilter & GO$p.adjust < p.adjustFilter),]  
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#####03  KEGG富集分析####

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

pvalueFilter <- 0.05        
p.adjustFilter <- 0.05      

rt <- read.table("input.txt", header=F, sep="\t", check.names=F)
genes <- unique(as.vector(rt[,1]))  
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt <- data.frame(genes, entrezID=entrezIDs)
gene <- entrezIDs[entrezIDs != "NA"] 
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG <- as.data.frame(kk)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) paste(rt$genes[match(strsplit(x,"/")[[1]], as.character(rt$entrezID))], collapse="/")))
KEGG <- KEGG[(KEGG$pvalue < pvalueFilter & KEGG$p.adjust < p.adjustFilter),]  
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

####04 LASSO####

set.seed(123)
library(glmnet)                  
inputFile="input.txt"      

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)


####05 SVM####

library(e1071)
library(kernlab)
library(caret)

set.seed(123)
inputFile="input.txt"       

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

Profile=rfe(x=data,
            y=as.numeric(as.factor(group)),
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")
pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()
featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)

####06 RF####

set.seed(123)
library(randomForest)
library(ggplot2)
inputFile="input.txt"     
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)
importance=importance(x=rf2)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>5])  
write.table(rfGenes, file="RF.gene.txt", sep="\t", quote=F, col.names=F, row.names=F)

####07 AUC####

library(glmnet)
library(pROC)

expFile <- "merge.normalize.txt"      
geneFile <- "interFeatureGenes.txt"   
rt <- read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
y <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(rt))
y <- ifelse(y == "Control", 0, 1)
geneRT <- read.table(geneFile, header=F, sep="\t", check.names=F)
bioCol <- c("#F05C3BFF","#5C88DAFF","#5CB85CFF", "#EEA236FF", "#9632B8FF", "#17BECFFF", "#BCBD22FF")
bioCol <- bioCol[1:nrow(geneRT)]
aucText <- c()
k <- 0
for (x in as.vector(geneRT[,1])) {
  k <- k + 1
  roc1 <- roc(y, as.numeric(rt[x,])) 
  if (k == 1) {
    pdf(file="ROC.genes.pdf", width=5.5, height=5)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
    aucText <- c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  } else {
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE, lwd=3)
    aucText <- c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}
legend("bottomright", aucText, lwd=3, bty="n", cex=1, col=bioCol[1:(ncol(rt)-1)])
dev.off()
rt <- rt[as.vector(geneRT[,1]), ]
rt <- as.data.frame(t(rt))
logit <- glm(y ~ ., family=binomial(link='logit'), data=rt)
pred <- predict(logit, newx=rt) 
roc1 <- roc(y, as.numeric(pred))  
ci1 <- ci.auc(roc1, method="bootstrap")  
ciVec <- as.numeric(ci1)
pdf(file="ROC.model.pdf", width=5.5, height=5)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Model", lwd=3)
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()


####07 MR####

library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
exposureFile="exposure.F.csv"        
outcomeID="ieu-a-1127"        
exposure_dat=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)
outcomeData=extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcomeID)
write.csv(outcomeData, file="outcome.csv", row.names=F) 
outcomeData$outcome=outcomeName
dat=harmonise_data(exposure_dat, outcomeData)
dat=dat[dat$pval.outcome>1e-5,]
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)
mrResult=mr(dat)
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

res_single=mr_singlesnp(dat)      
pdf(file="pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()



#####08 ssGSEA####

library(GSVA)
library(limma)
library(GSEABase)
library(pROC)
library(dplyr)
library(ggplot2)
library(reshape2)

expFile <- "merge.normalize.txt"  
geneFile <- "interFeatureGenes.txt"  
rt <- read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
geneRT <- read.table(geneFile, header=F, sep="\t", check.names=F)

immuneScore <- function(expFile, gmtFile, project) {
  rt <- read.table(expFile, header=T, sep="\t", check.names=F)
  rt <- as.matrix(rt)
  rownames(rt) <- rt[,1]
  exp <- rt[,2:ncol(rt)]
  geneSet <- getGmt(gmtFile, geneIdType=SymbolIdentifier())
  ssgseaScore <- gsva(exp, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
  normalize <- function(x) { return((x - min(x)) / (max(x) - min(x))) }
  ssgseaOut <- normalize(ssgseaScore)
  ssgseaOut <- rbind(id=colnames(ssgseaOut), ssgseaOut)
  write.table(ssgseaOut, file=paste0(project, ".score.txt"), sep="\t", quote=F, col.names=F)
}

immuneScore(expFile="merge.normalize.txt", gmtFile="cellMarker.gmt", project="EXP")
scoreCor <- function(riskFile, scoreFile, project) {
  data <- read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
  data <- t(data)
  risk <- read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  sameSample <- intersect(row.names(data), row.names(risk))
  data <- data[sameSample, , drop=F]
  risk <- risk[sameSample, , drop=F]
  rt <- cbind(data, risk[, c("PIK3R1", "group")])
  rt <- rt[, -ncol(rt)]

  immCell <- colnames(rt)
  rt1 <- rt[, immCell]
  data <- melt(rt1, id.vars=c("group"))
  colnames(data) <- c("Group", "Type", "Score")
  data$Group <- factor(data$Group, levels=c("PIK3R1_Low", "PIK3R1_High"))
  p <- ggboxplot(data, x="Type", y="Score", color="black", width=0.5, palette="npg", outlier.shape=NA, fill="Group") +
    theme_test() +
    stat_compare_means(aes(group=Group), method="wilcox.test", label="p.signif") +
    theme(text = element_text(size=10), legend.position="top", axis.text.x=element_text(angle=45, hjust=1))
  pdf(file=paste0(project, ".immCell.pdf"), width=10, height=7)
  print(p)
  dev.off()
}
scoreCor(riskFile="Group.txt", scoreFile="EXP.score.txt", project="EXP")
gene <- c("PIK3R1")
rt <- read.table("merge.normalize.txt", header=T, sep="\t", check.names=F)
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
data <- as.data.frame(t(exp))
data <- data[rowMeans(data) > 0,]
group <- data.frame(id=colnames(data), PIK3R1=data[gene,])
group$group <- ifelse(group$PIK3R1 > median(group$PIK3R1), "PIK3R1_High", "PIK3R1_Low")
write.table(group, "Group.txt", row.names=F, quote=F, sep="\t")
for (gene in gene) {
  cor_data <- cbind(exist_cib, gene=data[,gene])
  result <- data.frame(ColumnName=character(), Correlation=numeric(), PValue=numeric(), stringsAsFactors=FALSE)
  for (col in colnames(exist_cib)) {
    cor_test <- cor.test(exist_cib[[col]], data[,gene], method="pearson")
    result <- rbind(result, data.frame(Cell=col, Correlation=cor_test$estimate, PValue=cor_test$p.value))
  }
write.csv(result, paste0(gene, "_correlation_results.csv"), row.names=FALSE)
}





