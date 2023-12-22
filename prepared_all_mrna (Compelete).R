library(parallel)
library(limma)
library(xlsx)
setwd("D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/")
data_path="D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/"
pathOfSamplesAll=list.files(data_path, recursive = T, pattern = paste0("*.",".tsv"), full.names = T)
pathsSplitted = split(pathOfSamplesAll, cut(seq_along(pathOfSamplesAll), detectCores()-1))
names(pathsSplitted) = seq(detectCores()-1)
first_sampleExpression=read.delim(pathOfSamplesAll[1] , comment.char = '#', check.names = F,sep="\t")
summary_genes=data.frame(gene_id=first_sampleExpression$gene_id,
                         gene_name=first_sampleExpression$gene_name,
                         gene_type=first_sampleExpression$gene_type)

expressionMatrix =
  mclapply(pathsSplitted,
           function(paths)
           {
             tempExpressionMatrix = summary_genes
             for(j in paths)
             {
               
               temp = strsplit2(j, '/')
               fileName = temp[length(temp)]
               submitterID=strsplit2(fileName,'[.]')[1]
               sampleExpression=read.delim(j, comment.char = '#', check.names = F,sep="\t")[,c("gene_id","unstranded")]
               colnames(sampleExpression)=c("gene_id",submitterID)
               rownames(sampleExpression)=sampleExpression$gene_id
               sampleExpression=sampleExpression[tempExpressionMatrix$gene_id,]
               tempExpressionMatrix=cbind(tempExpressionMatrix,sampleExpression[submitterID])
             }
             tempExpressionMatrix$ENSEMBL_ID=NULL
             tempExpressionMatrix
           },
           mc.cores = 1)
expressionMatrix = Reduce(cbind, expressionMatrix)
write.csv(expressionMatrix,"D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/AllExpressionMatrix.csv")
expressionMatrix=read.csv("D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification//AllExpressionMatrix.csv")
#prepared_expressionMatrix=aggregate(.~ gene_name, data=subset(expressionMatrix, select = -c(gene_id,gene_type)), FUN=max)
#write.csv(prepared_expressionMatrix,"/media/lbb-admin/EXTERNAL_USB/thesis/TCGA/Trascriptomic/AllProjects/prepared_expressionMatrix.csv")



expressionMatrix=expressionMatrix[expressionMatrix$gene_type=="protein_coding",]
dim(expressionMatrix)

expressionMatrix$gene_name=gsub("\\..*", "",expressionMatrix$gene_name)
length(expressionMatrix$gene_name)
length(unique(expressionMatrix$gene_name))
prepared_expressionMatrix=aggregate(.~gene_name, data = subset(expressionMatrix, select = -c(gene_id, gene_type)), FUN = max)
length(prepared_expressionMatrix$gene_name)
length(unique(prepared_expressionMatrix$gene_name))
sum(is.na(prepared_expressionMatrix))
rownames(prepared_expressionMatrix)=prepared_expressionMatrix$gene_name
gene_names=rownames(prepared_expressionMatrix)
prepared_expressionMatrix=prepared_expressionMatrix[,!grepl("gene",colnames(prepared_expressionMatrix))]
prepared_expressionMatrix=apply(prepared_expressionMatrix, 2, as.numeric)
prepared_expressionMatrix=subset(prepared_expressionMatrix,select=-c(X))
rownames(prepared_expressionMatrix)=gene_names
RowMean=rowMeans(prepared_expressionMatrix)>1
prepared_expressionMatrix=prepared_expressionMatrix[RowMean,]
dim(prepared_expressionMatrix)
write.csv(prepared_expressionMatrix,"D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/prepared_expressionMatrix.csv")
prepared_expressionMatrix=read.csv("D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/prepared_expressionMatrix.csv")
head(colnames(prepared_expressionMatrix))
rownames(prepared_expressionMatrix)=prepared_expressionMatrix$X
prepared_expressionMatrix=subset(prepared_expressionMatrix,select=-c(X))
colnames(prepared_expressionMatrix)=gsub("\\.","-",gsub("X","",colnames(prepared_expressionMatrix)))

info=read.csv("D:/Uni/Multi Omics integration/geo-KIRC.csv")
info$file_name=substr(info$file_name,1,36)
table(info$sample_type)
sum(colnames(prepared_expressionMatrix) %in% info$file_name)

rownames(info)=info$file_name
colnames(prepared_expressionMatrix)=info[colnames(prepared_expressionMatrix),"cases"]
rownames(info)=info$cases

#If data have few or zero normal samples run this section code(4 # Lines)

sum(rowSums(prepared_expressionMatrix==0)>(0.1*dim(prepared_expressionMatrix)[2]))

Mymean<-rowMeans(prepared_expressionMatrix[rowSums(prepared_expressionMatrix==0)>(0.1*dim(prepared_expressionMatrix)[2]),])
hist(Mymean)

prepared_expressionMatrix=prepared_expressionMatrix[rowSums(prepared_expressionMatrix==0)<(0.1*dim(prepared_expressionMatrix)[2]),]

tumor_exp=prepared_expressionMatrix[,info[info$sample_type=="Primary Tumor","cases"]]
write.csv(tumor_exp,"D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/tumor_exp.csv")


Normal_exp=prepared_expressionMatrix[,info[info$sample_type=="Solid Tissue Normal","cases"]]
write.csv(Normal_exp,"D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/Normal_exp.csv")

#normalized_tumor_exp= log(tumor_exp+1)
#dim(prepared_expressionMatrix)
#write.csv(normalized_tumor_exp,"G:/TCGA_BRCA project(MultiOmics)/GDCdata/TCGA-BRCA/Transcriptome_Profiling/Gene_Expression_Quantification/normalized_tumor_exp.csv")
#dim()

#normalized_tumor_exp=log(tumor_exp+1)
#sds=apply(normalized_tumor_exp,1,sd)
#qs=quantile(sds)
#sum(sds<qs[2])
#filtered_tumor_exp=tumor_exp[sds>qs[2],]
#filtered_Normalized_tumor_exp=normalized_tumor_exp[sds>qs[2],]




library(DESeq2)

Normal_exp = read.csv("D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/Normal_exp.csv")
tumor_exp = read.csv("D:/Uni/Multi Omics integration/GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/tumor_exp.csv")


rownames(tumor_exp)=tumor_exp$X
rownames(Normal_exp)=Normal_exp$X
tumor_exp=subset(tumor_exp,select=-c(X))
Normal_exp=subset(Normal_exp,select=-c(X))
CountData=cbind(Normal_exp,tumor_exp)
ColInfo=data.frame(status=c(rep("normal",dim(Normal_exp)[2]),rep("Tumor",dim(tumor_exp)[2])))
ColInfo$status=factor(ColInfo$status, levels = c("normal","Tumor"))
dds <- DESeqDataSetFromMatrix(countData = CountData,colData = ColInfo,design= ~ status)
dds <- DESeq(dds)
res <- results(dds)
res_df=data.frame(res)
final_res=res_df[abs(res_df$log2FoldChange)>2 & res_df$padj <0.05,]










