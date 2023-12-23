library("minfi")
setwd("D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/DNA_Methylation/Masked_Intensities")
pathOfSamplesAll_27 = list.files(pattern = "\\.idat", recursive = TRUE)

file.copy(pathOfSamplesAll_27,to = "D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/Data/",overwrite = TRUE)
setwd("D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/Data/")

pathOfSamplesAll_27 = list.files(pattern = "ed\\.idat", recursive = FALSE)
rgset_27 = read.metharray(pathOfSamplesAll_27)

################

detP <- detectionP(rgset_27)
rowmean_detP=rowMeans(detP)
filtered_rgset_27=rgset_27[rowmean_detP<0.01,]
detP <- detectionP(filtered_rgset_27)
colmean_detP=colMeans(detP)
filtered_rgset_27=filtered_rgset_27[,colmean_detP<0.01]
MSet_27 <- preprocessRaw(filtered_rgset_27)
methyl=getM(MSet_27)
unmetyl=getUnmeth(MSet_27)
qc <- getQC(MSet_27)
plotQC(qc)
################
ratioSet = ratioConvert(Mset_27, what = "both", keepCN = TRUE)
gset <- mapToGenome(ratioSet)
beta <- getBeta(gset)
Mval=getM(gset)
head(beta)
densityPlot(beta)
min(beta)
max(beta)

min(Mval)
max(Mval)
write.csv(beta,file = "../beta.csv")
write.csv(beta,file = "../Mval.csv")

##############################
###############DMP
##############################
meth_info = read.csv("../meth27-LUAD.csv")
new_meth_info = data.frame(sample_type=meth_info$sample_type,
                           file_name = meth_info$file_name)
new_meth_info$sample_name=substr(new_meth_info$file_name,1,41)
write.csv(new_meth_info , file = "../new_meth_info.csv")

new_meth_info=unique(subset(new_meth_info,select=-c(file_name)))
rownames(new_meth_info)=new_meth_info$sample_name
new_new_meth_info=new_meth_info[colnames(Mval),]
colInfo=data.frame(status=new_new_meth_info$sample_type)
colInfo$status=factor(colInfo$status,levels = c("Solid Tissue Normal","Primary Tumor"))
design <- model.matrix(~status,colInfo)

# fit the actual linear model to the data
library(limma)
fit <- lmFit(Mval, design)
fit <- eBayes(fit)
result=topTable(fit,num=Inf)


###########################
# Reference: https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html

library(minfi)
library(IlluminaHumanMethylation27kmanifest)
library(IlluminaHumanMethylation27kanno.ilmn12.hg19)
library(limma)

setwd("D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/DNA_Methylation/Masked_Intensities")
pathOfSamplesAll_27 = list.files(pattern="\\.idat", recursive=TRUE)
# read.metharray return RGB channel sets, Green and Red, 
# get error because some green and red files are not
# in the same dir, so copy all of them in another folder called Data
file.copy(pathOfSamplesAll_27,to = "D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/Data/", overwrite = TRUE)
setwd("D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/Data/")
pathOfSamplesAll_27 = list.files("D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/Data/", pattern="\\.idat", recursive=FALSE)
rgset_27 <- read.metharray(pathOfSamplesAll_27)
# Error because we have duplicates, names only the last part (green and red) are different
# we need to get only one channel, red or grn
pathOfSamplesAll_27 = list.files("D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/Data/", pattern="ed\\.idat", recursive=FALSE)
# metharrayonly needs only grn or red, it will itself read the other one
rgset_27 <- read.metharray(pathOfSamplesAll_27)

#============ quality control ================================================
#remove props and samples that are not useful
# check if the probs are really working correctly, return P-Value for each prob for each sample
# derP = 27000 (number of props) * 150 (number of samples)

# first filter Props and then samples with high p-value, the order is because the number of props are more
detP <- detectionP(rgset_27)
rowmean_detP=rowMeans(detP)
filtered_rgset_27=rgset_27[rowmean_detP<0.01,]
detP <- detectionP(filtered_rgset_27)
# mean on columns/samples, return 150 means
colmean_detP=colMeans(detP)
# remove samples with p-values > 0.01
filtered_rgset_27=filtered_rgset_27[,colmean_detP<0.01]
MSet_27 <- preprocessRaw(filtered_rgset_27)
methyl=getM(MSet_27)
unmetyl=getUnmeth(MSet_27)
qc <- getQC(MSet_27)
plotQC(qc)

#============ M and U (Methylset) using preprocessRaw ==========================================================
#preprocessRaw dont apply normalization, later we will use another command with normalization
MSet_27 <- preprocessRaw(filtered_rgset_27)
#MSet_27
# apply getM to get methylation data and getUnmeth for unmethyl
methyl = getM(MSet_27)
unmethyl = getUnmeth(MSet_27)

# =============== quality control ====================================
qc <- getQC(MSet_27) # cluster to two groups, Good and Bad
# the ones with the samll methyl and unmethyl value concider as Bad
plotQC(qc)

# =================== calculate Ratios ===========================
#both : M-value and Beta
ratioSet <- ratioConvert(MSet_27, what="both", keepCN=TRUE)
gset <- mapToGenome(ratioSet)
beta <- getBeta(gset)
Mval=getM(gset)
write.csv(Mval, "D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/Mval.csv")
write.csv(beta, "D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/beta.csv")

head(beta)
# peaks in the plot should be around 0 and 1
densityPlot(beta)

# ============ DMP =====================================
# we want to DEG between normal and tumor, need two columns from meth_info: file_names and sample_type
meth_info = read.csv("D:/Uni/Multi Omics integration/GDCdata/TCGA-LUAD/meth27-LUAD.csv")
# file_name column in meth_info is sample names + "_Red.idat"
# sample names = colName(Mval)
# we should substr file_name to be like sample_names
new_meth_info=data.frame(sample_type=meth_info$sample_type,
                         file_name=meth_info$file_name)
new_meth_info$sample_name=substr(new_meth_info$file_name,1,41)
# duplicate rows because for one sample we have two rows(Grn and red)
new_meth_info=unique(subset(new_meth_info,select=-c(file_name)))

# rownames : sample names so we have rows in the same order in both MVal and meth_info
rownames(new_meth_info)=new_meth_info$sample_name

# new_new[colnames(Mval,)] to have rows in the same order as colnames in Mval
new_new_meth_info=new_meth_info[colnames(Mval),]

# status: normal or tumor
colInfo=data.frame(status=new_new_meth_info$sample_type)

# as.factor is like groupby in python
colInfo$status=factor(colInfo$status, levels = c("Solid Tissue Normal", "Primary Tumor"))
# one-hot classes to 0 and 1
design <- model.matrix(~status,colInfo)

# fit the actual linear model to the data
fit <- lmFit(Mval, design)
fit <- eBayes(fit)
#sort based on adjusted p-value
result=topTable(fit,num=Inf)
# here in results, for rows we have prop names
write.csv(result, file = "../result.csv")
