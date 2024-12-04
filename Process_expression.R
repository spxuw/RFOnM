library(GEOquery) # R 4.2.0
library(limma)
library(AnnotationDbi)
library(DESeq2)
library(illuminaHumanv4.db)

setwd("/code")

# load series and platform data from GEO

gset <- getGEO("GSE47460", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL14550", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)

colData <- data.frame(condition = gset@phenoData@data[["disease state:ch1"]],Age=gset@phenoData@data[["age:ch1"]],
                      Sex=gset@phenoData@data[["Sex:ch1"]],Smoke = gset@phenoData@data[["smoker?:ch1"]])
rownames(colData) <- colnames(ex)
colData$Sex <- factor(colData$Sex)

ex = ex[,!colData$condition=="Interstitial lung disease"]
colData = colData[!colData$condition=="Interstitial lung disease",]
colData$Age = as.numeric(colData$Age)

missing_index <- is.na(colData$Smoke)
colData$Smoke[missing_index] <- sample(colData$Smoke[!missing_index], sum(missing_index), replace = TRUE)

table(colData$condition)

Age = as.numeric(colData$Age)
Sex = as.character(colData$Sex)
Sex[Sex=="1-Male"] = 1
Sex[Sex=="2-Female"] = 2

group = (colData$condition)
group[group=="Chronic Obstructive Lung Disease"] = "COPD"
group = factor(group)

Smoke = colData$Smoke
Smoke[Smoke=="2-Ever (>100)"] = 2
Smoke[Smoke=="1-Current"] = 1
Smoke[Smoke=="3-Never"] = 3
Smoke = factor(Smoke)

Model <- "~0 + group + Sex + Smoke"
Model <- as.formula(Model)

design <- model.matrix(formula(Model))
fit <- lmFit(ex, design)

contrast.matrix <- makeContrasts(groupCOPD-groupControl, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results<-topTable(fit2, coef=1, number=Inf)


gpl <- getGEO("GPL14550")
annotationData <- Table(gpl)

results['Gene'] = annotationData$GENE_SYMBOL[match(rownames(results),annotationData$ID)]
expanded_results <- separate_rows(results, Gene, sep = " /// ")

expanded_results <- expanded_results[!is.na(expanded_results$Gene),]
expanded_results <- expanded_results[expanded_results$Gene != "", ]

expanded_results <- expanded_results[order(expanded_results$Gene, expanded_results$adj.P.Val), ]
expanded_results <- expanded_results[!duplicated(expanded_results$Gene), ]
expanded_results <- as.data.frame(expanded_results)
rownames(expanded_results) <- expanded_results$Gene

write.table(expanded_results,file = "../data/GEO/P_expression.csv",row.names = T,sep = ",")
