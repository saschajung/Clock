library(biomaRt)
library(sva)
library(plotly)
library(h2o)
library(DESeq2)

setwd("./.")

set.seed(100)

annot <- read.table("./SraRunTable.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
annot <- annot[,c("Run","age","Sex","ethnicity","disease")]
annot <- annot[order(annot$Run),]
colnames(annot) <- c("Sample","Age","Sex","Race","Disease")
annot$Age_Y <- annot$Age
annot$Age_Y <- gsub("YR$","",annot$Age_Y)
annot$Age_Y <- gsub("yr$","",annot$Age_Y)
idx_mos <- which(grepl("mos",annot$Age_Y))
annot$Age_Y[idx_mos] <- as.character(unname(sapply(annot$Age_Y[idx_mos],function(x){
  x <- gsub("mos","",x)
  x <- strsplit(x,"yr|ys")[[1]]
  return(round(as.numeric(x[1]) + as.numeric(x[2])/12))
})))
annot$Age_Y <- as.numeric(annot$Age_Y)

annot$Age_M <- annot$Age
annot$Age_M <- gsub("YR$","",annot$Age_M)
annot$Age_M <- gsub("yr$","",annot$Age_M)
idx_mos <- which(grepl("mos",annot$Age_M))
annot$Age_M[idx_mos] <- as.character(unname(sapply(annot$Age_M[idx_mos],function(x){
  x <- gsub("mos","",x)
  x <- strsplit(x,"yr|ys")[[1]]
  return(as.numeric(x[1]) + as.numeric(x[2])/12)
})))
annot$Age_M <- as.numeric(annot$Age_M)
annot$Age_M <- annot$Age_M*12

annot$Sex[which(annot$Sex == "Female" | annot$Sex == "female")] <- "F"
annot$Sex[which(annot$Sex == "Male" | annot$Sex == "male")] <- "M"

annot$Race[which(annot$Race == "caucasian")] <- "Caucasian"
annot$Race[which(annot$Race == "CaucasianSardinian")] <- "Caucasian"

data <- read.table("./featurecounts.Rmatrix.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data <- data[,c(1,order(colnames(data)[-1])+1)]
data <- data[order(data$Geneid),]
data$Geneid <- gsub("\\.[0-9]+","",data$Geneid)

ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")

ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = data$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)

data <- data[data$Geneid %in% ensTolen_agg$Group.1,]
all(data$Geneid == ensTolen_agg$Group.1)
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]

rapa <- read.table("./LIB049469_featurecounts.Rmatrix.txt",header = TRUE)

rapa$Geneid <- gsub("\\.[0-9]+","",rapa$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = rapa$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)

rapa <- rapa[rapa$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]

rapa <- rapa[order(rapa$Geneid),]
all(rapa$Geneid == ensTolen_agg$Group.1)

iso <- read.table("./LIB050518_featurecounts.Rmatrix.txt",header = TRUE)

iso$Geneid <- gsub("\\.[0-9]+","",iso$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = iso$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)

iso <- iso[iso$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]

iso <- iso[order(iso$Geneid),]
all(iso$Geneid == ensTolen_agg$Group.1)

ipsc <- read.table("./LIB047784_featurecounts.Rmatrix.txt",header = TRUE)

ipsc$Geneid <- gsub("\\.[0-9]+","",ipsc$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = ipsc$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)

ipsc <- ipsc[ipsc$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]

ipsc <- ipsc[order(ipsc$Geneid),]
ipsc <- ipsc[,c(1,6:9)]
all(ipsc$Geneid == ensTolen_agg$Group.1)

data_all_merged <- merge(data,rapa, by = "Geneid")
data_all_merged <- merge(data_all_merged,iso, by = "Geneid")
data_all_merged <- merge(data_all_merged,ipsc, by = "Geneid")
testinp <- data_all_merged[,2:(ncol(data_all_merged))]
rownames(testinp) <- data_all_merged[,1]
dds <- DESeqDataSetFromMatrix(countData = as.matrix(testinp),
                              colData = data.frame(Sample = colnames(data_all_merged)[2:(ncol(data_all_merged))],
                                                   Group = as.factor(c(rep(1,143),c(1,2,2,2,1,2,2,2,1,2,2,2,1,2,2,2),c(1,1,1,1,1,1,1,3),c(3,3,3,3))),
                                                   Batch = as.factor(c(rep("A",143),rep("B",ncol(rapa)-1),rep("C",ncol(iso)-1),rep("D",ncol(ipsc)-1)))),
                              design = ~ Batch)
dds <- estimateSizeFactors( dds, type = "poscounts")
sizeFactors(dds)
dds1 <- varianceStabilizingTransformation(dds, blind = F, fitType = "parametric")
mat <- assay(dds1)
mm <- model.matrix(~Group+0, colData(dds1))
mat <- limma::removeBatchEffect(mat, batch=dds1$Batch, design=mm)
assay(dds1) <- mat

inp_comp <- cbind.data.frame(t(assay(dds1)[,1:143]),annot[,c("Age_Y","Sex","Race","Disease")])
colnames(inp_comp) <- c(rownames(assay(dds1)),"Age","Sex","Race","Disease")
inp_comp <- inp_comp[which(inp_comp$Disease != "HGPS"),]
torem <- apply(inp_comp,2,function(x){all(x==x[1])})
torem[(length(torem)-3):length(torem)] <- FALSE
inp_comp <- inp_comp[,!torem]
tokeep <- apply(inp_comp[,1:ncol(inp_comp)],2,function(x){!any(unname(unlist(x))==0)})
tokeep[(length(tokeep)-3):length(tokeep)] <- TRUE
inp_comp <- inp_comp[,tokeep]

#### Read MBotC ontology data and prepare groups ####

mbotc <- read.table("./gene-SCPassociations.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")

all_procs <- setdiff(unique(mbotc$ProcessName),"Background genes")

genes_to_procs <- lapply(all_procs,function(x){unique(mbotc$Symbol[mbotc$ProcessName == x])})
names(genes_to_procs) <- all_procs

ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")

symbToEns <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = unique(mbotc$Symbol), mart= ensembl)

ensembl_to_procs <- lapply(genes_to_procs,function(x){intersect(colnames(inp_comp),unique(symbToEns$ensembl_gene_id[symbToEns$hgnc_symbol %in% x]))})
ensembl_to_procs <- ensembl_to_procs[unname(which(sapply(ensembl_to_procs,length) >= 3))]

level_to_procs <- sapply(names(ensembl_to_procs),function(x){unique(mbotc$ProcessLevel[mbotc$ProcessName == x])})

levels_to_consider <- c(3)

f <- function(age){if(age<=20){log(age+1)-log(21)}else if(age>=80){log(exp(20/7)+(age-80)/80)}else{(age-20)/(21)}}
f1 <- function(age){-log(1/age-1/101)-log(25)}
h <- function(age){if(age<=70){(f1(71)-f1(70))*(age-70)}else{3*(f1(71)-f1(70))*(age-70)}}
inp_comp$Age_scaled <- sapply(inp_comp$Age,h)

conn <- h2o.init(max_mem_size = "8G",
         min_mem_size = "8G",
         port = 43295)

models <- list()
preds_train <- list()
preds_test <- list()
models_pca <- list()
for(l in levels_to_consider){
  genesInLevel <- unique(do.call("c",ensembl_to_procs[names(level_to_procs)[level_to_procs == l]]))
  inp_comp_level <- inp_comp[,c(genesInLevel,"Age_scaled")]
  inp_comp_level_train <- as.h2o(inp_comp_level)
  for(proc in names(level_to_procs)[level_to_procs == l]){
    print(proc)
    genesToConsider <- intersect(ensembl_to_procs[[proc]],colnames(inp_comp))
    models_pca[[proc]] <- h2o.prcomp(training_frame = inp_comp_level_train,
                                x = genesToConsider,
                                transform = "STANDARDIZE",
                                k = length(genesToConsider),
                                impute_missing = TRUE,
                                seed = 100)
    trans <- predict(models_pca[[proc]],inp_comp_level_train)
    inp <- h2o.cbind(trans,inp_comp_level_train[,which(colnames(inp_comp_level_train) == "Age_scaled")])
    models[[proc]] <- h2o.glm(y = "Age_scaled",
                              training_frame = inp,
                              nfolds = 5,
                              fold_assignment = "Random",
                              family = "AUTO",
                              link = "identity",
                              seed = 100,
                              alpha = 0,
                              lambda_search = T,
                              standardize = T,
                              keep_cross_validation_predictions = F,
                              keep_cross_validation_fold_assignment = F,
                              )

    modpred <- predict(models[[proc]],inp)
    preds_train[[proc]] <- as.numeric(as.data.frame(modpred)$predict)

    rm(inp)
    rm(modpred)
    rm(trans)
    #gc()
    #gc()
    #h2o:::.h2o.garbageCollect()
    #h2o:::.h2o.garbageCollect()
    #h2o:::.h2o.garbageCollect()
  }
}

model_r2 <- as.data.frame(sapply(models,h2o.r2))
model_mae <- as.data.frame(sapply(models,h2o.mae))

mods_to_include <- unlist(lapply(models,function(x){!all(x@model[["coefficients_table"]][["coefficients"]][-1] == 0)}))

models <- models[names(mods_to_include)[which(mods_to_include == T)]]
preds_train <- preds_train[names(mods_to_include)[which(mods_to_include == T)]]

preds_train_sig <- preds_train[names(preds_train) %in%rownames(model_r2)[model_r2$`sapply(models, h2o.r2)` >= 0.0]]
preds_mat <- do.call("cbind",preds_train_sig)
rownames(preds_mat) <- rownames(inp_comp)

preds_mat_h2o <- as.h2o(cbind(preds_mat,inp_comp[,"Age_scaled",drop = FALSE]))

comb_model <- h2o.glm(y = "Age_scaled",
            training_frame = preds_mat_h2o,
            nfolds = 5,
            fold_assignment = "Random",
            lambda_search = T,
            standardize = T,
            alpha = 1,
            seed = 101,
        solver = "IRLSM"
            )

mod_var_imp <- as.data.frame(h2o.varimp(comb_model))

#Remove processes that have less then 10 genes

selectedProcs <- mod_var_imp$variable[which(mod_var_imp$relative_importance > 0)]
numGenes <- sapply(selectedProcs,function(x){length(ensembl_to_procs[[x]])})
priorProcs <- names(numGenes[which(numGenes >= 10)])




