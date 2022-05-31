library(biomaRt)
library(sva)
library(plotly)
library(h2o)
library(DESeq2)

set.seed(100)

#### Data normalization and transformation ####
#Incorporate WT&BFP data with train data, normalize, and freeze normalization arguments
testinp <- data_merged[,2:(ncol(data_merged))]
rownames(testinp) <- data_merged[,1]

dds <- DESeqDataSetFromMatrix(countData = as.matrix(testinp),
                              colData = data.frame(Sample = colnames(testinp),
                                                   Batch = as.factor(c(rep("A",(ncol(testinp)-length(idx_wt)-length(idx_bfp)-length(idx_test))),asl_ctrl_batches,asl_batches)),
                                                   Age = (c(train_metadata$Age_Y[1:(ncol(testinp)-length(idx_wt)-length(idx_bfp)-length(idx_test))],asl_ctrl_ages,rep(NA,length(idx_test)))),
                                                   Dox = as.factor(c(rep("F",(ncol(testinp)-length(idx_wt)-length(idx_bfp)-length(idx_test))),rep("F",length(idx_wt)),rep("T",length(idx_bfp)),rep("T",length(idx_test)))),
                                                   Group = as.factor(c(rep("NTg",133),asl_ctrl_gene,asl_gene)),
                                                   Passage = scale(c(rep(0,133),asl_ctrl_passage,asl_passage))),
                              design = ~ Group + Batch)

dds <- estimateSizeFactors( dds, type = "ratio")
sizeFactors(dds)
dds1 <- varianceStabilizingTransformation(dds, blind = T, fitType = "parametric")
mat <- assay(dds1)


mm <- model.matrix(~1, colData(dds))
mat <- limma::removeBatchEffect(mat, batch=dds$Batch, batch2 = NULL, design=mm)

inp_comp <- cbind.data.frame(t(mat[,1:(ncol(mat)-length(idx_test))]),dds$Age[1:(ncol(mat)-length(idx_test))])
colnames(inp_comp) <- c(rownames(mat),"Age")
torem <- apply(inp_comp,2,function(x){abs(var(x))<10^-6})
torem[length(torem)] <- FALSE
inp_comp <- inp_comp[,!torem]
tokeep <- apply(inp_comp[,1:ncol(inp_comp)],2,function(x){!any(unname(unlist(x))==0)})
tokeep[length(tokeep)] <- TRUE
inp_comp <- inp_comp[,tokeep]

#### Read MBotC ontology data and prepare groups ####

mbotc <- read.table("./Data/gene-SCPassociations.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")

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
h_inv <- function(x){
  if(x <= 0){
    return((-x + 70*log(7070/31) - 70*log(7171/30))/(log(7070/31) - log(7171/30)))
  }else{
    return((-x + 210*log(7070/31) - 210*log(7171/30))/(3*(log(7070/31) - log(7171/30))))
  }
}
inp_comp$Age_scaled <- sapply(inp_comp$Age,h)

conn <- h2o.init(max_mem_size = "8G",
                 min_mem_size = "8G",
                 port = 43295)

obs_weights <- rep(1,nrow(inp_comp))
batches <- dds$Batch[dds$Sample %in% rownames(inp_comp)]
for(l in levels(batches)){
  obs_weights[which(batches == l)] <- 1/length(which(batches == l))
}

inp_comp <- cbind(inp_comp,obs_weights)

models <- list()
preds_train <- list()
preds_test <- list()
models_pca <- list()
for(l in levels_to_consider){
  genesInLevel <- unique(do.call("c",ensembl_to_procs[names(level_to_procs)[level_to_procs == l]]))
  inp_comp_level <- inp_comp[,c(genesInLevel,"Age_scaled","obs_weights")]
  inp_comp_level_train <- as.h2o(inp_comp_level)
  for(proc in names(level_to_procs)[level_to_procs == l]){
    print(proc)
    genesToConsider <- intersect(ensembl_to_procs[[proc]],colnames(inp_comp))
  
    models[[proc]] <- h2o.glm(y = "Age_scaled",
                              x = genesToConsider,
                              training_frame = inp_comp_level_train,
                              nfolds = 5,
                              fold_assignment = "Random",
                              family = "AUTO",
                              link = "identity",
                              seed = 100,
                              weights_column = "obs_weights",
                              alpha = 0,
                              lambda_search = T,
                              standardize = T,
                              keep_cross_validation_predictions = F,
                              keep_cross_validation_fold_assignment = F,
    )
    
    modpred <- predict(models[[proc]],inp_comp_level_train)
    preds_train[[proc]] <- as.numeric(as.data.frame(modpred)$predict)
    
    #rm(inp)
    #rm(modpred)
    #rm(trans)
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
                      seed = 100,
                      solver = "IRLSM"
)

comb_model
mod_var_imp <- as.data.frame(h2o.varimp(comb_model))
#View(mod_var_imp)

final_pred_new <- h2o.predict(comb_model,preds_mat_h2o)
final_pred_new <- as.data.frame(final_pred_new)
rownames(final_pred_new) <- rownames(inp_comp)
final_pred_new$predict <- sapply(final_pred_new$predict,h_inv)
final_pred_new_train <- final_pred_new
inbag_pred <- final_pred_new$predict

pred_diff <- inp_comp$Age-inbag_pred
hist(pred_diff,100)

plot_ly(type = "scatter") %>% add_trace(x = inp_comp$Age, y = inbag_pred,color = dds$Group[dds$Sample %in% rownames(inp_comp)])

#### New Predictions ####
inp_comp_rapa <- t(mat[,(ncol(mat)-length(idx_test)+1):ncol(mat)])
inp_comp_rapa <- inp_comp_rapa[,which(colnames(inp_comp_rapa) %in% colnames(inp_comp))]
all(colnames(inp_comp_rapa) == colnames(inp_comp)[1:(ncol(inp_comp)-3)])

inp_comp_rapa_h2o <- as.h2o(inp_comp_rapa)

preds_new <- list()
for(n in names(models)){
  print(n)
  #pca_trans <- h2o.predict(models_pca[[n]],inp_comp_rapa_h2o)
  #preds_new[[n]] <- h2o.predict(models[[n]],pca_trans)
  preds_new[[n]] <- h2o.predict(models[[n]],inp_comp_rapa_h2o)
  #h2o.rm(pca_trans)
}

preds_new_mat <- do.call("cbind",lapply(preds_new,as.data.frame))
colnames(preds_new_mat) <- names(preds_new)

final_pred_new <- h2o.predict(comb_model,as.h2o(preds_new_mat))
final_pred_new <- as.data.frame(final_pred_new)
rownames(final_pred_new) <- rownames(inp_comp_rapa)
final_pred_new$predict <- sapply(final_pred_new$predict,h_inv)
final_pred_new_new <- final_pred_new

final_pred <- rbind(final_pred_new_train,final_pred_new_new)

m <- merge(meta_comb_new,final_pred,by.x = "SampleID",by.y = 0)

m$Agediff <- m$predict-m$Line_age
View(m)
write.table(m,file = "~/Documents/LIB054839_Preds.txt", sep = "\t", quote = F, row.names = F, col.names = T)

