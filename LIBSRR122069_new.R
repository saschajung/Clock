library(biomaRt)
library(sva)
library(plotly)
library(h2o)
library(DESeq2)
library(NOISeq)

set.seed(100)

h2o.removeAll()
h2o.removeAll()
h2o.shutdown(F)

f1 <- function(age){-log(1/age-1/101)-log(25)}
h <- function(age){if(age<=70){(f1(71)-f1(70))*(age-70)}else{3*(f1(71)-f1(70))*(age-70)}}
h_inv <- function(x){
  if(x <= 0){
    return((-x + 70*log(7070/31) - 70*log(7171/30))/(log(7070/31) - log(7171/30)))
  }else{
    return((-x + 210*log(7070/31) - 210*log(7171/30))/(3*(log(7070/31) - log(7171/30))))
  }
}

#### Data normalization and transformation ####
#Incorporate WT&BFP data with train data, normalize, and freeze normalization arguments
testinp <- data_merged[,2:(ncol(data_merged))]
rownames(testinp) <- data_merged[,1]

dds <- DESeqDataSetFromMatrix(countData = as.matrix(testinp),
                              colData = data.frame(Sample = colnames(testinp),
                                                   Batch = as.factor(c(rep("A",(ncol(testinp)-length(idx_wt)-length(idx_bfp)-length(idx_test))),asl_ctrl_batches,asl_batches)),
                                                   Age = sapply(c(train_metadata$Age_Y[1:(ncol(testinp)-length(idx_wt)-length(idx_bfp)-length(idx_test))],asl_ctrl_ages,rep(-1,length(idx_test))),h),
                                                   Dox = as.factor(c(rep("F",(ncol(testinp)-length(idx_wt)-length(idx_bfp)-length(idx_test))),rep("F",length(idx_wt)),rep("T",length(idx_bfp)),rep("T",length(idx_test)))),
                                                   Group = as.factor(c(rep("NTg",133),asl_ctrl_gene,asl_gene)),
                                                   Passage = scale(c(rep(0,133),asl_ctrl_passage,asl_passage))),
                              design = ~ Group + Batch)

mat <- DESeq2::counts(dds,normalized = F,replaced = F)
mat <- tmm(mat,refColumn = 4)
#mat <- log2(mat+1)
mat[1:5,1:5]

fleischersamps <- intersect(train_metadata$Sample,colnames(dds))
trainSamps <- c(fleischersamps,colnames(samples_new_df)[c(idx_wt)])
testSamps <- setdiff(colnames(dds),trainSamps)

cvs <- apply(mat,1,function(x){sd(x)/mean(x)})
cvs <- cvs[!is.na(cvs) & !is.infinite(cvs)]
cvs <- sort(cvs)

controlgenes <- names(cvs[1:10000])

av <- colData(dds)[fleischersamps,"Age"]
sampcors <- apply(mat[,fleischersamps],1,function(x){cor(x,av)})
sampcors <- abs(sampcors)
sampcors[is.na(sampcors)] <- 0
hist(sampcors)

idx <- which(sampcors < 0.05)

controlgenes <- intersect(names(idx),controlgenes)

control <- rep(0,nrow(mat))
names(control) <- rownames(mat)
control[controlgenes] <- 1
control <- control[rownames(mat)]

control_idx <- which(control == 1)

mat_ruv <- RUVnormalize::naiveRandRUV(Y = t(mat),
                                      cIdx = control_idx,
                                      k = 1)
mat_ruv <- t(mat_ruv)
mat.pca <- svd(mat_ruv,0,2)
plot_ly(x = mat.pca$v[,1], y = mat.pca$v[,2],type = "scatter",color = colData(dds)$Group) %>% layout(showlegend = F)

mat <- mat_ruv

genesToKeep <- rownames(mat)
for(b in levels(colData(dds)$Batch)){
  mat_b <- mat[,colData(dds)$Sample[colData(dds)$Batch == b]]
  vars <- apply(mat_b,1,var)
  genesToKeep <- intersect(genesToKeep,names(which(vars > 0)))
}
genesToKeep <- sort(unique(genesToKeep))
mat <- mat[genesToKeep,]
dim(mat)

inp_sva_train <- t(mat[,trainSamps])
inp_sva_test <- t(mat[,testSamps])
age_dat_scaled <- colData(dds)[rownames(inp_sva_train),"Age"]

inp_sva_train <- cbind(age_dat_scaled,inp_sva_train)
colnames(inp_sva_train)[1] <- "Age_scaled"
inp_sva_train[1:5,1:5]

#Read MBOTC and prepare processes
mbotc <- read.table("./Data/gene-SCPassociations.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")
all_procs <- setdiff(unique(mbotc$ProcessName),"Background genes")
genes_to_procs <- lapply(all_procs,function(x){unique(mbotc$Symbol[mbotc$ProcessName == x])})
names(genes_to_procs) <- all_procs
ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")
symbToEns <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = unique(mbotc$Symbol), mart= ensembl)
ensembl_to_procs <- lapply(genes_to_procs,function(x){intersect(colnames(inp_sva_train),unique(symbToEns$ensembl_gene_id[symbToEns$hgnc_symbol %in% x]))})
ensembl_to_procs <- ensembl_to_procs[unname(which(sapply(ensembl_to_procs,length) >= 3))]
level_to_procs <- sapply(names(ensembl_to_procs),function(x){unique(mbotc$ProcessLevel[mbotc$ProcessName == x])})
levels_to_consider <- c(3)
level_to_procs <- level_to_procs[level_to_procs %in% levels_to_consider] 

processes_to_consider <- c("WNT-Beta-catenin signaling pathway",
                           "Histone methylation and demethylation",
                           "Adherens junction organization",
                           "Translation initiation",
                           "Actin polymerization",
                           "Triacylglycerol transport by lipoproteins",
                           "Endoplasmic reticulum quality control system",
                           "Sodium transmembrane transport")

obs_weights <- rep(1,nrow(inp_sva_train))
genesInLevel <- unique(do.call("c",ensembl_to_procs[processes_to_consider]))

selectGenes <- list()
selectGenes_fraction <- rep(-1,length(processes_to_consider))
names(selectGenes_fraction) <-processes_to_consider
adjr2 <- rep(-1,length(processes_to_consider))
names(adjr2) <- processes_to_consider
for(proc in processes_to_consider){
  print(proc)
  inp_comp_level <- inp_sva_train
  g <- intersect(colnames(inp_comp_level),ensembl_to_procs[[proc]])
  inp_comp_level <- cbind(inp_sva_train[,c(g,"Age_scaled")])
  inp_comp_level <- as.data.frame(inp_comp_level)
  inp_comp_level[,-ncol(inp_comp_level)] <- apply(inp_comp_level[,-ncol(inp_comp_level)],2,scale)
  
  vargenes <- names(which(apply(inp_comp_level,2,function(x){!any(is.nan(x))})))
  inp_comp_level <- inp_comp_level[,vargenes]
  
  # Step 1: Define base intercept only model
  base.mod <- lm(Age_scaled ~ 1 , data=inp_comp_level,weights = obs_weights)  
  
  # Step 2: Full model with all predictors
  all.mod <- lm(Age_scaled ~ ., data= inp_comp_level,weights = obs_weights) 
  
  # Step 3: Perform step-wise algorithm. direction='both' implies both forward and backward stepwise
  stepMod <- step(base.mod, scope = list(lower = base.mod, upper = all.mod), direction = "both", trace = 0, steps = 5000, k = 2)  
  
  # Step 4: Get the shortlisted variable.
  shortlistedVars <- names(unlist(stepMod[[1]])) 
  selectGenes[[proc]] <- shortlistedVars[!shortlistedVars %in% "(Intercept)"] # remove intercept
  sumStepMod <- summary(stepMod)
  adjr2[proc] <- sumStepMod[["adj.r.squared"]]
  selectGenes_fraction[proc] <- length(shortlistedVars[!shortlistedVars %in% "(Intercept)"])/length(g)
}

x <- sapply(selectGenes,length)
y <- adjr2
plot_ly(x = x, y = y, color = names(y)) %>% layout(showlegend = F) #%>% add_lines(x = s$x, y = s$y)

genesInLevel <- unique(do.call("c",selectGenes[processes_to_consider]))

inp_comp_level <- as.data.frame(cbind(inp_sva_train[,c(genesInLevel,"Age_scaled")],obs_weights))

require(h2o)
conn <- h2o.init(max_mem_size="8G")
inp_comp_level_h2o <- as.h2o(inp_comp_level)

model <- h2o.glm(y = "Age_scaled",
                 training_frame = inp_comp_level_h2o,
                 nfolds = 10,
                 fold_assignment = "Random",
                 family = "AUTO",
                 link = "identity",
                 lambda_search = T,
                 standardize = T,
                 weights_column = "obs_weights",
                 alpha = 0,
                 seed = 100,
                 max_active_predictors = ncol(inp_comp_level_h2o),
                 solver = "IRLSM"
)
model

inp_comp_asl <- as.data.frame(cbind(inp_sva_test[,genesInLevel]))
inp_comp_asl_h2o <- as.h2o(inp_comp_asl)

preds <- as.data.frame(h2o.predict(model,inp_comp_asl_h2o))
aslpreds <- sapply(preds$predict, h_inv)
aslpreds<- data.frame(cbind(SampleID=rownames(inp_comp_asl),Age=as.numeric(aslpreds)), stringsAsFactors=FALSE)
aslpreds$Age <- as.numeric(aslpreds$Age)
aslpreds <- merge(meta_comb_new,aslpreds,by="SampleID")
aslpreds <- aslpreds[order(aslpreds$Gene),]
aslpreds$Agediff <- aslpreds$Age - aslpreds$Line_age

preds_inbag <- as.data.frame(h2o.predict(model,inp_comp_level_h2o))
preds_inbag <- sapply(preds_inbag$predict, h_inv)
preds_inbag<- data.frame(cbind(SampleID=rownames(inp_comp_level),Age=as.numeric(preds_inbag)), stringsAsFactors=FALSE)
preds_inbag$Age <- as.numeric(preds_inbag$Age)
preds_inbag <- merge(meta_comb_new,preds_inbag,by="SampleID")
preds_inbag <- preds_inbag[order(preds_inbag$Gene),]
preds_inbag$Agediff <- preds_inbag$Age - preds_inbag$Line_age

preds_total <- rbind(preds_inbag,aslpreds)
View(preds_total)

#coeff <- model@model$coefficients_table

#write.table(coeff,file = "./LIBSRR122069_Coefficients_final.txt", sep = "\t", quote = F, row.names = F, col.names = T)
#write.table(preds_total,file = "./LIBSRR122069_Preds_final.txt", sep = "\t", quote = F, row.names = F, col.names = T)

coeff <- model@model$coefficients_table
coeff_vec <- coeff$coefficients
names(coeff_vec) <- coeff$names

mat_scaled <- mat[,setdiff(colnames(mat),fleischersamps)]

libToSymb <- c("Fibroblast"="Fibroblast")

meta <- meta_comb_new

mat_scaled <- mat_scaled[,meta$SampleID]

all(meta$SampleID == colnames(mat_scaled))

meta$LineID <- libToSymb[meta$Line]
meta$SampleID <- paste0(meta$LineID,gsub("S","",meta$SampleID))

colnames(mat_scaled) <- meta$SampleID

genes <- setdiff(unique(meta$Gene),c("NTg","mTagBFP2"))
idx_ntg <- meta$SampleID[meta$Gene == "NTg"]
idx_bfp <- meta$SampleID[meta$Gene == "mTagBFP2"]
process_activity <- list()
for(p in processes_to_consider){
  print(p)
  genes_in_proc <- ensembl_to_procs[[p]]
  coeffs_model <- coeff_vec[genes_in_proc]
  coeffs_model <- coeffs_model[!is.na(coeffs_model)]
  
  sample_list <- list()
  for(g in genes){
    idx_g <- meta$SampleID[meta$Gene == g]
    mat_scaled_proc <- mat_scaled[names(coeffs_model),c(idx_g,idx_ntg,idx_bfp)]
    
    tmp <- as.data.frame((t(mat_scaled_proc) %*% coeffs_model))
    tmp$Gene <- meta$Gene[match(rownames(tmp),meta$SampleID)]
    tmp$Line <- gsub("[0-9]+","",rownames(tmp))
    colnames(tmp)[1] <- "Val"
    sample_list[[g]] <- tmp
  }
  
  process_activity[[p]] <- Reduce(function(x, y) join(x,y,by = c("Line","Gene"), type = "full", match = "first"),sample_list)
  #process_activity[[p]] <- ddply(process_activity[[p]],c("Line","Gene"),transform,Replicate=seq(from=1,by=1,length.out=length(Line)))
  #process_activity[[p]]$LineRep <- paste0(process_activity[[p]]$Line,"_",process_activity[[p]]$Replicate)
  #process_activity[[p]]$Replicate <- NULL
  #process_activity[[p]]$Line <- NULL
  process_activity[[p]] <- aggregate(process_activity[[p]]$Val,by = list(process_activity[[p]]$Line,process_activity[[p]]$Gene),FUN = mean)
  colnames(process_activity[[p]]) <- c("Line","Gene","Val")
  process_activity[[p]] <- reshape(process_activity[[p]], idvar = "Gene", timevar = "Line", direction = "wide")
  colnames(process_activity[[p]]) <- gsub("Val\\.","",colnames(process_activity[[p]]))
  #process_activity[[p]] <- process_activity[[p]][,c("Gene","A","A","B","B","C","C")]
  rownames(process_activity[[p]]) <- process_activity[[p]]$Gene
  process_activity[[p]]$Gene <- NULL
  #for(i in 1:ncol(process_activity[[p]])){
  #  process_activity[[p]][is.na(process_activity[[p]][,i]),i] <- process_activity[[p]]["NTg",i]
  #}
  #process_activity[[p]]$A[is.na(process_activity[[p]]$A)] <- process_activity[[p]]["NTg","A"]
  #process_activity[[p]]$B[is.na(process_activity[[p]]$B)] <- process_activity[[p]]["NTg","B"]
  #process_activity[[p]]$C[is.na(process_activity[[p]]$C)] <- process_activity[[p]]["NTg","C"]
}


mat_fleischer <- mat[,fleischersamps]

process_activity_fleischer <- list()
for(p in processes_to_consider){
  print(p)
  genes_in_proc <- ensembl_to_procs[[p]]
  coeffs_model <- coeff_vec[genes_in_proc]
  coeffs_model <- coeffs_model[!is.na(coeffs_model)]
  
  mat_fleischer_proc <- mat_fleischer[names(coeffs_model),]
  tmp <- as.data.frame((t(mat_fleischer_proc) %*% coeffs_model))
  colnames(tmp)[1] <- p
  
  process_activity_fleischer[[p]] <- tmp
  
}

process_activity_fleischer <- do.call("cbind",process_activity_fleischer)

save(libToSymb,process_activity,process_activity_fleischer,file = "~/Documents/ProcessActivity_LIBSRR122069.RData")





