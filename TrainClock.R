library(DescTools)
library(sva)
library(psych)
library(h2o)
library(plotly)
library(biomaRt)

setwd("./")

samples_new_train <- readRDS("../TrainingData_raw.rds")
samples_new_df <- readRDS("../AllSamples_raw.rds")
meta_comb_new <- readRDS("../AllSamples_meta.rds")
meta_comb_train <- readRDS("../TrainingData_meta.rds")
rownames(samples_new_train) <- samples_new_df$Geneid

horvath_transform <- function(age) {
  ifelse(age <= 20,
         log(age + 1) - log(21),
         (age - 20) / 21)
}

horvath_inverse_transform <- function(transformed_age) {
  ifelse(transformed_age <= 0,
         exp(transformed_age + log(21)) - 1,
         21 * transformed_age + 20)
}

computeTPM <- function(counts, lengths) {
  
  rpk <- counts / (lengths / 1e3)         # reads per kilobase
  tpm <- rpk / sum(rpk) * 1e6             # scale to 1 million
  return(tpm)
}

#Remove samples with less than 10 million counts

#Compute weights for ages
create_age_weights <- function(
    age_vector,
    nbins = 20,
    method = c("inverse_freq", "target_uniform"),
    binning = c("quantile", "equal_width")
) {
  method  <- match.arg(method)
  binning <- match.arg(binning)
  
  # --- Bin ages ---
  if (binning == "quantile") {
    bins <- cut(age_vector,
                breaks = unique(quantile(age_vector, probs = seq(0, 1, length.out = nbins + 1))),
                include.lowest = TRUE)
  } else {
    bins <- cut(age_vector, breaks = nbins, include.lowest = TRUE)
  }
  
  # counts per bin
  bin_counts <- table(bins)
  
  # --- weight calculation ---
  if (method == "inverse_freq") {
    # weight = 1 / (count in bin)
    w <- 1 / bin_counts[bins]
  } else if (method == "target_uniform") {
    # weight = target_prob / actual_prob
    target_prob <- 1 / length(bin_counts)
    actual_prob <- bin_counts / sum(bin_counts)
    w <- target_prob / actual_prob[bins]
  }
  
  # normalize weights to mean = 1 (h2o-friendly)
  w <- as.numeric(w)
  w <- w / mean(w)
  
  return(w)
}

weights_new <- create_age_weights(meta_comb_train$Line_age,nbins = 10,method = "target_uniform",binning = "equal_width")
weights_new <- round(weights_new*10)
plot_ly(x = meta_comb_train$Line_age,y = weights_new)

genesToFilter_df <- read.table("./genes_to_filter.txt",header = T, sep = "\t")
genesToKeep <- setdiff(rownames(samples_new_train),genesToFilter_df$ensembl_id)
samples_new_train <- samples_new_train[genesToKeep,]

#TPM and Log-normalize the training data
samples_new_train <- apply(samples_new_train, 2, computeTPM, lengths = unname(gl[rownames(samples_new_train)]))
samples_new_train <- log(samples_new_train+1)
samples_new_train[1:5,1:5]

saveRDS(rownames(samples_new_train),"../GenesForTPMNorm.rds")

rmeans <- rowMeans(samples_new_train)
samples_new_train <- samples_new_train[which(rmeans >= 0.5),]

#Perform surrogate variable analysis on a subset of samples
meta_comb_train$Age <- horvath_transform(meta_comb_train$Line_age)
meta_comb_train$Dox <- as.numeric(meta_comb_train$Dox)

set.seed(100)
trainSamples <- sample(meta_comb_train$SampleID,0.66*nrow(meta_comb_train),replace = F)
testSamples <- setdiff(meta_comb_train$SampleID,trainSamples)

meta_comb_train_train <- meta_comb_train[which(meta_comb_train$SampleID %in% trainSamples),]
meta_comb_train_test <- meta_comb_train[which(meta_comb_train$SampleID %in% testSamples),]

samples_new_train_train <- samples_new_train[,meta_comb_train_train$SampleID]

mod <- model.matrix(~Age + Gene + Dox,meta_comb_train_train) #Keep age as the covariate of interest
mod0 <- model.matrix(~Gene + Dox,meta_comb_train_train)
nsv_be <- num.sv(samples_new_train_train, mod, vfilter = NULL, method = "be", seed = 100)
nsv_be

out <- sva(dat = as.matrix(samples_new_train_train),
           mod = mod,
           mod0 = mod0,
           vfilter = NULL,
           n.sv = nsv_be)

save(meta_comb_train_train,samples_new_train_train,out,mod,mod0,nsv_be,file = "./SVAModel.RData")

samples_new_train_corrected <- fsva(as.matrix(samples_new_train_train),mod,out,as.matrix(samples_new_train),"exact")
samples_new_train_corrected <- samples_new_train_corrected$new

libtype <- meta_comb_train$LIBPrep
names(libtype) <- meta_comb_train$SampleID
libtype <- as.factor(libtype)
inp_sva_train <- cbind(meta_comb_train$Age,t(samples_new_train_corrected[,]))
colnames(inp_sva_train)[1] <- "Age"

inp_sva_train <- cbind(rep(1,length(weights_new)),inp_sva_train)
colnames(inp_sva_train)[1] <- "Weight"

inp_sva_train <- as.data.frame(inp_sva_train)
inp_sva_train <- cbind(unname(libtype[rownames(inp_sva_train)]),inp_sva_train)
colnames(inp_sva_train)[1] <- "LibType"
inp_sva_train[1:5,1:5]

#Read MBOTC and prepare processes
mbotc <- read.table("./Data/gene-SCPassociations.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")
all_procs <- setdiff(unique(mbotc$ProcessName),"Background genes")
genes_to_procs <- lapply(all_procs,function(x){unique(mbotc$Symbol[mbotc$ProcessName == x])})
names(genes_to_procs) <- all_procs
ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")
symbToEns <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = unique(mbotc$Symbol), mart= ensembl)
ensembl_to_procs <- lapply(genes_to_procs,function(x){unique(symbToEns$ensembl_gene_id[symbToEns$hgnc_symbol %in% x])})
ensembl_to_procs <- lapply(ensembl_to_procs,function(x){intersect(colnames(inp_sva_train),x)})
ensembl_to_procs <- ensembl_to_procs[unname(which(sapply(ensembl_to_procs,length) >= 3))]

hpa_cts <- read.table("../rna_single_cell_type.tsv", sep = "\t", header = T)
hpa_cts$Gene.name <- NULL
hpa_cts_mat <- reshape(hpa_cts,timevar = "Cell.type", idvar = "Gene",direction = "wide")
rownames(hpa_cts_mat) <- hpa_cts_mat$Gene
hpa_cts_mat$Gene <- NULL

spec <- sapply(ensembl_to_procs,function(x){
  m <- colMeans(hpa_cts_mat[x,],na.rm = T)
  m <- m/max(m)
  return(m["nTPM.Fibroblasts"])
})
names(spec) <- gsub("\\.nTPM\\.Fibroblasts","",names(spec))
priorProcs <- c("WNT-Beta-catenin signaling pathway","Actin polymerization","Adherens junction organization","Histone methylation and demethylation","Endoplasmic reticulum quality control system","Translation initiation","Sodium transmembrane transport")

procs_to_remove <- setdiff(sort(names(which(spec <= mean(spec)))),priorProcs)

ensembl_to_procs <- ensembl_to_procs[setdiff(names(ensembl_to_procs),procs_to_remove)]

require(h2o)
conn <- h2o.init(max_mem_size="12G",nthreads = 4)

genesToConsider <- unique(do.call("c",ensembl_to_procs[priorProcs]))

model <- h2o.glm(y = "Age",
                 x = c(genesToConsider,"LibType"),
                 training_frame = inp_sva_train_h2o,
                 validation_frame = inp_sva_test_h2o,
                 nfolds = 5,
                 fold_assignment = "Random",
                 family = "AUTO",
                 link = "identity",
                 lambda_search = T,
                 standardize = T,
                 weights_column = "Weight",
                 keep_cross_validation_predictions = F,
                 alpha = 0,
                 seed = 100,
                 solver = "IRLSM"
)

h2o.saveModel(model, path = "./")

background_df <- data.frame(Genes = colnames(inp_sva_train_h2o)[-c(1:3)],stringsAsFactors = F)
write.table(background_df,file = "./BackgroundGenes.txt",sep = "\t",quote = F, col.names = T, row.names = F)
save(trainSamples,testSamples,file = "./TrainTestSamples.RData")

