library(sva)
library(h2o)

horvath_inverse_transform <- function(transformed_age) {
  ifelse(transformed_age <= 0,
         exp(transformed_age + log(21)) - 1,
         21 * transformed_age + 20)
}

computeTPM <- function(counts, lengths) {
  
  rpk <- counts / (1.0*lengths / 1e3)         # reads per kilobase
  tpm <- rpk / sum(1.0*rpk) * 1e6             # scale to 1 million
  return(tpm)
}

# INPUT DATAFRAME - RAW COUNTS, ENSEMBL IDs AS ROWNAMES
samples_new_apply <- readRDS("./AllSamples_raw.rds")
rownames(samples_new_apply) <- samples_new_apply$Geneid
samples_new_apply$Geneid <- NULL
# INPUT METADATA - DATAFRAME, REQUIRES ONE COLUMN called 'LIBPrep'. Factor/character with two values: 'Total' and 'mRNA'
meta_comb_apply <- readRDS("./AllSamples_meta.rds")

orig.sva <- new.env()
load("./SVAModel.RData",envir = orig.sva)
gl <- readRDS("./GeneLengths.rds")
gTPM <- readRDS("./GenesForTPMNorm.rds")

samples_new_apply <- samples_new_apply[gTPM,]

samples_new_apply <- apply(samples_new_apply, 2, computeTPM, lengths = unname(gl[rownames(samples_new_apply)]))
samples_new_apply <- log(samples_new_apply+1)

samples_new_apply <- samples_new_apply[rownames(orig.sva$samples_new_train_train),]
samples_new_apply[1:5,1:5]

samples_new_apply_corrected <- fsva(as.matrix(orig.sva$samples_new_train_train),
                                    model.matrix(~Age + Gene + Dox,orig.sva$meta_comb_train_train),
                                    orig.sva$out,
                                    as.matrix(samples_new_apply),
                                    "exact")

samples_new_apply_corrected <- samples_new_apply_corrected$new
samples_new_apply_corrected[1:5,1:5]

libtype <- meta_comb_apply$LIBPrep
names(libtype) <- meta_comb_apply$SampleID
libtype <- as.factor(libtype)
samples_new_apply_corrected <- as.data.frame(t(samples_new_apply_corrected))
inp_sva_test <- cbind(unname(libtype[rownames(samples_new_apply_corrected)]),samples_new_apply_corrected)
colnames(inp_sva_test)[1] <- "LibType"
inp_sva_test[1:5,1:5]

conn <- h2o.init(max_mem_size="12G",nthreads = 4)
samples_new_apply_corrected_h2o <- as.h2o(inp_sva_test)
# THIS MUST MATCH THE MODEL THAT IS SAVED BY 'TrainClock.R'
# H2O GIVES EACH SAVED MODEL A NEW NAME --> THIS NEEDS TO BE ADAPTED!
model <- h2o.loadModel("./GLM_model_R_1792823899860_1")

preds_apply <- h2o.predict(model,newdata = samples_new_apply_corrected_h2o)
preds_apply <- as.data.frame(preds_apply)
preds_apply$realpredage <- horvath_inverse_transform(preds_apply$predict)

meta_comb_apply_wPreds <- meta_comb_apply
meta_comb_apply_wPreds$Preds <- preds_apply$realpredage
View(meta_comb_apply_wPreds)

saveRDS(meta_comb_apply_wPreds,file = "./Applied.rds")





