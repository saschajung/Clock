library(sva)
library(h2o)

genesToFilter_df <- read.table("./genes_to_filter.txt",header = T, sep = "\t")
genesToKeep <- setdiff(rownames(samples_new_apply),genesToFilter_df$ensembl_id)
samples_new_apply <- samples_new_apply[genesToKeep,]

samples_new_apply <- apply(samples_new_apply, 2, computeTPM, lengths = unname(gl[rownames(samples_new_apply)]))
samples_new_apply <- log(samples_new_apply+1)

samples_new_apply <- samples_new_apply[rownames(samples_new_train_corrected),]
samples_new_apply[1:5,1:5]

samples_new_apply_corrected <- fsva(as.matrix(samples_new_train_train),
                                    model.matrix(~Age + Gene + Dox,meta_comb_train_train),
                                    out,
                                    as.matrix(samples_new_apply),
                                    "exact")

samples_new_apply_corrected <- samples_new_apply_corrected$new
samples_new_apply_corrected[1:5,1:5]

samples_new_apply_corrected <- samples_new_apply

libtype <- meta_comb_apply$LIBPrep
names(libtype) <- meta_comb_apply$SampleID
libtype <- as.factor(libtype)
samples_new_apply_corrected <- as.data.frame(t(samples_new_apply_corrected))
inp_sva_test <- cbind(unname(libtype[rownames(samples_new_apply_corrected)]),samples_new_apply_corrected)
colnames(inp_sva_test)[1] <- "LibType"
inp_sva_test[1:5,1:5]

samples_new_apply_corrected_h2o <- as.h2o(inp_sva_test)

preds_apply <- h2o.predict(model,newdata = samples_new_apply_corrected_h2o)
preds_apply <- as.data.frame(preds_apply)
preds_apply$realpredage <- horvath_inverse_transform(preds_apply$predict) 

meta_comb_apply_wPreds <- meta_comb_apply
meta_comb_apply_wPreds$Preds <- preds_apply$realpredage
View(meta_comb_apply_wPreds)

saveRDS(meta_comb_apply_wPreds,file = "./Applied.rds")





