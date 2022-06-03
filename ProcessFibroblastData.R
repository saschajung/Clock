library(biomaRt)
library(plyr)
library(stringr)
#### Read and process training data ####
#Read data
train_data <- read.table("./Data/Fleischer_featurecounts.Rmatrix.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)

train_data <- train_data[,c(1,order(colnames(train_data)[-1])+1)]
train_data <- train_data[order(train_data$Geneid),]

#Cleanup GeneId
train_data$Geneid <- gsub("\\.[0-9]+","",train_data$Geneid)
ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = train_data$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
train_data <- train_data[train_data$Geneid %in% ensTolen_agg$Group.1,]
all(train_data$Geneid == ensTolen_agg$Group.1)
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]

#Add and process metadata
train_metadata <- read.table("./Data/Fleischer_metadata.txt", sep = ',', header = TRUE, stringsAsFactors = FALSE)
train_metadata <- train_metadata[,c("Run","AGE","sex","ETHNICITY","disease")]
train_metadata <- train_metadata[order(train_metadata$Run),]
colnames(train_metadata) <- c("Sample","Age","Sex","Race","Disease")
train_metadata$Age_Y <- train_metadata$Age
train_metadata$Age_Y <- gsub("YR$","",train_metadata$Age_Y)
train_metadata$Age_Y <- gsub("yr$","",train_metadata$Age_Y)
idx_mos <- which(grepl("mos",train_metadata$Age_Y))
train_metadata$Age_Y[idx_mos] <- as.character(unname(sapply(train_metadata$Age_Y[idx_mos],function(x){
  x <- gsub("mos","",x)
  x <- strsplit(x,"yr|ys")[[1]]
  return(round(as.numeric(x[1]) + as.numeric(x[2])/12))
})))
train_metadata$Age_Y <- as.numeric(train_metadata$Age_Y)
train_metadata$Age_M <- train_metadata$Age
train_metadata$Age_M <- gsub("YR$","",train_metadata$Age_M)
train_metadata$Age_M <- gsub("yr$","",train_metadata$Age_M)
idx_mos <- which(grepl("mos",train_metadata$Age_M))
train_metadata$Age_M[idx_mos] <- as.character(unname(sapply(train_metadata$Age_M[idx_mos],function(x){
  x <- gsub("mos","",x)
  x <- strsplit(x,"yr|ys")[[1]]
  return(as.numeric(x[1]) + as.numeric(x[2])/12)
})))
train_metadata$Age_M <- as.numeric(train_metadata$Age_M)
train_metadata$Age_M <- train_metadata$Age_M*12
train_metadata$Sex[which(train_metadata$Sex == "Female" | train_metadata$Sex == "female")] <- "F"
train_metadata$Sex[which(train_metadata$Sex == "Male" | train_metadata$Sex == "male")] <- "M"
train_metadata$Race[which(train_metadata$Race == "caucasian")] <- "Caucasian"
train_metadata$Race[which(train_metadata$Race == "CaucasianSardinian")] <- "Caucasian"

hgps_meta <- data.frame(SampleID = train_metadata$Sample[train_metadata$Disease == "HGPS"],
                             GeneID = 0,
                             Replicate = 1,
                             Gene = "HGPS",
                             gene_symbol = "HGPS",
                             Line = "Fibroblast",
                             Passage = -1,
                             Line_age = train_metadata$Age_Y[train_metadata$Disease == "HGPS"],
                             Dox = 0,
                             LIB = "train",
                             Comment = ""
)

hgps <- train_data[,c("Geneid",hgps_meta$SampleID)]

#### Read ASL data ####
#Read data
ASL_20M <- read.table("./Data/all_data_featurecounts(1).txt",header = T)#read.table("./Data/20M_full_featurecounts.txt",header = TRUE)
ASL_meta <- read.csv("./Data/all_metadata.csv",header = TRUE)#read.csv("./Data/ASL_metadata_full.csv",header = TRUE)
ASL_meta <- ASL_meta[!is.na(ASL_meta$Line_age),]
ASL_meta$Gene[grepl("iPS",ASL_meta$Comment) & (ASL_meta$Gene == "NTg")] <- "NTg_iPS"
ASL_20M <- ASL_20M[,c("Geneid",ASL_meta$SampleID)]

#Cleanup GeneId
ASL_20M$Geneid <- gsub("\\.[0-9]+","",ASL_20M$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = ASL_20M$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
ASL_20M <- ASL_20M[ASL_20M$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
ASL_20M <- ASL_20M[order(ASL_20M$Geneid),]
all(ASL_20M$Geneid == ensTolen_agg$Group.1)

#Read additional fibro data
#SRR303074 (MRC5 cell line/ 14 weeks foetus, mRNA)
SRR303074 <- read.table("./Data/1_SRR303074_featurecounts.txt", header = T, sep = "\t")
gene_length <- SRR303074$Length
names(gene_length) <- gsub("\\.[0-9]+","",SRR303074$Geneid)
SRR303074_meta <- read.table("./Data/1_SraRunTable.txt", header = T, sep = ",")

SRR303074[ ,2:6] <- list(NULL)
colnames(SRR303074)[-1] <- str_extract(colnames(SRR303074)[-1],"SRR[0-9]+")
SRR303074$Geneid <- gsub("\\.[0-9]+","",SRR303074$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR303074$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR303074 <- SRR303074[SRR303074$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR303074 <- SRR303074[order(SRR303074$Geneid),]
all(SRR303074$Geneid == ensTolen_agg$Group.1)

SRR303074_meta <- SRR303074_meta[,c("Run","Cell_Line","cell_state","source_name")]
SRR303074_meta <- ddply(SRR303074_meta,c("cell_state"),transform,Replicate=seq(from=1,by=1,length.out=length(cell_state)))
SRR303074_meta$Gene <- paste0("NTg_",SRR303074_meta$source_name)

SRR303074_meta <- data.frame(SampleID = SRR303074_meta$Run,
                             GeneID = 0,
                             Replicate = SRR303074_meta$Replicate,
                             Gene = SRR303074_meta$Gene,
                             gene_symbol = SRR303074_meta$Gene,
                             Line = SRR303074_meta$Cell_Line,
                             Passage = -1,
                             Line_age = 0,
                             Dox = 0,
                             LIB = "LIBSRR303074",
                             Comment = SRR303074_meta$cell_state
)

commonsamps <- intersect(SRR303074_meta$Sample,colnames(SRR303074))
SRR303074_meta <- SRR303074_meta[SRR303074_meta$Sample %in% commonsamps,]
SRR303074 <- SRR303074[,c("Geneid",commonsamps)]

all(SRR303074_meta$Sample == colnames(SRR303074)[-1])

#SRR121037 (BJ cell line/ newborn/ mRNA)
SRR121037 <- read.table("./Data/2_SRR121037_featurecounts.txt", header = T, sep = "\t")
SRR121037_meta <- read.table("./Data/2_SraRunTable.txt", header = T, sep = ",")

SRR121037[ ,2:6] <- list(NULL)
colnames(SRR121037)[-1] <- str_extract(colnames(SRR121037)[-1],"SRR[0-9]+")
SRR121037$Geneid <- gsub("\\.[0-9]+","",SRR121037$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR121037$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR121037 <- SRR121037[SRR121037$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR121037 <- SRR121037[order(SRR121037$Geneid),]
all(SRR121037$Geneid == ensTolen_agg$Group.1)

SRR121037_meta <- SRR121037_meta[,c("Run","Cell_type","source_name","develpmental_stage","treatment")]
SRR121037_meta <- ddply(SRR121037_meta,c("develpmental_stage","treatment"),transform,Replicate=seq(from=1,by=1,length.out=length(develpmental_stage)))
SRR121037_meta$Treated <- "Treated"
SRR121037_meta$Treated[SRR121037_meta$treatment == "none"] <- "Untreated"
SRR121037_meta$Gene <- paste0("NTg_",SRR121037_meta$develpmental_stage,"_",SRR121037_meta$Treated)

SRR121037_meta <- data.frame(SampleID = SRR121037_meta$Run,
                             GeneID = 0,
                             Replicate = SRR121037_meta$Replicate,
                             Gene = SRR121037_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR121037_meta$Cell_type,
                             Passage = -1,
                             Line_age = 0,
                             Dox = 0,
                             LIB = "LIBSRR121037",
                             Comment = SRR121037_meta$treatment
)

commonsamps <- intersect(SRR121037_meta$Sample,colnames(SRR121037))
SRR121037_meta <- SRR121037_meta[SRR121037_meta$Sample %in% commonsamps,]
SRR121037 <- SRR121037[,c("Geneid",commonsamps)]

all(SRR121037_meta$Sample == colnames(SRR121037)[-1])

#SRR315202 (MRC5/HFFcell line/ foetus/ total RNA)
SRR315202 <- read.table("./Data/4_SRR315202_featurecounts.txt", header = T, sep = "\t")
SRR315202_meta <- read.table("./Data/4_SraRunTable.txt", header = T, sep = ",")

SRR315202[ ,2:6] <- list(NULL)
colnames(SRR315202)[-1] <- str_extract(colnames(SRR315202)[-1],"SRR[0-9]+")
SRR315202$Geneid <- gsub("\\.[0-9]+","",SRR315202$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR315202$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR315202 <- SRR315202[SRR315202$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR315202 <- SRR315202[order(SRR315202$Geneid),]
all(SRR315202$Geneid == ensTolen_agg$Group.1)

SRR315202_meta <- SRR315202_meta[,c("Run","agent","Cell_Line","population_doublings")]
SRR315202_meta <- ddply(SRR315202_meta,c("Cell_Line","population_doublings"),transform,Replicate=seq(from=1,by=1,length.out=length(population_doublings)))
SRR315202_meta$Cell_Line[SRR315202_meta$Cell_Line != "HFF"] <- "MRC5"
SRR315202_meta$agent <- gsub(" ","",SRR315202_meta$agent)
SRR315202_meta$Gene <- paste0("NTg_",SRR315202_meta$agent)

SRR315202_meta <- data.frame(SampleID = SRR315202_meta$Run,
                             GeneID = 0,
                             Replicate = SRR315202_meta$Replicate,
                             Gene = SRR315202_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR315202_meta$Cell_Line,
                             Passage = SRR315202_meta$population_doublings,
                             Line_age = 0,
                             Dox = 0,
                             LIB = "LIBSRR315202",
                             Comment = SRR315202_meta$agent
)

commonsamps <- intersect(SRR315202_meta$Sample,colnames(SRR315202))
SRR315202_meta <- SRR315202_meta[SRR315202_meta$Sample %in% commonsamps,]
SRR315202 <- SRR315202[,c("Geneid",commonsamps)]

all(SRR315202_meta$Sample == colnames(SRR315202)[-1])

#SRR275112 (MRC5/HFF/BJ/IMR90/Wi38 cell line/ foetus/ total RNA)
SRR275112 <- read.table("./Data/5_SRR275112_featurecounts.txt", header = T, sep = "\t")
SRR275112_meta <- read.table("./Data/5_SraRunTable.txt", header = T, sep = ",")

SRR275112[ ,2:6] <- list(NULL)
colnames(SRR275112)[-1] <- str_extract(colnames(SRR275112)[-1],"SRR[0-9]+")
SRR275112$Geneid <- gsub("\\.[0-9]+","",SRR275112$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR275112$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR275112 <- SRR275112[SRR275112$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR275112 <- SRR275112[order(SRR275112$Geneid),]
all(SRR275112$Geneid == ensTolen_agg$Group.1)

SRR275112_meta <- SRR275112_meta[,c("Run","Cell_type","age_category","population_doublings")]
SRR275112_meta$age_category[SRR275112_meta$age_category == "-"] <- ""
SRR275112_meta$Cell_type[SRR275112_meta$Cell_type == "BJ fibroblasts"] <- "BJ"
SRR275112_meta$Cell_type[SRR275112_meta$Cell_type == "HFF fibroblasts"] <- "HFF"
SRR275112_meta$Cell_type[SRR275112_meta$Cell_type == "IMR-90 fibroblasts"] <- "IMR90"
SRR275112_meta$Cell_type[SRR275112_meta$Cell_type == "MRC-5 fibroblasts"] <- "MRC5"
SRR275112_meta$Cell_type[SRR275112_meta$Cell_type == "MRC-5"] <- "MRC5"
SRR275112_meta$Cell_type[SRR275112_meta$Cell_type == "Wi-38 fibroblasts"] <- "WI38"
SRR275112_meta <- ddply(SRR275112_meta,c("Cell_type","age_category"),transform,Replicate=seq(from=1,by=1,length.out=length(age_category)))

SRR275112_meta <- data.frame(SampleID = SRR275112_meta$Run,
                             GeneID = 0,
                             Replicate = SRR275112_meta$Replicate,
                             Gene = "NTg",
                             gene_symbol = "NTg",
                             Line = SRR275112_meta$Cell_type,
                             Passage = SRR275112_meta$population_doublings,
                             Line_age = 0,
                             Dox = 0,
                             LIB = "LIBSRR275112",
                             Comment = SRR275112_meta$age_category
)

commonsamps <- intersect(SRR275112_meta$Sample,colnames(SRR275112))
SRR275112_meta <- SRR275112_meta[SRR275112_meta$Sample %in% commonsamps,]
SRR275112 <- SRR275112[,c("Geneid",commonsamps)]

all(SRR275112_meta$Sample == colnames(SRR275112)[-1])

#ERR137716 (dermal fibroblasts/ mRNA)
ERR137716 <- read.table("./Data/6_ERR137716_featurecounts.txt", header = T, sep = "\t")
ERR137716_meta <- read.table("./Data/6_SraRunTable.txt", header = T, sep = ",")

ERR137716[ ,2:6] <- list(NULL)
colnames(ERR137716)[-1] <- str_extract(colnames(ERR137716)[-1],"ERR[0-9]+")
ERR137716$Geneid <- gsub("\\.[0-9]+","",ERR137716$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = ERR137716$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
ERR137716 <- ERR137716[ERR137716$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
ERR137716 <- ERR137716[order(ERR137716$Geneid),]
all(ERR137716$Geneid == ensTolen_agg$Group.1)

ERR137716_meta <- ERR137716_meta[,c("Run","Age","Cell_type","sampling_site","sex","stimulus")]
ERR137716_meta$stimulus[ERR137716_meta$stimulus == "UV exposed"] <- "UV"
ERR137716_meta$Cell_type <- gsub(" ","",ERR137716_meta$Cell_type)

ERR137716_meta <- ddply(ERR137716_meta,c("Cell_type","sampling_site","sex","stimulus","Age"),transform,Replicate=seq(from=1,by=1,length.out=length(Age)))

ERR137716_meta$Gene <- "NTg"
ERR137716_meta$Gene[ERR137716_meta$stimulus == "UV"] <- "NTg_UV"

ERR137716_meta <- data.frame(SampleID = ERR137716_meta$Run,
                             GeneID = 0,
                             Replicate = ERR137716_meta$Replicate,
                             Gene = ERR137716_meta$Gene,
                             gene_symbol = "NTg",
                             Line = ERR137716_meta$Cell_type,
                             Passage = -1,
                             Line_age = ERR137716_meta$Age,
                             Dox = 0,
                             LIB = "LIBERR137716",
                             Comment = paste0(ERR137716_meta$sampling_site,";",ERR137716_meta$sex)
)

commonsamps <- intersect(ERR137716_meta$Sample,colnames(ERR137716))
ERR137716_meta <- ERR137716_meta[ERR137716_meta$Sample %in% commonsamps,]
ERR137716 <- ERR137716[,c("Geneid",commonsamps)]

all(ERR137716_meta$Sample == colnames(ERR137716)[-1])

#SRR101478 (skin fibroblasts/ total RNA)
SRR101478 <- read.table("./Data/8_SRR101478_featurecounts.txt", header = T, sep = "\t")
SRR101478_meta <- read.table("./Data/8_SraRunTable.txt", header = T, sep = ",")

SRR101478[ ,2:6] <- list(NULL)
colnames(SRR101478)[-1] <- str_extract(colnames(SRR101478)[-1],"SRR[0-9]+")
SRR101478$Geneid <- gsub("\\.[0-9]+","",SRR101478$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR101478$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR101478 <- SRR101478[SRR101478$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR101478 <- SRR101478[order(SRR101478$Geneid),]
all(SRR101478$Geneid == ensTolen_agg$Group.1)

SRR101478_meta <- SRR101478_meta[,c("Run","tissue","time_point")]
SRR101478_meta$time_point[SRR101478_meta$time_point == "early timepoint (E)"] <- "EarlyTime"
SRR101478_meta$time_point[SRR101478_meta$time_point == "late timepoint (L)"] <- "LateTime"
SRR101478_meta$tissue[SRR101478_meta$tissue == "primary skin fibroblast"] <- "Skin"

SRR101478_meta <- ddply(SRR101478_meta,c("tissue","time_point"),transform,Replicate=seq(from=1,by=1,length.out=length(time_point)))

SRR101478_meta$Gene <- "NTg"

SRR101478_meta <- data.frame(SampleID = SRR101478_meta$Run,
                             GeneID = 0,
                             Replicate = SRR101478_meta$Replicate,
                             Gene = SRR101478_meta$Gene,
                             gene_symbol = "NTg",
                             Line = paste0(SRR101478_meta$tissue,"_",SRR101478_meta$time_point),
                             Passage = -1,
                             Line_age = -1,
                             Dox = 0,
                             LIB = "LIBSRR101478",
                             Comment = SRR101478_meta$time_point
)

commonsamps <- intersect(SRR101478_meta$Sample,colnames(SRR101478))
SRR101478_meta <- SRR101478_meta[SRR101478_meta$Sample %in% commonsamps,]
SRR101478 <- SRR101478[,c("Geneid",commonsamps)]

all(SRR101478_meta$Sample == colnames(SRR101478)[-1])

#SRR173636 (HFF/MRC5/ newborn/ total RNA)
SRR173636 <- read.table("./Data/9_SRR173636_featurecounts.txt", header = T, sep = "\t")
SRR173636_meta <- read.table("./Data/9_SraRunTable.txt", header = T, sep = ",")

SRR173636[ ,2:6] <- list(NULL)
colnames(SRR173636)[-1] <- str_extract(colnames(SRR173636)[-1],"SRR[0-9]+")
SRR173636$Geneid <- gsub("\\.[0-9]+","",SRR173636$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR173636$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR173636 <- SRR173636[SRR173636$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR173636 <- SRR173636[order(SRR173636$Geneid),]
all(SRR173636$Geneid == ensTolen_agg$Group.1)

SRR173636_meta <- SRR173636_meta[,c("Run","Cell_Line","treatment","population_doublings")]
SRR173636_meta$treatment[SRR173636_meta$treatment == "100 nM rotenone"] <- "rotenone"

SRR173636_meta$Gene <- "NTg"
SRR173636_meta$Gene[SRR173636_meta$treatment != "none"] <- paste0(SRR173636_meta$Gene[SRR173636_meta$treatment != "none"],"_",SRR173636_meta$treatment[SRR173636_meta$treatment != "none"])

SRR173636_meta$Cell_Line[SRR173636_meta$Cell_Line == "MRC-5"] <- "MRC5"

SRR173636_meta <- ddply(SRR173636_meta,c("Cell_Line","population_doublings","treatment"),transform,Replicate=seq(from=1,by=1,length.out=length(treatment)))

SRR173636_meta <- data.frame(SampleID = SRR173636_meta$Run,
                             GeneID = 0,
                             Replicate = SRR173636_meta$Replicate,
                             Gene = SRR173636_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR173636_meta$Cell_Line,
                             Passage = SRR173636_meta$population_doublings,
                             Line_age = 0,
                             Dox = 0,
                             LIB = "LIBSRR173636",
                             Comment = SRR173636_meta$treatment
)

commonsamps <- intersect(SRR173636_meta$Sample,colnames(SRR173636))
SRR173636_meta <- SRR173636_meta[SRR173636_meta$Sample %in% commonsamps,]
SRR173636 <- SRR173636[,c("Geneid",commonsamps)]

all(SRR173636_meta$Sample == colnames(SRR173636)[-1])

#SRR176804 (foreskin/ mRNA)
SRR176804 <- read.table("./Data/10_SRR176804_featurecounts.txt", header = T, sep = "\t")
SRR176804_meta <- read.table("./Data/10_SraRunTable.txt", header = T, sep = ",")

SRR176804[ ,2:6] <- list(NULL)
colnames(SRR176804)[-1] <- str_extract(colnames(SRR176804)[-1],"SRR[0-9]+")
SRR176804$Geneid <- gsub("\\.[0-9]+","",SRR176804$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR176804$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR176804 <- SRR176804[SRR176804$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR176804 <- SRR176804[order(SRR176804$Geneid),]
all(SRR176804$Geneid == ensTolen_agg$Group.1)

SRR176804_meta <- SRR176804_meta[,c("Run","Cell_type","source")]
SRR176804_meta$Cell_type[SRR176804_meta$Cell_type == "proliferative fibroblast"] <- "prolif"
SRR176804_meta$Cell_type[SRR176804_meta$Cell_type == "queiscent fibroblast"] <- "quiescent"
SRR176804_meta$Cell_type[SRR176804_meta$Cell_type == "Rapamycin treated fibroblast"] <- "rapamycin"

SRR176804_meta$source <- "foreskin"

SRR176804_meta$Gene <- paste0("NTg_",SRR176804_meta$Cell_type)

SRR176804_meta <- ddply(SRR176804_meta,c("Cell_type"),transform,Replicate=seq(from=1,by=1,length.out=length(Cell_type)))

SRR176804_meta <- data.frame(SampleID = SRR176804_meta$Run,
                             GeneID = 0,
                             Replicate = SRR176804_meta$Replicate,
                             Gene = SRR176804_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR176804_meta$Cell_type,
                             Passage = -1,
                             Line_age = -1,
                             Dox = 0,
                             LIB = "LIBSRR176804",
                             Comment = SRR176804_meta$Cell_type
)

commonsamps <- intersect(SRR176804_meta$Sample,colnames(SRR176804))
SRR176804_meta <- SRR176804_meta[SRR176804_meta$Sample %in% commonsamps,]
SRR176804 <- SRR176804[,c("Geneid",commonsamps)]

all(SRR176804_meta$Sample == colnames(SRR176804)[-1])

#SRR627541 (BJ/ mRNA)
SRR627541 <- read.table("./Data/12_SRR627541_featurecounts.txt", header = T, sep = "\t")
SRR627541_meta <- read.table("./Data/12_SraRunTable.txt", header = T, sep = ",")

SRR627541[ ,2:6] <- list(NULL)
colnames(SRR627541)[-1] <- str_extract(colnames(SRR627541)[-1],"SRR[0-9]+")
SRR627541$Geneid <- gsub("\\.[0-9]+","",SRR627541$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR627541$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR627541 <- SRR627541[SRR627541$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR627541 <- SRR627541[order(SRR627541$Geneid),]
all(SRR627541$Geneid == ensTolen_agg$Group.1)

SRR627541_meta <- SRR627541_meta[SRR627541_meta$LibrarySelection == "cDNA",c("Run","Cell_Line","Condition")]
SRR627541_meta$Cell_Line <- "BJ"
SRR627541_meta$Condition[SRR627541_meta$Condition == "Proliferation\\, normal conditions"] <- "prolif"
SRR627541_meta$Condition[SRR627541_meta$Condition == "Quiescence induced by serum depletion"] <- "quiescence"
SRR627541_meta$Condition[SRR627541_meta$Condition == "pre-senescence; 5 dys after RASG12V induction"] <- "presenes"
SRR627541_meta$Condition[SRR627541_meta$Condition == "Senescence; 14 dys after RASG12V induction"] <- "senes"
SRR627541_meta$Condition[SRR627541_meta$Condition == "Transformed cells (induced by RASG12V in the background of stable p53 and p16INK4A kds)"] <- "transformed"

SRR627541_meta$Gene <- paste0("NTg_",SRR627541_meta$Condition)

SRR627541_meta <- ddply(SRR627541_meta,c("Condition"),transform,Replicate=seq(from=1,by=1,length.out=length(Condition)))

SRR627541_meta <- data.frame(SampleID = SRR627541_meta$Run,
                             GeneID = 0,
                             Replicate = SRR627541_meta$Replicate,
                             Gene = SRR627541_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR627541_meta$Cell_Line,
                             Passage = -1,
                             Line_age = -1,
                             Dox = 0,
                             LIB = "LIBSRR627541",
                             Comment = SRR627541_meta$Condition
)

commonsamps <- intersect(SRR627541_meta$Sample,colnames(SRR627541))
SRR627541_meta <- SRR627541_meta[SRR627541_meta$Sample %in% commonsamps,]
SRR627541 <- SRR627541[,c("Geneid",commonsamps)]

all(SRR627541_meta$Sample == colnames(SRR627541)[-1])

#SRR810099 (BJ/ mRNA)
SRR810099 <- read.table("./Data/13_SRR810099_featurecounts.txt", header = T, sep = "\t")
SRR810099_meta <- read.table("./Data/13_SraRunTable.txt", header = T, sep = ",")

SRR810099[ ,2:6] <- list(NULL)
colnames(SRR810099)[-1] <- str_extract(colnames(SRR810099)[-1],"SRR[0-9]+")
SRR810099$Geneid <- gsub("\\.[0-9]+","",SRR810099$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR810099$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR810099 <- SRR810099[SRR810099$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR810099 <- SRR810099[order(SRR810099$Geneid),]
all(SRR810099$Geneid == ensTolen_agg$Group.1)

SRR810099_meta <- SRR810099_meta[SRR810099_meta$LibrarySelection == "cDNA",c("Run","Cell_Line","Condition")]
SRR810099_meta$Cell_Line <- "BJ"
SRR810099_meta$Condition[SRR810099_meta$Condition == "Proliferation\\, normal conditions"] <- "prolif"
SRR810099_meta$Condition[SRR810099_meta$Condition == "Quiescence induced by serum depletion"] <- "quiescence"
SRR810099_meta$Condition[SRR810099_meta$Condition == "pre-senescence; 5 dys after RASG12V induction"] <- "presenes"
SRR810099_meta$Condition[SRR810099_meta$Condition == "Senescence; 14 dys after RASG12V induction"] <- "senes"
SRR810099_meta$Condition[SRR810099_meta$Condition == "Transformed cells (induced by RASG12V in the background of stable p53 and p16INK4A kds)"] <- "transformed"
SRR810099_meta$Condition[SRR810099_meta$Condition == "Control sample"] <- "none"
SRR810099_meta$Condition[SRR810099_meta$Condition == "Nutlin-3a\\, 2h"] <- "Nutlin3a_2h"
SRR810099_meta$Condition[SRR810099_meta$Condition == "Nutlin-3a\\, 4h"] <- "Nutlin3a_4h"
SRR810099_meta$Condition[SRR810099_meta$Condition == "Nutlin-3a\\, 6h"] <- "Nutlin3a_6h"
SRR810099_meta$Condition[SRR810099_meta$Condition == "Nutlin-3a\\, 19h"] <- "Nutlin3a_19h"

SRR810099_meta$Gene <- paste0("NTg_",SRR810099_meta$Condition)
SRR810099_meta$Gene[SRR810099_meta$Condition == "none"] <- "NTg"

SRR810099_meta <- ddply(SRR810099_meta,c("Condition"),transform,Replicate=seq(from=1,by=1,length.out=length(Condition)))

SRR810099_meta <- data.frame(SampleID = SRR810099_meta$Run,
                             GeneID = 0,
                             Replicate = SRR810099_meta$Replicate,
                             Gene = SRR810099_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR810099_meta$Cell_Line,
                             Passage = -1,
                             Line_age = -1,
                             Dox = 0,
                             LIB = "LIBSRR810099",
                             Comment = SRR810099_meta$Condition
)

commonsamps <- intersect(SRR810099_meta$Sample,colnames(SRR810099))
SRR810099_meta <- SRR810099_meta[SRR810099_meta$Sample %in% commonsamps,]
SRR810099 <- SRR810099[,c("Geneid",commonsamps)]

all(SRR810099_meta$Sample == colnames(SRR810099)[-1])

#SRR700080 (mixed/ total RNA)
SRR700080 <- read.table("./Data/14_SRR700080_featurecounts.txt", header = T, sep = "\t")
SRR700080_meta <- read.table("./Data/14_SraRunTable.txt", header = T, sep = ",")

SRR700080[ ,2:6] <- list(NULL)
colnames(SRR700080)[-1] <- str_extract(colnames(SRR700080)[-1],"SRR[0-9]+")
SRR700080$Geneid <- gsub("\\.[0-9]+","",SRR700080$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR700080$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR700080 <- SRR700080[SRR700080$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR700080 <- SRR700080[order(SRR700080$Geneid),]
all(SRR700080$Geneid == ensTolen_agg$Group.1)

SRR700080_meta <- SRR700080_meta[,c("Run","Passage","source_name","treatment")]
SRR700080_meta$Passage[SRR700080_meta$Passage == "P17-P20"] <- "17"
SRR700080_meta$Passage[SRR700080_meta$Passage == "P2-P5 (after reprogramming)"] <- "2"
SRR700080_meta$Passage[SRR700080_meta$Passage == "P1-P3 (after differentiation)"] <- "1"
SRR700080_meta$Passage <- as.numeric(SRR700080_meta$Passage)

SRR700080_meta$treatment[SRR700080_meta$treatment == "non-IR"] <- "none"
SRR700080_meta$treatment[SRR700080_meta$treatment == "IR 5Gy"] <- "IR_5Gy"

SRR700080_meta$source_name[SRR700080_meta$source_name == "iPS cells"] <- "iPS"
SRR700080_meta$source_name[SRR700080_meta$source_name == "Neural progenitor cells"] <- "NPC"

SRR700080_meta$Gene <- paste0("NTg_",SRR700080_meta$source_name,"_",SRR700080_meta$treatment)
SRR700080_meta$Gene[SRR700080_meta$source_name == "Fibroblasts" & SRR700080_meta$treatment == "none"] <- "NTg"

SRR700080_meta <- ddply(SRR700080_meta,c("Passage","source_name","treatment"),transform,Replicate=seq(from=1,by=1,length.out=length(treatment)))

SRR700080_meta <- data.frame(SampleID = SRR700080_meta$Run,
                             GeneID = 0,
                             Replicate = SRR700080_meta$Replicate,
                             Gene = SRR700080_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR700080_meta$source_name,
                             Passage = -1,
                             Line_age = -1,
                             Dox = 0,
                             LIB = "LIBSRR700080",
                             Comment = SRR700080_meta$treatment
)

commonsamps <- intersect(SRR700080_meta$Sample,colnames(SRR700080))
SRR700080_meta <- SRR700080_meta[SRR700080_meta$Sample %in% commonsamps,]
SRR700080 <- SRR700080[,c("Geneid",commonsamps)]

all(SRR700080_meta$Sample == colnames(SRR700080)[-1])

#SRR178756 (mixed/ total RNA)
SRR178756 <- read.table("./Data/15_SRR178756_featurecounts.txt", header = T, sep = "\t")
SRR178756_meta <- read.table("./Data/15_SraRunTable.txt", header = T, sep = ",")

SRR178756[ ,2:6] <- list(NULL)
colnames(SRR178756)[-1] <- str_extract(colnames(SRR178756)[-1],"SRR[0-9]+")
SRR178756$Geneid <- gsub("\\.[0-9]+","",SRR178756$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR178756$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR178756 <- SRR178756[SRR178756$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR178756 <- SRR178756[order(SRR178756$Geneid),]
all(SRR178756$Geneid == ensTolen_agg$Group.1)

SRR178756_meta <- SRR178756_meta[,c("Run","agent","Cell_Line","Genotype","Passages")]
SRR178756_meta$Passages[SRR178756_meta$Passages == "15-20"] <- "15"
SRR178756_meta$Passages <- as.numeric(SRR178756_meta$Passages)

SRR178756_meta$Gene <- paste0("NTg_",SRR178756_meta$agent)
SRR178756_meta$Gene[SRR178756_meta$Cell_Line %in% c("C1","C2") & SRR178756_meta$agent == "DMSO"] <- "NTg"

SRR178756_meta <- ddply(SRR178756_meta,c("Cell_Line","agent","Genotype"),transform,Replicate=seq(from=1,by=1,length.out=length(agent)))

SRR178756_meta <- data.frame(SampleID = SRR178756_meta$Run,
                             GeneID = 0,
                             Replicate = SRR178756_meta$Replicate,
                             Gene = SRR178756_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR178756_meta$Cell_Line,
                             Passage = SRR178756_meta$Passages,
                             Line_age = -1,
                             Dox = 0,
                             LIB = "LIBSRR178756",
                             Comment = paste0(SRR178756_meta$agent,";",SRR178756_meta$Genotype)
)

commonsamps <- intersect(SRR178756_meta$Sample,colnames(SRR178756))
SRR178756_meta <- SRR178756_meta[SRR178756_meta$Sample %in% commonsamps,]
SRR178756 <- SRR178756[,c("Geneid",commonsamps)]

all(SRR178756_meta$Sample == colnames(SRR178756)[-1])

#SRR942655 (dermal fibroblast/ mRNA)
SRR942655 <- read.table("./Data/17_SRR942655_featurecounts.txt", header = T, sep = "\t")
SRR942655_meta <- read.table("./Data/17_SraRunTable.txt", header = T, sep = ",")

SRR942655[ ,2:6] <- list(NULL)
colnames(SRR942655)[-1] <- str_extract(colnames(SRR942655)[-1],"SRR[0-9]+")
SRR942655$Geneid <- gsub("\\.[0-9]+","",SRR942655$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR942655$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR942655 <- SRR942655[SRR942655$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR942655 <- SRR942655[order(SRR942655$Geneid),]
all(SRR942655$Geneid == ensTolen_agg$Group.1)

SRR942655_meta <- SRR942655_meta[,c("Run","source_name","Passage","Donor","Cell_type")]
SRR942655_meta$Passage[SRR942655_meta$Passage == "between 4 and 7"] <- "4"
SRR942655_meta$Passage <- as.numeric(SRR942655_meta$Passage)

SRR942655_meta$Treatment <- SRR942655_meta$source_name
SRR942655_meta$Treatment[SRR942655_meta$Treatment == "NHDF\\, no treatment"] <- "none"
SRR942655_meta$Treatment[SRR942655_meta$Treatment == "NHDF\\, TGFbeta"] <- "Tgfb"
SRR942655_meta$Treatment[SRR942655_meta$Treatment == "NHDF\\, exudate"] <- "exudate"
SRR942655_meta$Treatment[SRR942655_meta$Treatment == "NHDF\\, TGFbeta + exudate"] <- "Tgfb,exudate"


SRR942655_meta$Gene <- paste0("NTg_",SRR942655_meta$Treatment)
SRR942655_meta$Gene[SRR942655_meta$Gene == "NTg_none"] <- "NTg"

SRR942655_meta$Cell_type <- "NHDF"

SRR942655_meta <- ddply(SRR942655_meta,c("Cell_type","Donor","Treatment"),transform,Replicate=seq(from=1,by=1,length.out=length(Treatment)))

SRR942655_meta <- data.frame(SampleID = SRR942655_meta$Run,
                             GeneID = 0,
                             Replicate = SRR942655_meta$Replicate,
                             Gene = SRR942655_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR942655_meta$Cell_type,
                             Passage = SRR942655_meta$Passage,
                             Line_age = -1,
                             Dox = 0,
                             LIB = "LIBSRR942655",
                             Comment = SRR942655_meta$Treatment
)

commonsamps <- intersect(SRR942655_meta$Sample,colnames(SRR942655))
SRR942655_meta <- SRR942655_meta[SRR942655_meta$Sample %in% commonsamps,]
SRR942655 <- SRR942655[,c("Geneid",commonsamps)]

all(SRR942655_meta$Sample == colnames(SRR942655)[-1])

#ERR260891 (dermal fibroblast/ Total (ribo-depleted))
ERR260891 <- read.table("./Data/18_ERR260891_featurecounts.txt", header = T, sep = "\t")
ERR260891_meta <- read.table("./Data/18_SraRunTable.txt", header = T, sep = ",")

ERR260891[ ,2:6] <- list(NULL)
colnames(ERR260891)[-1] <- str_extract(colnames(ERR260891)[-1],"ERR[0-9]+")
ERR260891$Geneid <- gsub("\\.[0-9]+","",ERR260891$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = ERR260891$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
ERR260891 <- ERR260891[ERR260891$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
ERR260891 <- ERR260891[order(ERR260891$Geneid),]
all(ERR260891$Geneid == ensTolen_agg$Group.1)

ERR260891_meta <- ERR260891_meta[,c("Run","sample_name")]

ERR260891_meta$Age <- str_extract(ERR260891_meta$sample_name,"[0-9]+y")
ERR260891_meta$Age <- gsub("y","",ERR260891_meta$Age)
ERR260891_meta$Age <- as.numeric(ERR260891_meta$Age)

ERR260891_meta$Treatment <- str_extract(ERR260891_meta$sample_name,"[a-zA-Z][a-zA-Z]+[0-9]*")
ERR260891_meta$Treatment[is.na(ERR260891_meta$Treatment)] <- "none"

ERR260891_meta$Donor <- str_extract(ERR260891_meta$sample_name,"[0-9]+_")
ERR260891_meta$Donor <- gsub("_","",ERR260891_meta$Donor)

ERR260891_meta$Gene <- paste0("NTg_",ERR260891_meta$Treatment)
ERR260891_meta$Gene[ERR260891_meta$Gene == "NTg_none"] <- "NTg"

ERR260891_meta$Cell_type <- "Fibroblast"

ERR260891_meta <- ddply(ERR260891_meta,c("Donor","Treatment","Age"),transform,Replicate=seq(from=1,by=1,length.out=length(Treatment)))

ERR260891_meta <- data.frame(SampleID = ERR260891_meta$Run,
                             GeneID = 0,
                             Replicate = ERR260891_meta$Replicate,
                             Gene = ERR260891_meta$Gene,
                             gene_symbol = "NTg",
                             Line = ERR260891_meta$Cell_type,
                             Passage = -1,
                             Line_age = ERR260891_meta$Age,
                             Dox = 0,
                             LIB = "LIBERR260891",
                             Comment = ERR260891_meta$Treatment
)

commonsamps <- intersect(ERR260891_meta$Sample,colnames(ERR260891))
ERR260891_meta <- ERR260891_meta[ERR260891_meta$Sample %in% commonsamps,]
ERR260891 <- ERR260891[,c("Geneid",commonsamps)]

all(ERR260891_meta$Sample == colnames(ERR260891)[-1])

#SRR539844 (mixed/ mRNA)
SRR539844 <- read.table("./Data/20_SRR539844_featurecounts.txt", header = T, sep = "\t")
SRR539844_meta <- read.table("./Data/20_SraRunTable.txt", header = T, sep = ",")

SRR539844[ ,2:6] <- list(NULL)
colnames(SRR539844)[-1] <- str_extract(colnames(SRR539844)[-1],"SRR[0-9]+")
SRR539844$Geneid <- gsub("\\.[0-9]+","",SRR539844$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR539844$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR539844 <- SRR539844[SRR539844$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR539844 <- SRR539844[order(SRR539844$Geneid),]
all(SRR539844$Geneid == ensTolen_agg$Group.1)

SRR539844_meta <- SRR539844_meta[,c("Run","Age","Cell_type","gender","tissue")]

SRR539844_meta$Cell_type[SRR539844_meta$Cell_type == "Induced Pluripotent Stem Cells"] <- "iPS"
SRR539844_meta$Cell_type[SRR539844_meta$Cell_type == "Human embryonic stem cells"] <- "ESC"
SRR539844_meta$Cell_type[SRR539844_meta$Cell_type == "fibroblasts"] <- "Fibroblast"

SRR539844_meta$AgeNum <- 0
SRR539844_meta$AgeNum[SRR539844_meta$Age == "50 years old"] <- 50

SRR539844_meta$Gene <- paste0("NTg_",SRR539844_meta$Cell_type)
SRR539844_meta$Gene[SRR539844_meta$Gene == "NTg_Fibroblast"] <- "NTg"

SRR539844_meta <- ddply(SRR539844_meta,c("Cell_type","gender","tissue"),transform,Replicate=seq(from=1,by=1,length.out=length(tissue)))

SRR539844_meta <- data.frame(SampleID = SRR539844_meta$Run,
                             GeneID = 0,
                             Replicate = SRR539844_meta$Replicate,
                             Gene = SRR539844_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR539844_meta$Cell_type,
                             Passage = -1,
                             Line_age = SRR539844_meta$AgeNum,
                             Dox = 0,
                             LIB = "LIBSRR539844",
                             Comment = SRR539844_meta$gender
)

commonsamps <- intersect(SRR539844_meta$Sample,colnames(SRR539844))
SRR539844_meta <- SRR539844_meta[SRR539844_meta$Sample %in% commonsamps,]
SRR539844 <- SRR539844[,c("Geneid",commonsamps)]

all(SRR539844_meta$Sample == colnames(SRR539844)[-1])

#SRR357505 (D551/FM1/ embryo/ total RNA)
SRR357505 <- read.table("./Data/21_SRR357505_featurecounts.txt", header = T, sep = "\t")
SRR357505_meta <- read.table("./Data/21_SraRunTable.txt", header = T, sep = ",")

SRR357505[ ,2:6] <- list(NULL)
colnames(SRR357505)[-1] <- str_extract(colnames(SRR357505)[-1],"SRR[0-9]+")
SRR357505$Geneid <- gsub("\\.[0-9]+","",SRR357505$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR357505$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR357505 <- SRR357505[SRR357505$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR357505 <- SRR357505[order(SRR357505$Geneid),]
all(SRR357505$Geneid == ensTolen_agg$Group.1)

SRR357505_meta <- SRR357505_meta[,c("Run","Assay.Type","Cell_type","source_name","treatment")]
SRR357505_meta <- SRR357505_meta[SRR357505_meta$Assay.Type == "RNA-Seq",]

SRR357505_meta$Cell_type[SRR357505_meta$Cell_type == "FM1-derived iPSC"] <- "iPS"
SRR357505_meta$Cell_type[SRR357505_meta$Cell_type == "D551-derived iPSC"] <- "iPS"

SRR357505_meta$AgeNum <- 0

SRR357505_meta$Gene <- "NTg_Treated"
SRR357505_meta$Gene[SRR357505_meta$treatment == "none"] <- "NTg"

SRR357505_meta <- ddply(SRR357505_meta,c("Cell_type","treatment"),transform,Replicate=seq(from=1,by=1,length.out=length(treatment)))

SRR357505_meta <- data.frame(SampleID = SRR357505_meta$Run,
                             GeneID = 0,
                             Replicate = SRR357505_meta$Replicate,
                             Gene = SRR357505_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR357505_meta$Cell_type,
                             Passage = -1,
                             Line_age = SRR357505_meta$AgeNum,
                             Dox = 0,
                             LIB = "LIBSRR357505",
                             Comment = SRR357505_meta$treatment
)

commonsamps <- intersect(SRR357505_meta$Sample,colnames(SRR357505))
SRR357505_meta <- SRR357505_meta[SRR357505_meta$Sample %in% commonsamps,]
SRR357505 <- SRR357505[,c("Geneid",commonsamps)]

all(SRR357505_meta$Sample == colnames(SRR357505)[-1])

#SRR203849 (mixed/ mRNA) EXCLUDE!!!! PLUS AND MINUS STRAND SEPERATED
SRR203849 <- read.table("./Data/22_SRR203849_featurecounts.txt", header = T, sep = "\t")
SRR203849_meta <- read.table("./Data/22_SraRunTable.txt", header = T, sep = ",")

SRR203849[ ,2:6] <- list(NULL)
colnames(SRR203849)[-1] <- str_extract(colnames(SRR203849)[-1],"SRR[0-9]+")
SRR203849$Geneid <- gsub("\\.[0-9]+","",SRR203849$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR203849$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR203849 <- SRR203849[SRR203849$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR203849 <- SRR203849[order(SRR203849$Geneid),]
all(SRR203849$Geneid == ensTolen_agg$Group.1)

SRR203849_meta <- SRR203849_meta[,c("Run","source_name","treatment","Cell_Line")]

SRR203849_meta$source_name[SRR203849_meta$source_name == "reference hESC line"] <- "ESC"
SRR203849_meta$source_name[SRR203849_meta$source_name == "reference hIPSC line"] <- "iPS"
SRR203849_meta$source_name[SRR203849_meta$source_name == "reprogramming intermediates of hiF-T cells"] <- "reprogInt_hiF-T_cells"
SRR203849_meta$source_name[SRR203849_meta$source_name == "human inducible fibroblasts-like cells (hiF)"] <- "hiF"
SRR203849_meta$source_name[SRR203849_meta$source_name == "immortalized human inducible fibroblasts-like cells (hiF-T)"] <- "hiF-T"
SRR203849_meta$source_name[SRR203849_meta$source_name == "reprogramming intermediates of hiF cells"] <- "reprogInt_hiF_cells"
SRR203849_meta$source_name[SRR203849_meta$source_name == "hIPSCs from reprogrammed hiF cells (hIPSC)"] <- "hiF-iPS"
SRR203849_meta$source_name[SRR203849_meta$source_name == "hIPSCs from reprogrammed hiF-T cells (hIPSC-T)"] <- "hiF-T-iPS"

SRR203849_meta$AgeNum <- 0
SRR203849_meta$AgeNum[!SRR203849_meta$source_name %in% c("ESC","iPS","BJ")] <- -1

SRR203849_meta$treatment[SRR203849_meta$treatment == "no DOX"] <- "none"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-2"] <- "DOX"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-5"] <- "DOX"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-8"] <- "DOX"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-10"] <- "DOX"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-14"] <- "DOX"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-20"] <- "DOX"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-24"] <- "DOX"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-19"] <- "DOX"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-5; LSD1 inhibition"] <- "DOX;LSD1_inhib"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-5; LSD1 control"] <- "DOX;LSD1_control"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-13; LSD1 inhibition"] <- "DOX;LSD1_inhib"
SRR203849_meta$treatment[SRR203849_meta$treatment == "DOX days 0-13; LSD1 control"] <- "DOX;LSD1_control"
SRR203849_meta$treatment[SRR203849_meta$treatment == ""] <- "none"

SRR203849_meta$Gene <- paste0("NTg_",SRR203849_meta$treatment)
SRR203849_meta$Gene[SRR203849_meta$Gene == "NTg_none"] <- "NTg"

SRR203849_meta <- ddply(SRR203849_meta,c("source_name","treatment"),transform,Replicate=seq(from=1,by=1,length.out=length(treatment)))

SRR203849_meta <- data.frame(SampleID = SRR203849_meta$Run,
                             GeneID = 0,
                             Replicate = SRR203849_meta$Replicate,
                             Gene = SRR203849_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR203849_meta$source_name,
                             Passage = -1,
                             Line_age = SRR203849_meta$AgeNum,
                             Dox = 0,
                             LIB = "LIBSRR203849",
                             Comment = SRR203849_meta$treatment
)

commonsamps <- intersect(SRR203849_meta$Sample,colnames(SRR203849))
SRR203849_meta <- SRR203849_meta[SRR203849_meta$Sample %in% commonsamps,]
SRR203849 <- SRR203849[,c("Geneid",commonsamps)]

all(SRR203849_meta$Sample == colnames(SRR203849)[-1])

#SRR180214 (mixed/ total RNA)
SRR180214 <- read.table("./Data/23_SRR180214_featurecounts.txt", header = T, sep = "\t")
SRR180214_meta <- read.table("./Data/23_SraRunTable.txt", header = T, sep = ",")

SRR180214[ ,2:6] <- list(NULL)
colnames(SRR180214)[-1] <- str_extract(colnames(SRR180214)[-1],"SRR[0-9]+")
SRR180214$Geneid <- gsub("\\.[0-9]+","",SRR180214$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR180214$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR180214 <- SRR180214[SRR180214$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR180214 <- SRR180214[order(SRR180214$Geneid),]
all(SRR180214$Geneid == ensTolen_agg$Group.1)

SRR180214_meta <- SRR180214_meta[,c("Run","Assay.Type","Cell_type")]
SRR180214_meta <- SRR180214_meta[SRR180214_meta$Assay.Type == "RNA-Seq",]

SRR180214_meta$Cell_type[SRR180214_meta$Cell_type == "breast epithelial"] <- "Epithelial"

SRR180214_meta$AgeNum <- -1

SRR180214_meta$Gene <- "NTg"

SRR180214_meta <- ddply(SRR180214_meta,c("Cell_type"),transform,Replicate=seq(from=1,by=1,length.out=length(Cell_type)))

SRR180214_meta <- data.frame(SampleID = SRR180214_meta$Run,
                             GeneID = 0,
                             Replicate = SRR180214_meta$Replicate,
                             Gene = SRR180214_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR180214_meta$Cell_type,
                             Passage = -1,
                             Line_age = SRR180214_meta$AgeNum,
                             Dox = 0,
                             LIB = "LIBSRR180214",
                             Comment = ""
)

commonsamps <- intersect(SRR180214_meta$Sample,colnames(SRR180214))
SRR180214_meta <- SRR180214_meta[SRR180214_meta$Sample %in% commonsamps,]
SRR180214 <- SRR180214[,c("Geneid",commonsamps)]

all(SRR180214_meta$Sample == colnames(SRR180214)[-1])

#SRR519079 (mixed/ mRNA)  3-seq???? Maybe exclude
SRR519079 <- read.table("./Data/24_SRR519079_featurecounts.txt", header = T, sep = "\t")
SRR519079_meta <- read.table("./Data/24_SraRunTable.txt", header = T, sep = ",")

SRR519079[ ,2:6] <- list(NULL)
colnames(SRR519079)[-1] <- str_extract(colnames(SRR519079)[-1],"SRR[0-9]+")
SRR519079$Geneid <- gsub("\\.[0-9]+","",SRR519079$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR519079$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR519079 <- SRR519079[SRR519079$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR519079 <- SRR519079[order(SRR519079$Geneid),]
all(SRR519079$Geneid == ensTolen_agg$Group.1)

SRR519079_meta <- SRR519079_meta[,c("Run","age_group","sex","source_name","treatment","PATIENT_ID")]

SRR519079_meta$Cell_type <- "Skin"

SRR519079_meta$AgeNum <- -1

SRR519079_meta$Gene <- paste0("NTg_",SRR519079_meta$treatment)
SRR519079_meta$Gene[SRR519079_meta$Gene %in% c("NTg_none","NTg_")] <- "NTg"

SRR519079_meta <- ddply(SRR519079_meta,c("source_name"),transform,Replicate=seq(from=1,by=1,length.out=length(source_name)))

SRR519079_meta <- data.frame(SampleID = SRR519079_meta$Run,
                             GeneID = 0,
                             Replicate = SRR519079_meta$Replicate,
                             Gene = SRR519079_meta$Gene,
                             gene_symbol = "NTg",
                             Line = SRR519079_meta$Cell_type,
                             Passage = -1,
                             Line_age = -1,
                             Dox = 0,
                             LIB = "LIBSRR519079",
                             Comment = SRR519079_meta$source_name
)

commonsamps <- intersect(SRR519079_meta$Sample,colnames(SRR519079))
SRR519079_meta <- SRR519079_meta[SRR519079_meta$Sample %in% commonsamps,]
SRR519079 <- SRR519079[,c("Geneid",commonsamps)]

all(SRR519079_meta$Sample == colnames(SRR519079)[-1])

#SRR150962 (Interventions/ total RNA)
SRR150962 <- read.table("./Data/3_SRR150962_featurecounts.txt", header = T, sep = "\t")
SRR150962_meta <- read.table("./Data/3_SRR150962_RunTable.txt", header = T, sep = ",")

SRR150962[ ,2:6] <- list(NULL)
colnames(SRR150962)[-1] <- str_extract(colnames(SRR150962)[-1],"SRR[0-9]+")
SRR150962$Geneid <- gsub("\\.[0-9]+","",SRR150962$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR150962$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR150962 <- SRR150962[SRR150962$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR150962 <- SRR150962[order(SRR150962$Geneid),]
all(SRR150962$Geneid == ensTolen_agg$Group.1)

SRR150962_meta <- SRR150962_meta[,c("Run","DONOR_AGE","Cell_Line","clinical_condition","percent_oxygen","replicate_line","inhouse_cell_line_name","Passage","treatments")]
SRR150962_meta$clinical_condition[SRR150962_meta$clinical_condition == "Normal"] <- "NTg"
SRR150962_meta$clinical_condition[SRR150962_meta$clinical_condition == "SURF1_Mutation"] <- "SURF1"
SRR150962_meta$clinical_condition <- paste0(SRR150962_meta$clinical_condition,"_",SRR150962_meta$treatments,"_Ox",SRR150962_meta$percent_oxygen)
SRR150962_meta$clinical_condition[SRR150962_meta$clinical_condition == "NTg_Control_Ox21"] <- "NTg"

SRR150962_meta$Passage[is.na(SRR150962_meta$Passage)] <- 0

SRR150962_meta <- data.frame(SampleID = SRR150962_meta$Run,
                             GeneID = 0,
                             Replicate = SRR150962_meta$replicate_line+1,
                             Gene = SRR150962_meta$clinical_condition,
                             gene_symbol = SRR150962_meta$clinical_condition,
                             Line = SRR150962_meta$inhouse_cell_line_name,
                             Passage = SRR150962_meta$Passage,
                             Line_age = SRR150962_meta$DONOR_AGE,
                             Dox = 0,
                             LIB = "LIBSRR150962",
                             Comment = SRR150962_meta$treatments
)

commonsamps <- intersect(SRR150962_meta$Sample,colnames(SRR150962))
SRR150962_meta <- SRR150962_meta[SRR150962_meta$Sample %in% commonsamps,]
SRR150962 <- SRR150962[,c("Geneid",commonsamps)]

all(SRR150962_meta$Sample == colnames(SRR150962)[-1])

#ERR270009 (mRNA)
ERR270009 <- read.table("./Data/7_ERR270009_featurecounts.txt", header = T, sep = "\t")
ERR270009_meta <- read.table("./Data/7_ERR270009_Run_Table_complete.txt", header = T, sep = "\t")

ERR270009[ ,2:6] <- list(NULL)
colnames(ERR270009)[-1] <- str_extract(colnames(ERR270009)[-1],"ERR[0-9]+")
ERR270009$Geneid <- gsub("\\.[0-9]+","",ERR270009$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = ERR270009$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
ERR270009 <- ERR270009[ERR270009$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
ERR270009 <- ERR270009[order(ERR270009$Geneid),]
all(ERR270009$Geneid == ensTolen_agg$Group.1)

library(plyr)
ERR270009_meta <- ERR270009_meta[,c("Comment.ENA_RUN.","Characteristics.individual.","Characteristics.age.","Factor.Value.compound.","Factor.Value.sampling.time.point.")]
colnames(ERR270009_meta) <- c("Run","Individual","Age","Compound","TimePoint")
ERR270009_meta <- ERR270009_meta[!duplicated(ERR270009_meta),]
ERR270009_meta <- ddply(ERR270009_meta,c("Individual","Compound","TimePoint"),transform,Replicate=seq(from=1,by=1,length.out=length(Individual)))
ERR270009_meta <- ERR270009_meta[ERR270009_meta$Age != "not available",]
ERR270009_meta$Age <- as.numeric(ERR270009_meta$Age)
ERR270009_meta$Gene <- paste0(ERR270009_meta$Compound,"_",ERR270009_meta$TimePoint)
ERR270009_meta$Gene[ERR270009_meta$Gene == "none_0"] <- "NTg"

ERR270009_meta <- data.frame(SampleID = ERR270009_meta$Run,
                             GeneID = 0,
                             Replicate = ERR270009_meta$Replicate,
                             Gene = ERR270009_meta$Gene,
                             gene_symbol = ERR270009_meta$Gene,
                             Line = ERR270009_meta$Individual,
                             Passage = 0,
                             Line_age = ERR270009_meta$Age,
                             Dox = 0,
                             LIB = "LIBERR270009",
                             Comment = ""
)

commonsamps <- intersect(ERR270009_meta$Sample,colnames(ERR270009))
ERR270009_meta <- ERR270009_meta[ERR270009_meta$Sample %in% commonsamps,]
ERR270009 <- ERR270009[,c("Geneid",commonsamps)]

all(ERR270009_meta$Sample == colnames(ERR270009)[-1])

#SRR108366 (mRNA)
SRR108366 <- read.table("./Data/19_SRR108366_featurecounts.txt", header = T, sep = "\t")
SRR108366_meta <- read.table("./Data/19_SRR108366_RunTable.txt", header = T, sep = ",")

SRR108366[ ,2:6] <- list(NULL)
colnames(SRR108366)[-1] <- str_extract(colnames(SRR108366)[-1],"SRR[0-9]+")
SRR108366$Geneid <- gsub("\\.[0-9]+","",SRR108366$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR108366$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR108366 <- SRR108366[SRR108366$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR108366 <- SRR108366[order(SRR108366$Geneid),]
all(SRR108366$Geneid == ensTolen_agg$Group.1)

SRR108366_meta <- SRR108366_meta[,c("Run","Library.Name","Age","tissue")]
SRR108366_meta$Library.Name <- gsub("[0-9]","",SRR108366_meta$Library.Name)
SRR108366_meta <- ddply(SRR108366_meta,"Library.Name",transform,Replicate=seq(from=1,by=1,length.out=length(Library.Name)))
SRR108366_meta$Treated = "NTg"
SRR108366_meta$Treated[grepl("_t",SRR108366_meta$Library.Name)] <- "Treated"
SRR108366_meta$Line <- gsub("_[tuy]","",SRR108366_meta$Library.Name)

SRR108366_meta <- data.frame(SampleID = SRR108366_meta$Run,
                             GeneID = 0,
                             Replicate = 1,
                             Gene = SRR108366_meta$Treated,
                             gene_symbol = SRR108366_meta$Treated,
                             Line = SRR108366_meta$Line,
                             Passage = 0,
                             Line_age = SRR108366_meta$Age,
                             Dox = 0,
                             LIB = "LIBSRR108366",
                             Comment = ""
)

commonsamps <- intersect(SRR108366_meta$Sample,colnames(SRR108366))
SRR108366_meta <- SRR108366_meta[SRR108366_meta$Sample %in% commonsamps,]
SRR108366 <- SRR108366[,c("Geneid",commonsamps)]

all(SRR108366_meta$Sample == colnames(SRR108366)[-1])

#SRR309684 (mRNA)
SRR309684 <- read.table("./Data/26_SRR309684_featurecounts.txt", header = T, sep = "\t")
SRR309684_meta <- read.table("./Data/26_SraRunTable.txt", header = T, sep = ",")

SRR309684[ ,2:6] <- list(NULL)
colnames(SRR309684)[-1] <- str_extract(colnames(SRR309684)[-1],"SRR[0-9]+")
SRR309684$Geneid <- gsub("\\.[0-9]+","",SRR309684$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR309684$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR309684 <- SRR309684[SRR309684$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR309684 <- SRR309684[order(SRR309684$Geneid),]
all(SRR309684$Geneid == ensTolen_agg$Group.1)

SRR309684_meta <- SRR309684_meta[,c("Run","disease","Age","cell_subtype","pmi")]
SRR309684_meta <- ddply(SRR309684_meta,c("Age","disease","cell_subtype","pmi"),transform,Replicate=seq(from=1,by=1,length.out=length(pmi)))
SRR309684_meta$Treated = "NTg_Disease"
SRR309684_meta$Treated[SRR309684_meta$disease == "Not reported"] <- "NTg"

SRR309684_meta <- data.frame(SampleID = SRR309684_meta$Run,
                             GeneID = 0,
                             Replicate = 1,
                             Gene = SRR309684_meta$Treated,
                             gene_symbol = SRR309684_meta$Treated,
                             Line = "Fibroblast",
                             Passage = 0,
                             Line_age = SRR309684_meta$Age,
                             Dox = 0,
                             LIB = "LIBSRR309684",
                             Comment = paste0(SRR309684_meta$disease,";pmi:",SRR309684_meta$pmi)
)

commonsamps <- intersect(SRR309684_meta$Sample,colnames(SRR309684))
SRR309684_meta <- SRR309684_meta[SRR309684_meta$Sample %in% commonsamps,]
SRR309684 <- SRR309684[,c("Geneid",commonsamps)]

all(SRR309684_meta$Sample == colnames(SRR309684)[-1])

#SRR122069 (Total)
SRR122069 <- read.table("./Data/27_SRR122069_featurecounts.txt", header = T, sep = "\t")
SRR122069_meta <- read.table("./Data/27_SraRunTable.txt", header = T, sep = ",")

SRR122069[ ,2:6] <- list(NULL)
colnames(SRR122069)[-1] <- str_extract(colnames(SRR122069)[-1],"SRR[0-9]+")
SRR122069$Geneid <- gsub("\\.[0-9]+","",SRR122069$Geneid)
ensTolen <- getBM(attributes=c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = SRR122069$Geneid, mart= ensembl)
ensTolen_agg <- aggregate(ensTolen$transcript_length,by = list(ensTolen$ensembl_gene_id),FUN = median)
SRR122069 <- SRR122069[SRR122069$Geneid %in% ensTolen_agg$Group.1,]
ensTolen_agg <- ensTolen_agg[order(ensTolen_agg$Group.1),]
SRR122069 <- SRR122069[order(SRR122069$Geneid),]
all(SRR122069$Geneid == ensTolen_agg$Group.1)

SRR122069_meta <- SRR122069_meta[SRR122069_meta$Assay.Type == "RNA-Seq",c("Run","Cell_type","source_name")]
SRR122069_meta$Annot <- c("D0",
                          "D3",
                          "D7",
                          "D11",
                          "D15",
                          "D20",
                          "D28",
                          "p1",
                          "p2",
                          "p3",
                          "iPS",
                          "iPS_TIG108-4f3_clone6",
                          "iPS_TIG118-4f1_clone7",
                          "iPS_TKCBV5-6_clone1",
                          "iPS_451F3_clone5")

SRR122069_meta <- ddply(SRR122069_meta,c("Annot"),transform,Replicate=seq(from=1,by=1,length.out=length(Annot)))
SRR122069_meta$Treated = paste0("NTg_",SRR122069_meta$Annot)
SRR122069_meta$Treated[SRR122069_meta$Treated == "NTg_D0"] <- "NTg"

SRR122069_meta <- data.frame(SampleID = SRR122069_meta$Run,
                             GeneID = 0,
                             Replicate = SRR122069_meta$Replicate,
                             Gene = SRR122069_meta$Treated,
                             gene_symbol = SRR122069_meta$Treated,
                             Line = "Fibroblast",
                             Passage = -1,
                             Line_age = -1,
                             Dox = 0,
                             LIB = "LIBSRR122069",
                             Comment = SRR122069_meta$Annot
)

commonsamps <- intersect(SRR122069_meta$Sample,colnames(SRR122069))
SRR122069_meta <- SRR122069_meta[SRR122069_meta$Sample %in% commonsamps,]
SRR122069 <- SRR122069[,c("Geneid",commonsamps)]

all(SRR122069_meta$Sample == colnames(SRR122069)[-1])

#Exclude SRR700080,SRR519079,SRR203849 //SRR539844,SRR357505

#Combine all new data (ASL, SRRs and ERRs)
meta_comb_new <- rbind.data.frame(ASL_meta,
                                  ERR137716_meta,
                                  ERR260891_meta,
                                  ERR270009_meta,
                                  SRR101478_meta,
                                  SRR108366_meta,
                                  SRR121037_meta,
                                  SRR150962_meta,
                                  SRR173636_meta,
                                  SRR176804_meta,
                                  SRR178756_meta,
                                  SRR180214_meta,
                                  #SRR203849_meta,
                                  SRR275112_meta,
                                  SRR303074_meta,
                                  SRR315202_meta,
                                  SRR357505_meta,
                                  #SRR519079_meta,
                                  SRR539844_meta,
                                  SRR627541_meta,
                                  SRR700080_meta,
                                  SRR810099_meta,
                                  SRR942655_meta,
                                  SRR309684_meta,
                                  SRR122069_meta)
meta_comb_new$LIBPrep <- "mRNA"
meta_comb_new$LIBPrep[meta_comb_new$LIB %in% c("LIBSRR122069","LIB047784","LIB048130","LIB049469","LIB050518","LIBSRR315202","LIBSRR275112","LIBSRR101478","LIBSRR173636","LIBSRR700080","LIBSRR178756","LIBERR260891","LIBSRR357505","LIBSRR180214","LIBSRR150962")] <- "Total"
meta_comb_new$Dox[is.na(meta_comb_new$Dox)] <- 0

meta_comb_new$CT <- "Fibroblast"
meta_comb_new$CT[meta_comb_new$Line %in% c("F3","B6")] <- "iPS"
meta_comb_new$CT[meta_comb_new$Line %in% c("EC")] <- "EC"
meta_comb_new$CT[meta_comb_new$Line %in% c("Skin","Skin_EarlyTime","Skin_LateTime")] <- "Skin"

samples_new_list <- list(ASL_20M,
                         ERR137716,
                         ERR260891,
                         ERR270009,
                         SRR101478,
                         SRR108366,
                         SRR121037,
                         SRR150962,
                         SRR173636,
                         SRR176804,
                         SRR178756,
                         SRR180214,
                         #SRR203849,
                         SRR275112,
                         SRR303074,
                         SRR315202,
                         SRR357505,
                         #SRR519079,
                         SRR539844,
                         SRR627541,
                         SRR700080,
                         SRR810099,
                         SRR942655,
                         SRR309684,
                         SRR122069)
samples_new_df <- Reduce(function(x, y) merge(x, y, by = "Geneid"), samples_new_list)
samples_new_df[1:5,1:5]

meta_comb_new$SampleID <- colnames(samples_new_df)[-1]

all(colnames(samples_new_df)[-1] == meta_comb_new$SampleID)

meta_comb_new$Passage[is.na(meta_comb_new$Passage)] <- 1


#idx_wt <- which(colnames(samples_new_df) %in% meta_comb_new$SampleID[which(meta_comb_new$Gene == "NTg" & (meta_comb_new$Line_age > -1) & (! meta_comb_new$Line %in% excludedWTLines) & (meta_comb_new$Passage <= maxPassagingNumber))])
#wt_libs <- unique(meta_comb_new$LIB[meta_comb_new$SampleID %in% colnames(samples_new_df)[idx_wt]])


sort(unique(wt_libs))
##### TEST #####
meta_comb_new <- meta_comb_new[meta_comb_new$LIB %in% c("LIBSRR122069"),]
#meta_comb_new <- meta_comb_new[startsWith(meta_comb_new$SampleID,"C"),]
samples_new_df <- samples_new_df[,c("Geneid",meta_comb_new$SampleID)]

excludedWTLines <- c("ESC","iPS","H1","H9","EC","Epithelial","rapamycin","Parkin-expressing MRC5 fibroblasts","Skin_EarlyTime","Skin_LateTime","prolif","quiescent","","NPC")
maxPassagingNumber <- 20

#Now finalize the data
idx_wt <- which(colnames(samples_new_df) %in% meta_comb_new$SampleID[which(meta_comb_new$Gene == "NTg" & (meta_comb_new$Line_age > -1) & (! meta_comb_new$Line %in% excludedWTLines) & (meta_comb_new$Passage <= maxPassagingNumber))])
idx_bfp <- which(colnames(samples_new_df) %in% meta_comb_new$SampleID[which(meta_comb_new$Gene == "mTagBFP2" & (meta_comb_new$Line_age > -1) & (! meta_comb_new$Line %in% excludedWTLines) & (meta_comb_new$Passage <= maxPassagingNumber))])
idx_test <- which(colnames(samples_new_df) %in% meta_comb_new$SampleID[which(!(meta_comb_new$Gene %in% c("mTagBFP2","NTg")) | meta_comb_new$Line_age == -1 | (meta_comb_new$Line %in% excludedWTLines) | (meta_comb_new$Passage > maxPassagingNumber))])
data_merged <- merge(train_data[,c(1,which(colnames(train_data) %in% train_metadata$Sample[which(train_metadata$Disease == "Normal")]))],samples_new_df[,c(1,idx_wt,idx_bfp,idx_test)], by = "Geneid")

asl_ctrl_batches <- meta_comb_new[match(colnames(samples_new_df[c(idx_wt, idx_bfp)]), meta_comb_new$SampleID), "LIB"] #get the batches for the asl control samples (wt,bfp)
asl_batches <- meta_comb_new[match(colnames(samples_new_df[idx_test]), meta_comb_new$SampleID), "LIB"] #get batches for the asl test samples
asl_ctrl_ages <- meta_comb_new[match(colnames(samples_new_df[c(idx_wt, idx_bfp)]), meta_comb_new$SampleID), "Line_age"] #get ages for asl control samples
asl_ages <- meta_comb_new[match(colnames(samples_new_df[idx_test]), meta_comb_new$SampleID), "Line_age"] #get ages for asl control samples
asl_ctrl_libprep <- meta_comb_new[match(colnames(samples_new_df[c(idx_wt, idx_bfp)]), meta_comb_new$SampleID), "LIBPrep"]
asl_libprep <- meta_comb_new[match(colnames(samples_new_df[idx_test]), meta_comb_new$SampleID), "LIBPrep"] #get batches for the asl test samples
asl_ctrl_dox <- meta_comb_new[match(colnames(samples_new_df[c(idx_wt, idx_bfp)]), meta_comb_new$SampleID), "Dox"]
asl_dox <- meta_comb_new[match(colnames(samples_new_df[idx_test]), meta_comb_new$SampleID), "Dox"] #get batches for the asl test samples
asl_ctrl_line <- meta_comb_new[match(colnames(samples_new_df[c(idx_wt, idx_bfp)]), meta_comb_new$SampleID), "CT"]
asl_line <- meta_comb_new[match(colnames(samples_new_df[idx_test]), meta_comb_new$SampleID), "CT"] #get batches for the asl test samples
asl_ctrl_gene <- meta_comb_new[match(colnames(samples_new_df[c(idx_wt, idx_bfp)]), meta_comb_new$SampleID), "Gene"]
asl_gene <- meta_comb_new[match(colnames(samples_new_df[idx_test]), meta_comb_new$SampleID), "Gene"] #get batches for the asl test samples
asl_ctrl_passage <- meta_comb_new[match(colnames(samples_new_df[c(idx_wt, idx_bfp)]), meta_comb_new$SampleID), "Passage"]
asl_passage <- meta_comb_new[match(colnames(samples_new_df[idx_test]), meta_comb_new$SampleID), "Passage"] #get batches for the asl test samples
  
  
  
