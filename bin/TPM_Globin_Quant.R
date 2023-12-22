# Load in the packages

library(tidyverse)
library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
# Read in the counts
bovine_counts <- fread(args[1])
reindeer_counts <- fread(args[2])
human_counts <- fread(args[3])


# Read in annotation
bovine_annotation <- fread(args[4]) %>% filter(V3 == "gene") %>% 
select(V4,V5,V9) %>% separate(., V9, into = c("ensemble", "gene_version", "gene_name"), sep = ";")
bovine_annotation$ensemble <- gsub("gene_id ", "", bovine_annotation$ensemble)
bovine_annotation$ensemble <- gsub('"', "", bovine_annotation$ensemble)
bovine_annotation$ensemble <- gsub(' ', "", bovine_annotation$ensemble)
bovine_annotation$gene_name <- gsub("gene_name ", "", bovine_annotation$gene_name)
bovine_annotation$gene_name <- gsub("gene_source ", "", bovine_annotation$gene_name)
bovine_annotation$gene_name <- gsub('"', "", bovine_annotation$gene_name)
bovine_annotation$gene_name <- gsub(' ', "", bovine_annotation$gene_name)
bovine_annotation$length = abs(bovine_annotation$V4 - bovine_annotation$V5)
bovine_annotation <- bovine_annotation %>% filter(ensemble %in% bovine_counts$Geneid)
bovine_annotation <- bovine_annotation %>% select(3,5,6)


# Reindeer annotation
reindeer_annotation <- fread(args[5]) %>% filter(V3 == "gene") %>%
select(V4,V5,V9) %>% separate(., V9, into = c("gene_id", "id", "name", "annotation"), sep = ";")
reindeer_annotation$length <- abs(reindeer_annotation$V4 - reindeer_annotation$V5)
reindeer_annotation$gene_id <- gsub("gene_id ", "", reindeer_annotation$gene_id)
reindeer_annotation$gene_id <- gsub('"', "", reindeer_annotation$gene_id)
reindeer_annotation$gene_id<- gsub(' ', "", reindeer_annotation$gene_id)
reindeer_annotation <- reindeer_annotation %>% select(gene_id, length, annotation)


# Human annotation
human_annotation <- fread(args[6])
human_annotation <- human_annotation %>% filter(V3 == "gene")
human_annotation <- human_annotation %>% separate(., V9, into = c("gene_id", "gene_type", "gene_name", "level"), sep = ";") %>% select(1:11)
human_annotation$gene_name <- gsub("gene_name ", "", human_annotation$gene_name)
human_annotation$gene_id <- gsub("gene_id ", "", human_annotation$gene_id)
human_annotation$gene_id <- gsub('"', "", human_annotation$gene_id)
human_annotation$length <- abs(human_annotation$V4 - human_annotation$V5)
human_annotation <- human_annotation %>% select(gene_id, length)


# Now join the counts with annotation
bovine_pre_TPM <- bovine_counts
rownames(bovine_pre_TPM) <- bovine_pre_TPM$Geneid
bovine_pre_TPM <- bovine_pre_TPM %>% select(-1)
rownames(bovine_pre_TPM)
dim(bovine_annotation)
dim(bovine_counts)
all(rownames(bovine_pre_TPM) == bovine_annotation$ensemble)
bovine_annotation_simple <- bovine_annotation %>% select(-c(1,2)) %>% as.matrix()
bovine_pre_TPM = cbind(bovine_pre_TPM, bovine_annotation_simple)
bovine_pre_TPM <- as.matrix(bovine_pre_TPM)
bovine_TPM <- bovine_pre_TPM[,c(-length(colnames(bovine_pre_TPM)))] / bovine_pre_TPM[,length(colnames(bovine_pre_TPM))]
bovine_TPM <- t( t(bovine_TPM) * 1e6 / colSums(bovine_TPM) )

head(bovine_TPM)
bovine_TPM <- cbind(bovine_TPM, bovine_annotation) %>% select(1:6)
head(bovine_TPM)
rownames(bovine_TPM) <- bovine_counts$Geneid

head(human_annotation)
# Human 
human_annotation_simple <- human_annotation %>% select(gene_id, length)
rownames(human_annotation_simple) <- human_annotation_simple$gene_id
head(human_annotation_simple)
GeneSym <- human_annotation_simple$gene_id
human_annotation_simple <- human_annotation_simple %>% as.matrix()
human_pre_TPM <- cbind(human_counts, human_annotation_simple)
head(human_pre_TPM)

human_pre_TPM <- human_pre_TPM %>% select(-c(Geneid, gene_id))
#human_pre_TPM = cbind(human_pre_TPM, human_annotation_simple)
human_pre_TPM <- sapply(human_pre_TPM, function(x) as.numeric(as.character(x)))
human_pre_TPM <- as.matrix(human_pre_TPM)
rownames(human_pre_TPM) <- GeneSym
human_TPM <- human_pre_TPM[,c(-length(colnames(human_pre_TPM)))] / human_pre_TPM[,length(colnames(human_pre_TPM))]
human_TPM <- t( t(human_TPM) * 1e6 / colSums(human_TPM) )
human_TPM <- as.data.frame(human_TPM)
human_TPM <- cbind(human_TPM, human_annotation_simple)
human_TPM$Gene_id <- rownames(human_TPM)
head(human_TPM)


# Reindeer
reindeer_annotation_simple <- reindeer_annotation %>% select(gene_id, length)
rownames(reindeer_annotation_simple) <- reindeer_annotation_simple$gene_id
head(reindeer_annotation_simple)
GeneSym <- reindeer_annotation_simple$gene_id
reindeer_annotation_simple <- reindeer_annotation_simple %>% as.matrix()
reindeer_pre_TPM <- cbind(reindeer_counts, reindeer_annotation_simple)
head(reindeer_pre_TPM)
reindeer_pre_TPM <- reindeer_pre_TPM %>% select(-c(Geneid, gene_id))
reindeer_pre_TPM <- sapply(reindeer_pre_TPM, function(x) as.numeric(as.character(x)))
reindeer_pre_TPM <- as.matrix(reindeer_pre_TPM)
rownames(reindeer_pre_TPM) <- GeneSym
reindeer_TPM <- reindeer_pre_TPM[,c(-length(colnames(reindeer_pre_TPM)))] / reindeer_pre_TPM[,length(colnames(reindeer_pre_TPM))]
reindeer_TPM <- t( t(reindeer_TPM) * 1e6 / colSums(reindeer_TPM) )
reindeer_TPM <- as.data.frame(reindeer_TPM)
reindeer_TPM <- cbind(reindeer_TPM, reindeer_annotation)
reindeer_TPM$Gene_id <- rownames(reindeer_TPM)
head(reindeer_TPM)
head(reindeer_annotation)

# Now get final data frames in order
reindeer_TPM <- reindeer_TPM %>% select(1:7)
human_TPM <- human_TPM %>% select(1:5)
# bovine in order as ensemble hence why extra column

reindeer_genes <- c("ANN11336","ANN11337","ANN11338","ANN11339","ANN11340","ANN28434","ANN28435","ANN28436", "ANN29413", "ANN29414")
bovine_genes <- c("ENSBTAG00000037815",
"ENSBTAG00000037751",
"ENSBTAG00000037644",
"ENSBTAG00000052199",
"ENSBTAG00000038748",
"ENSBTAG00000047356",
"ENSBTAG00000038921",
"ENSBTAG00000038757",
"ENSBTAG00000054844",
"ENSBTAG00000051412",
"ENSBTAG00000016567")
human_genes <- c("HBA1", "HBA2", "HBB") #HBA1, HBA2, HBB





head(bovine_TPM)
rownames(bovine_TPM)


# bovine
tpm.bov <- bovine_TPM[bovine_genes, -c(5:6)]
tpm.bov
interest_matrix_bov_sum <- colSums(tpm.bov)
interest_matrix_bov_sum
proportion_bov <-  interest_matrix_bov_sum / colSums(bovine_TPM[,1:4])
proportion_bov
percentage_bov <- proportion_bov * 100
percentage_bov <- percentage_bov[-5]
per_df_bov <- as.data.frame(percentage_bov)
colnames(per_df_bov) <- "percentage"
per_df_bov$type <- "Bovine"


# Reindeer
interest_matrix <- reindeer_TPM[reindeer_genes,]
interest_matrix_sum <- colSums(interest_matrix[,1:4])

proportion <- interest_matrix_sum / colSums(reindeer_TPM[,1:4])
percentage <- proportion * 100
per_df <- as.data.frame(percentage)
per_df$type <- "Reindeer"



# human
tpm.hs <- human_TPM[human_genes, -5]
tpm.hs
interest_matrix_hs_sum <- colSums(tpm.hs)
interest_matrix_hs_sum
proportion_hs <-  interest_matrix_hs_sum / colSums(human_TPM[,1:4])
proportion_hs
percentage_hs <- proportion_hs * 100
per_df_hs <- as.data.frame(percentage_hs)
colnames(per_df_hs) <- "percentage"
proportion_hs
per_df_hs
per_df_hs$type <- "Human"



master_df <- rbind(per_df,per_df_bov)
master_df$type <- factor(master_df$type, levels = c("Bovine", "Reindeer"))
ggplot(master_df, aes(x = type, y = percentage, fill = type)) +
  geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.00), size = 3) +
  labs(x = "Species", y = "Percentage (%) of reads mapping to haemoglobin related genes",
       fill = "") +
       theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major.y = element_line(color = "gray", linetype = "dashed")) # add horizontal gridlines
ggsave(args[7], width = 12, height = 12, dpi = 600)

master_df <- rbind(per_df, per_df_bov, per_df_hs)
master_df$type <- factor(master_df$type, levels = c("Bovine", "Reindeer", "Human"))
ggplot(master_df, aes(x = type, y = percentage, fill = type)) +
  geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.00), size = 3) +
  labs(x = "Species", y = "Percentage (%) of reads mapping to haemoglobin related genes",
       fill = "") +
       theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major.y = element_line(color = "gray", linetype = "dashed")) # add horizontal gridlines
ggsave(args[8], width = 12, height = 12, dpi = 600)