###############################
#' TODO: validation of explanatory factor prioritized genomic variations
#'  
###############################################################################
# Working directory
setwd("/~")


###
# 1. Data preparation --------------------
#' 
#'###############################################################################
### (1) Explanatory factor rank
library(xlsx)

factor_rank_lists <- read.xlsx("./Results/Explanatory factor rank.xlsx", sheetIndex = 1, check.names = FALSE)
dim(factor_rank_lists) # 691 8  
factor_rank_lists$final_rank <- gsub("\\_.*", "", factor_rank_lists$'log_weight_AVG_SUM&WEIGHT')
factor_rank_lists$final_rank_class <- gsub(".*\\_", "", factor_rank_lists$'log_weight_AVG_SUM&WEIGHT')
table(factor_rank_lists$final_rank_class)


### (2) Training set sample ID
load(file = "./Results/training_cancer_types.RData")


### (3) mutation status matrix of quality-controlled solid tumor RNA-seq samples 
# 8223 common samples with RNA-seq
load(file = "./Data/TCGA/Solid_tumor_mutation_status_mat.RData") # 487 8223


### (4) GISTIC 2.0 CNV matrix of quality-controlled solid tumor RNA-seq samples
GISTIC_results <- read.csv("./Data/TCGA/all_thresholded.by_genes_whitelisted.tsv", sep = "\t", check.names = FALSE) # 25128 10716

# driver CNV matrix
GISTIC_results <- GISTIC_results[which(GISTIC_results$'Gene Symbol' %in% rownames(mutation_status_mat)),] # 487 10176
colnames(GISTIC_results)[4:10716] <- substr(colnames(GISTIC_results)[4:10716], 1, 16)
rownames(GISTIC_results) <- GISTIC_results$'Gene Symbol'

# common samples
common_samples <- Reduce(intersect, list(training_cancer_types$ID, colnames(mutation_status_mat), colnames(GISTIC_results))) # 7972
GISTIC_results <- GISTIC_results[, common_samples] # 487 7972
mutation_status_mat <- mutation_status_mat[, common_samples] # 487 7972


### (5) Association of driver SNVs and CNVs with CD8A expression
load(file = "./Results/SNV_CD8A_wilcox_res.RData")
load(file = "./Results/CNV_CD8A_wilcox_res.RData")



###
# 2. Mutation frequency comparison --------------------
#' top-ranked and bottom-ranked
#'###############################################################################
#'######
# (1) pancancer SNV mutation frequency --------------------
#'######
storage.mode(mutation_status_mat) <- "numeric"
mutation_freq <- rowSums(mutation_status_mat)/ncol(mutation_status_mat)
pancancer_SNV_freq <- data.frame(Gene = names(mutation_freq), SNV_freq = mutation_freq)
save(pancancer_SNV_freq, file = "./Results/pancancer_SNV_freq.RData")


#'######
# (2) pancancer CNV mutation frequency --------------------
#'######
# class(GISTIC_results); mode(GISTIC_results)
mutation_freq <- apply(GISTIC_results, 1, function(x){
  tmp <- length(which(x != 0))/length(x)
  return(tmp)
})
pancancer_CNV_freq <- data.frame(Gene = names(mutation_freq), CNV_freq = mutation_freq)
save(pancancer_CNV_freq, file = "./Results/pancancer_CNV_freq.RData")



#'######
# (3) Mutation frequency comparison --------------------
#'######
top_ranked_SNV <- factor_rank_lists$final_rank[which(factor_rank_lists$final_rank_class == "SNV")[1:50]]
bottom_ranked_SNV <- factor_rank_lists$final_rank[which(factor_rank_lists$final_rank_class == "SNV")[(204-49):204]]

pancancer_SNV_freq$final_rank <- NA
pancancer_SNV_freq$final_rank[which(pancancer_SNV_freq$Gene %in% top_ranked_SNV)] <- "Top50"
pancancer_SNV_freq$final_rank[which(pancancer_SNV_freq$Gene %in% bottom_ranked_SNV)] <- "Bottom50"
sub_pancancer_SNV_freq <- na.omit(pancancer_SNV_freq)
sub_pancancer_SNV_freq$final_rank <- factor(sub_pancancer_SNV_freq$final_rank, levels = c("Top50", "Bottom50"))

  
## boxplot
library(ggpubr)
library(patchwork)
p_snv <- ggviolin(sub_pancancer_SNV_freq, x = "final_rank", y = "SNV_freq", color = "final_rank",
              fill = "final_rank", 
              add = "boxplot", add.params = list(fill = "white", color = "black")) + 
  stat_compare_means() + # Add global p-value
  labs(y = "SNV Mutation frequency", x = "") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.5), legend.position = "none")




# CNV
top_ranked_CNV <- factor_rank_lists$final_rank[which(factor_rank_lists$final_rank_class == "CNV")[1:50]]
bottom_ranked_CNV <- factor_rank_lists$final_rank[which(factor_rank_lists$final_rank_class == "CNV")[(487-49):487]]

pancancer_CNV_freq$final_rank <- NA
pancancer_CNV_freq$final_rank[which(pancancer_CNV_freq$Gene %in% top_ranked_CNV)] <- "Top50"
pancancer_CNV_freq$final_rank[which(pancancer_CNV_freq$Gene %in% bottom_ranked_CNV)] <- "Bottom50"
sub_pancancer_CNV_freq <- na.omit(pancancer_CNV_freq)
sub_pancancer_CNV_freq$final_rank <- factor(sub_pancancer_CNV_freq$final_rank, levels = c("Top50", "Bottom50"))


p_cnv <- ggviolin(sub_pancancer_CNV_freq, x = "final_rank", y = "CNV_freq", color = "final_rank",
                  fill = "final_rank", 
                  add = "boxplot", add.params = list(fill = "white", color = "black")) + 
  stat_compare_means() + # Add global p-value
  labs(y = "CNV Mutation frequency", x = "") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.5), legend.position = "none")



###
# 3. Association of driver SNVs and CNVs with CD8A expression --------------------
#' top-ranked and bottom-ranked
#'###############################################################################
idx_top <- which(SNV_CD8A_wilcox_res$SNV_gene %in% top_ranked_SNV)
idx_bottom <- which(SNV_CD8A_wilcox_res$SNV_gene %in% bottom_ranked_SNV)

length(which(SNV_CD8A_wilcox_res$'p.adj'[idx_top] < 0.05)) # 44 
length(which(SNV_CD8A_wilcox_res$'p.adj'[idx_bottom] < 0.05)) # 26


top_ranked_CNV <- gsub("-", "_", top_ranked_CNV)
idx_top <- which(CNV_CD8A_wilcox_res$CNV_gene %in% top_ranked_CNV)
idx_bottom <- which(CNV_CD8A_wilcox_res$CNV_gene %in% bottom_ranked_CNV)

length(which(CNV_CD8A_wilcox_res$'p.adj'[idx_top] < 0.05)) # 32 
length(which(CNV_CD8A_wilcox_res$'p.adj'[idx_bottom] < 0.05)) # 25


# Mix SNV and CNV
top_CD8A_wilcox_res <- data.frame(t(sapply(1:50, function(i){
  if(factor_rank_lists$final_rank_class[i] == "CNV"){
    tmp <- CNV_CD8A_wilcox_res[match(factor_rank_lists$final_rank[i], CNV_CD8A_wilcox_res$CNV_gene), ]
  }else{
    tmp <- SNV_CD8A_wilcox_res[match(factor_rank_lists$final_rank[i], SNV_CD8A_wilcox_res$SNV_gene), ]
  }
  colnames(tmp)[1] <- "Gene"
  return(tmp)
})))
length(which(top_CD8A_wilcox_res$'p.adj' < 0.05))# 40 

bottom_CD8A_wilcox_res <- data.frame(t(sapply((nrow(factor_rank_lists) - 49):nrow(factor_rank_lists), function(i){
  if(factor_rank_lists$final_rank_class[i] == "CNV"){
    tmp <- CNV_CD8A_wilcox_res[match(factor_rank_lists$final_rank[i], CNV_CD8A_wilcox_res$CNV_gene), ]
  }else{
    tmp <- SNV_CD8A_wilcox_res[match(factor_rank_lists$final_rank[i], SNV_CD8A_wilcox_res$SNV_gene), ]
  }
  colnames(tmp)[1] <- "Gene"
  return(tmp)
})))

length(which(bottom_CD8A_wilcox_res$'p.adj' < 0.05)) # 25


### Comparison of significance proportion, Fisher's exact test
sum_top <- as.data.frame(table(top_CD8A_wilcox_res$'p.adj' < 0.05))
sum_top$Rank_class <- "Top50"

sum_bottom <- as.data.frame(table(bottom_CD8A_wilcox_res$'p.adj' < 0.05))
sum_bottom$Rank_class <- "Bottom50"

CD8A_cor_count <- rbind.data.frame(sum_top, sum_bottom)
colnames(CD8A_cor_count)[1] <- "Signif"
CD8A_cor_count$Signif <- ifelse(CD8A_cor_count$Signif == TRUE, "adj.p < 0.05", "adj.p > 0.05")

CD8A_cor_count_fisher <- chisq.test(matrix(CD8A_cor_count$Freq, nrow = 2, byrow = TRUE), correct = TRUE)


CD8A_cor_count$percent <- CD8A_cor_count$Freq/rep(tapply(CD8A_cor_count$Freq, CD8A_cor_count$Rank_class, sum), each = 2)
CD8A_cor_count$Rank_class <- factor(CD8A_cor_count$Rank_class, levels = c("Top50", "Bottom50"))

p_CD8_ratio_cor <- ggplot(data = CD8A_cor_count, aes(x = Rank_class, y = percent, fill = Signif)) +  
  geom_bar(stat = "identity") +
  theme_set(theme_classic2()) + 
  labs(title = paste0("Fisher's exact test, p = ", signif(CD8A_cor_count_fisher$p.value, 2))) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 8), plot.title = element_text(size = 10)) +
  xlab("")




# Results Visualization
library(patchwork)
pdf(file = "./Results/Comb Mutation frequency comparison.pdf", height = 4, width = 10)
p_snv + p_cnv + p_CD8_ratio_cor
dev.off()
