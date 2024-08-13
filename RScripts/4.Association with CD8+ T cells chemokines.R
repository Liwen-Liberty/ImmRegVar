###############################
#' TODO: Calculation of association of driver SNVs/CNVs with CD8+ T cells chemokines
#' 
###############################################################################
# Working directory
setwd("/~")

###
# 1. Data preparation --------------------
#' 
#'###############################################################################
### (1) TCGA solid tumor expression
PanCancer.exp <- readRDS(file = "./Data/TCGA/PanCancer.exp.rds") # 17957 10176
colnames(PanCancer.exp) <- substr(colnames(PanCancer.exp), 1, 16)

### (2) Training set sample ID
load(file = "./Results/training_cancer_types.RData")


### (3) mutation status matrix of quality-controlled solid tumor RNA-seq samples 
# 8223 common samples with RNA-seq
load(file = "./Data/TCGA/Solid_tumor_mutation_status_mat.RData") # 487 8223


### (4) GISTIC 2.0 CNV matrix of quality-controlled solid tumor RNA-seq samples
GISTIC_results <- read.csv("./Data/TCGA/all_thresholded.by_genes_whitelisted.tsv", 
                           sep = "\t", check.names = FALSE) # 25128 10716

# driver CNV matrix
GISTIC_results <- GISTIC_results[which(GISTIC_results$'Gene Symbol' %in% rownames(mutation_status_mat)),] # 487 10176
colnames(GISTIC_results)[4:10716] <- substr(colnames(GISTIC_results)[4:10716], 1, 16)
rownames(GISTIC_results) <- GISTIC_results$'Gene Symbol'

# common samples
common_samples <- Reduce(intersect, list(training_cancer_types$ID, colnames(mutation_status_mat), colnames(GISTIC_results))) # 7972
GISTIC_results <- GISTIC_results[, common_samples] # 487 7972
mutation_status_mat <- mutation_status_mat[, common_samples] # 487 7972
PanCancer.exp <- PanCancer.exp[, common_samples]

# CNV class
GISTIC_results_refine <- sapply(GISTIC_results, function(x) {
  res <- gsub("1|2", "Amplified", (gsub("-1|-2", "Deleted", x)))
  return(res)}
)
rownames(GISTIC_results_refine) <- rownames(GISTIC_results) 


### (5) Explanatory factor rank
library(xlsx)

factor_rank_lists <- read.xlsx("./Results/Explanatory factor rank.xlsx", sheetIndex = 1, check.names = FALSE)
dim(factor_rank_lists) # 691 8  
factor_rank_lists$final_rank <- gsub("\\_.*", "", factor_rank_lists$'log_weight_AVG_SUM&WEIGHT')
factor_rank_lists$final_rank_class <- gsub(".*\\_", "", factor_rank_lists$'log_weight_AVG_SUM&WEIGHT')


chemokines <- c("CXCL9", "CXCL10", "CCL5", "CCL3")
all(is.element(chemokines, rownames(PanCancer.exp))) #TRUE



###
# 2. Calculation of association of driver SNVs/CNVs with CD8+ T cells chemokines --------------------
#' 
#'###############################################################################
library(coin)

top_idx <- 1:50
bottom_idx <- (nrow(factor_rank_lists) - 49): nrow(factor_rank_lists)

comb_chemokines_wilcox_res_list <- lapply(list(Top50 = top_idx, Bottom50 = bottom_idx), function(idx){
  chemokines_wilcox_res_list <- lapply(idx, function(i){
    candidate_gene <- factor_rank_lists$final_rank[i]
    
    if(factor_rank_lists$final_rank_class[i] == "SNV"){
      data_merge <- data.frame(t(PanCancer.exp[chemokines, ]), 
                               candidate_gene = mutation_status_mat[candidate_gene, ])
      data_merge$candidate_gene <- factor(data_merge$candidate_gene, levels = c("1", "0"))
      
    }else{
      data_merge <- data.frame(t(PanCancer.exp[chemokines, ]), 
                               candidate_gene = GISTIC_results_refine[candidate_gene, ])
      
      tmp <- sort(table(data_merge$candidate_gene))
      CNV_class <- names(tmp)[2:3]
      if(!is.element("0", CNV_class)) {CNV_class[1] <- "0"} 
      
      data_merge <- subset(data_merge, data_merge[, "candidate_gene"] %in% CNV_class)
      data_merge$candidate_gene <- factor(data_merge$candidate_gene, levels = c(setdiff(CNV_class, "0"), "0"))
      
    }
    
    # each chemokine
    wilcox_res <- lapply(chemokines, function(k){
      tmp_result <- wilcox_test(as.formula(paste0(k, " ~ candidate_gene")), data = data_merge)
      refine.res <- data.frame(Gene = candidate_gene, Chemokine = k, feature_class = factor_rank_lists$final_rank_class[i],
                               Zvalue = statistic(tmp_result, type = "standardized"), p.value = pvalue(tmp_result))
      return(refine.res)
    })
    
    
    chemokine_wilcox_res <- do.call(rbind.data.frame, wilcox_res)
    colnames(chemokine_wilcox_res)[4] <- "Zvalue"
    chemokine_wilcox_res$direction <- gsub("[0-9]", "", rownames(chemokine_wilcox_res))
    
    return(chemokine_wilcox_res)
  })
  comb_chemokine_wilcox_res <- do.call(rbind.data.frame, chemokines_wilcox_res_list)
  comb_chemokine_wilcox_res$p.adj <- p.adjust(comb_chemokine_wilcox_res$p.value, method = "BH")
  return(comb_chemokine_wilcox_res)
})

save(comb_chemokines_wilcox_res_list,  file = "./Results/comb_chemokines_wilcox_res_list.RData")
length(which(comb_chemokines_wilcox_res_list$Top50$'p.adj' < 0.05)) # 172
length(which(comb_chemokines_wilcox_res_list$Bottom50$'p.adj' < 0.05)) # 164


signif_num_list <- lapply(list(Top50 = top_idx, Bottom50 = bottom_idx), function(idx){
  signif_num <- sapply(idx, function(i){
    candidate_gene <- factor_rank_lists$final_rank[i]
    tmp <- subset(comb_chemokines_wilcox_res_list$Top50, Gene == candidate_gene)
    return(length(which(tmp$'p.adj' < 0.05)))
  })
  return(signif_num)
}) 



###
# 3.  Results Visualization --------------------
#' 
#'###############################################################################
library(ggpubr)
signif_num_df <- rbind.data.frame(data.frame(Cor_num = signif_num_list$Top50, 
                                             Gene = factor_rank_lists$final_rank[top_idx],
                                             Mut_class = factor_rank_lists$final_rank_class[top_idx],
                                             location_y = rep(10:6, each = 10),
                                             location_x = rep(1:10, 5)), 
                                  data.frame(Cor_num = signif_num_list$Bottom50, 
                                             Gene = factor_rank_lists$final_rank[bottom_idx],
                                             Mut_class = factor_rank_lists$final_rank_class[bottom_idx],
                                             location_y = rep(5:1, each = 10),
                                             location_x = rep(1:10, 5)))
signif_num_df$Cor_num <- factor(signif_num_df$Cor_num, levels = as.character(0:4))


p1 <- ggplot(signif_num_df, aes(y = location_y, x = location_x, label = Gene)) + 
  geom_point(aes(fill = Cor_num), color = "gray60", shape = 22, size = 7) + 
  geom_text(aes(color = Mut_class), vjust = -1.2) + 
  scale_fill_manual(values = list("0" = "white", "1" = "#BCC6DD", "2" = "#98A3CA", "3" = "#8092C4", 
                                  "4" = "#455D99", "CNV"= "#f87669", "SNV"= "#2fa1dd"),  
                    ) + 
  labs(y = "Genomic variations", x = "Rank split") +  
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1), legend.position = 'top')

pdf(file = "./Results/Comparison of corelation with CD8A_text.pdf", height = 8, width = 10)
p1
dev.off()