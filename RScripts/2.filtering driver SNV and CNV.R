###############################
#' TODO: preprocess of filtering driver SNVs and CNVs in TCGA solid tumor 
#' Requirement: Correlated with CD8 + T cell infiltration score
#'  
###############################################################################
# Working directory
setwd("/~")

###################################################################

###
# 1. mutation status matrix of quality-controlled solid tumor RNA-seq samples --------------------
#' 8223 common samples with RNA-seq
#'###############################################################################

load(file = "./Data/TCGA/Solid_tumor_mutation_status_mat.RData") # 487 8223
rownames(mutation_status_mat) <- gsub("-|\\ ", "_", rownames(mutation_status_mat))


###
# 2. CD8 T cell infiltration score --------------------
#' 
#'###############################################################################
load(file = "./Data/TCGA/Solid_tumor_CD8T_score_df.RData")



###  
#' 3. Wilcoxon's rank-sum test of CD8+ T cell infiltration score between the samples with genomic variations or not
#' adjustp-value (method = "BH")
#' 
################################
library(coin)

identical(colnames(mutation_status_mat), CD8T_score_df$PatientID)
data.merge <- cbind.data.frame(CD8T_score = CD8T_score_df$scale_CD8T_score, 
                               t(mutation_status_mat))

# Each SNV
SNV_wilcox_res <- lapply(rownames(mutation_status_mat), function(g){
  tmp.data <- data.merge[, c("CD8T_score", g)]
  tmp.data[, g] <- factor(tmp.data[, g], levels = c("1", "0"))
  tmp.result <- wilcox_test(as.formula(paste0("CD8T_score ~ ", g)), data = tmp.data)
  refine.res <- data.frame(SNV_gene = g, Zvalue = statistic(tmp.result, type = "standardized"), p.value = pvalue(tmp.result))
  return(refine.res)
})


SNV_wilcox_res <- do.call(rbind.data.frame, wilcox_res)
colnames(SNV_wilcox_res)[2] <- "Zvalue"


save(SNV_wilcox_res, file = "./Results/SNV_wilcox_res.RData")


length(which(SNV_wilcox_res$'p.value' < 0.05)) # 204

SNV_wilcox_res <- subset(SNV_wilcox_res, p.value < 0.05)
write.csv(SNV_wilcox_res, file = "./Results/SNV_wilcox_res.csv", quote = FALSE, row.names = FALSE)



###
# 4. GISTIC 2.0 CNV matrix of quality-controlled solid tumor RNA-seq samples --------------------
#' 
#'###############################################################################
GISTIC_results <- read.csv("./Data/TCGA/all_thresholded.by_genes_whitelisted.tsv", 
                           sep = "\t", check.names = FALSE) # 25128 10716

length(intersect(rownames(mutation_status_mat), GISTIC_results$'Gene Symbol')) # 487 个driver gene

# driver CNV matrix
GISTIC_results <- GISTIC_results[which(GISTIC_results$'Gene Symbol' %in% rownames(mutation_status_mat)),]
colnames(GISTIC_results)[4:10716] <- substr(colnames(GISTIC_results)[4:10716], 1, 16)

common_samples <- intersect(colnames(mutation_status_mat), colnames(GISTIC_results)) # 7973
rownames(GISTIC_results) <- gsub("-|\\ ", "_", GISTIC_results$'Gene Symbol')



###
# 5. CNV matrix of quality-controlled solid tumor RNA-seq samples --------------------
#' 
#'###############################################################################
a1 <- sapply(GISTIC_results[, -(1:3)], function(x) {
  res <- gsub("1|2", "Amplified", (gsub("-1|-2", "Deleted", x)))
  return(res)}
)
rownames(a1) <- rownames(GISTIC_results) 
GISTIC_results_refine <- a1[, common_samples]


load("./Data/TCGA/Solid_tumor_CD8T_score_df.RData")
# identical(colnames(GISTIC_results_refine), CD8T_score_df$PatientID)
common_samples <- intersect(CD8T_score_df$PatientID, colnames(GISTIC_results_refine)) # 7973


data.merge <- cbind.data.frame(CD8T_score = CD8T_score_df$scale_CD8T_score[match(common_samples, CD8T_score_df$PatientID)], 
                               t(GISTIC_results_refine[, match(common_samples, colnames(GISTIC_results_refine))]))


# Each CNV
library(coin)

wilcox_res <- lapply(rownames(GISTIC_results_refine), function(g){
  # print(paste0(g, "----------"))
  tmp.data <- data.merge[, c("CD8T_score", g)]
  tmp <- sort(table(data.merge[, g]))
  CNV_class <- names(tmp)[2:3]
  if(!is.element("0", CNV_class)) {CNV_class[1] <- "0"} 
  
  tmp.data <- subset(tmp.data, tmp.data[, g] %in% CNV_class)
  tmp.data[, g] <- factor(tmp.data[, g], levels = c(setdiff(CNV_class, "0"), "0"))

  tmp.result <- wilcox_test(as.formula(paste0("CD8T_score ~ ", g)), data = tmp.data)
  refine.res <- data.frame(CNV_gene = g, Zvalue = statistic(tmp.result, type = "standardized"), p.value = pvalue(tmp.result))

  return(refine.res)
})


CNV_wilcox_res <- do.call(rbind.data.frame, wilcox_res)
colnames(CNV_wilcox_res)[2] <- "Zvalue"
save(CNV_wilcox_res, file = "./Results/CNV_wilcox_res.RData")

length(which(CNV_wilcox_res$'p.value' < 0.05)) # 487
