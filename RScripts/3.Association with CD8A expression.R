###############################
#' TODO: Calculation of association of driver SNVs/CNVs with CD8A expression
#'      
###############################################################################
# Working directory
setwd("/~")

###
# 1.  Expression Data of TCGA solid tumor --------------------
#' 
#'###############################################################################
PanCancer.exp <- readRDS(file = "./Data/TCGA/PanCancer.exp.rds") # 17957 10176
colnames(PanCancer.exp) <- substr(colnames(PanCancer.exp), 1, 16)



###
# 2. mutation status matrix of quality-controlled solid tumor RNA-seq samples --------------------
#' 
#'###############################################################################
load(file = "./Data/TCGA/Solid_tumor_mutation_status_mat.RData") # 487 8223
rownames(mutation_status_mat) <- gsub("-|\\ ", "_", rownames(mutation_status_mat))



###
# 3. GISTIC 2.0 CNV matrix of quality-controlled solid tumor RNA-seq samples --------------------
#' 
#'###############################################################################
GISTIC_results <- read.csv("./Data/TCGA/all_thresholded.by_genes_whitelisted.tsv", 
                           sep = "\t", check.names = FALSE) # 25128 10716
GISTIC_results$'Gene Symbol' <- gsub("-|\\ ", "_", GISTIC_results$'Gene Symbol')
length(intersect(rownames(mutation_status_mat), GISTIC_results$'Gene Symbol')) # 487 个driver 交集 


# driver CNV matrix
driver_CNV_matrix <- GISTIC_results[which(GISTIC_results$'Gene Symbol' %in% rownames(mutation_status_mat)),]
colnames(driver_CNV_matrix)[4:10716] <- substr(colnames(driver_CNV_matrix)[4:10716], 1, 16)
rownames(driver_CNV_matrix) <- driver_CNV_matrix$'Gene Symbol'



###
# 4. common samples with matched SNV, CNV and expression data --------------------
#' 
#'###############################################################################
# identical(colnames(mutation_status_mat), colnames(driver_CNV_matrix))
# identical(colnames(driver_CNV_matrix), colnames(PanCancer.exp))

common_samples <- Reduce(intersect, list(colnames(mutation_status_mat), colnames(driver_CNV_matrix), 
                                         colnames(PanCancer.exp))) # 7973

mutation_status_mat <- mutation_status_mat[, common_samples]
driver_CNV_matrix <- driver_CNV_matrix[, common_samples]
PanCancer.exp <- PanCancer.exp[, common_samples]



###  
#' 5. Calculation of association of driver SNVs/CNVs with CD8A expression --------------------
#' Wilcoxon's rank-sum test 
#'###############################################################################

### (1) SNVs
load("./Results/SNV_wilcox_res.RData")
SNV_wilcox_res <- subset(SNV_wilcox_res, p.value < 0.05) 
candidate_SNV <- SNV_wilcox_res$SNV_gene # 204
  

data_merge <- cbind.data.frame(CD8A = as.numeric(PanCancer.exp["CD8A", ]), 
                               t(mutation_status_mat[candidate_SNV, ]))

# each SNV
library(coin)
wilcox_res <- lapply(candidate_SNV, function(g){
  tmp.data <- data_merge[, c("CD8A", g)]
  tmp.data[, g] <- factor(tmp.data[, g], levels = c("1", "0"))
  tmp.result <- wilcox_test(as.formula(paste0("CD8A ~ ", g)), data = tmp.data)
  refine.res <- data.frame(SNV_gene = g, Zvalue = statistic(tmp.result, type = "standardized"), p.value = pvalue(tmp.result))
  return(refine.res)
})

SNV_CD8A_wilcox_res <- do.call(rbind.data.frame, wilcox_res)
colnames(SNV_CD8A_wilcox_res)[2] <- "Zvalue"
SNV_CD8A_wilcox_res$p.adj <- p.adjust(SNV_CD8A_wilcox_res$p.value, method = "BH")

length(which(SNV_CD8A_wilcox_res$'p.value' < 0.05)) # 159
length(which(SNV_CD8A_wilcox_res$'p.adj' < 0.05)) # 155

save(SNV_CD8A_wilcox_res, file = "./Results/SNV_CD8A_wilcox_res.RData")

SNV_CD8A_wilcox_res <- subset(SNV_CD8A_wilcox_res, p.adj < 0.05)
write.csv(SNV_CD8A_wilcox_res, file = "./Results/SNV_CD8A_wilcox_res.csv", quote = FALSE, row.names = FALSE)


### (2) CNVs
load("./Results/CNV_wilcox_res.RData")
CNV_wilcox_res <- subset(CNV_wilcox_res, p.value < 0.05) 
candidate_CNV <- CNV_wilcox_res$CNV_gene # 487

# combine CNV class
driver_CNV_matrix_refine <- sapply(driver_CNV_matrix, function(x) {
  res <- gsub("1|2", "Amplified", (gsub("-1|-2", "Deleted", x)))
  return(res)}
)
rownames(driver_CNV_matrix_refine) <- rownames(driver_CNV_matrix) 


data_merge <- cbind.data.frame(CD8A = as.numeric(PanCancer.exp["CD8A", ]), 
                               t(driver_CNV_matrix_refine))

# each CNV
library(coin)
wilcox_res <- lapply(rownames(driver_CNV_matrix_refine), function(g){
  tmp.data <- data_merge[, c("CD8A", g)]
  tmp <- sort(table(data_merge[, g]))
  CNV_class <- names(tmp)[2:3]
  if(!is.element("0", CNV_class)) {CNV_class[1] <- "0"} 
  
  tmp.data <- subset(tmp.data, tmp.data[, g] %in% CNV_class)
  tmp.data[, g] <- factor(tmp.data[, g], levels = c(setdiff(CNV_class, "0"), "0"))
  
  tmp.result <- wilcox_test(as.formula(paste0("CD8A ~ ", g)), data = tmp.data)
  refine.res <- data.frame(CNV_gene = g, Zvalue = statistic(tmp.result, type = "standardized"), p.value = pvalue(tmp.result))
  
  return(refine.res)
})

CNV_CD8A_wilcox_res <- do.call(rbind.data.frame, wilcox_res)
colnames(CNV_CD8A_wilcox_res)[2] <- "Zvalue"

CNV_CD8A_wilcox_res$p.adj <- p.adjust(CNV_CD8A_wilcox_res$p.value, method = "BH")

length(which(CNV_CD8A_wilcox_res$'p.value' < 0.05)) # 271
length(which(CNV_CD8A_wilcox_res$'p.adj' < 0.05)) # 250
save(CNV_CD8A_wilcox_res, file = "./Results/CNV_CD8A_wilcox_res.RData")

CNV_CD8A_wilcox_res <- subset(CNV_CD8A_wilcox_res, p.adj < 0.05)
write.csv(CNV_CD8A_wilcox_res, file = "./Results/CNV_CD8A_wilcox_res.csv", quote = FALSE, row.names = FALSE)

