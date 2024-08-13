# Working directory
setwd("/~")

###
# 1. Expression Data preparation --------------------
#' 
#'###############################################################################
PanCancer.exp <- readRDS(file = "./Data/TCGA/PanCancer.exp.rds") # 17957 10176
colnames(PanCancer.exp) <- substr(colnames(PanCancer.exp), 1, 16)



###
# 2. CD8 T cell infiltration score calculation --------------------
#' 
#'###############################################################################
CD8_Tcell_signature_geneset <- readRDS(file = "./Results/CD8_Tcell_signature_geneset.rds")

library(GSVA)   
CD8T_ssGSEAScore <- gsva(expr = as.matrix(PanCancer.exp), gset.idx.list = CD8_Tcell_signature_geneset, 
                    method = "ssgsea", ssgsea.norm = FALSE) # 1 8223 matrix


CD8T_score_df <- data.frame(PatientID = colnames(CD8T_ssGSEAScore), 
                            CD8T_score = CD8T_ssGSEAScore[1,],
                            scale_CD8T_score = scale(CD8T_ssGSEAScore[1,]),
                            log2_CD8T_score = log2(CD8T_ssGSEAScore[1,])
)
save(CD8T_score_df, file = "./Data/TCGA/Solid_tumor_CD8T_score_df.RData")
write.csv(CD8T_score_df, file = "./Data/TCGA/Solid_tumor_CD8T_score_df.csv", quote = FALSE,
          row.names = FALSE)
