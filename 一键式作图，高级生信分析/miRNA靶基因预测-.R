BiocManager::install("miRWalk")
rm(list = ls())
library(multiMiR)
library(edgeR)
library(openxlsx)
db.ver = multimir_dbInfoVersions()
db.ver

db.tables = multimir_dbTables()
db.tables

# Load data
setwd("D:/Desktop/")
data = read.xlsx("Diff_miRNA.xlsx",2)
multimir_results1 <- get_multimir(org     = 'hsa',
                                  mirna   = "hsa-miR-877-3p",
                                  predicted.site = "all",
                                  table   = "all",# search all validated tables "mirecords", "mirtarbase", and "tarbase"
                                  summary = FALSE)
write.xlsx(multimir_results1@data, file="hsa-miR-877-3p-靶基因预测.xlsx", asTable = FALSE, overwrite = TRUE)

# write.csv(multimir_results1@summary,file="miRNAdisease.drug summary.csv")
sessionInfo()