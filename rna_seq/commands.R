setwd('./GSE63310_RAW/')

files <- c("GSM1545535_10_6_5_11.txt","GSM1545536_9_6_5_11.txt",   "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",   "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",   "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt",   "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma", version = "3.8")

library(limma)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")

library(edgeR)

	
