library(GAPIT)

setwd("F://projects/19_Barley_population/GWAS/MMV/")

myG <- read.table("MMV.hmp.txt", header = F, sep = "\t", comment.char = "")
myY <- read.table("../../rawdata/MMV_average.txt", header = F, sep = "\t")

myGAPIT_MLM <- GAPIT(
  Y=myY,
  G=myG,
  model="MLM",
  PCA.total=6,
  SNP.MAF = 0.01,
  file.output=T
)
