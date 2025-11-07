library(ggplot2)
library(IRanges)

setwd("~/proj/junks/50807_Gs_wild_rice/03_ssweep/s_region/perl/SOR")

my_PSR <- read.table("../PSR_2.txt", col.names = c("chr", "id", "start", "end"))
my_SOR <- my_PSR %>% 
  filter(id!=56)

myidx <- read.table("~/proj/Genome/O.sativa/IRGSP_build5/IRGSPb5.fa.masked.fai",sep = "\t")
myidx$chr <- 1:12

ggplot() +
  geom_segment(data=myidx, aes(x = 1, xend = V2, y = chr, yend = chr), color="grey", size = 3) +
  geom_segment(data = my_SOR, aes(x = start, xend = end, y = chr, yend = chr), color="red", size = 3) +
  scale_y_reverse() +
  theme_void()
