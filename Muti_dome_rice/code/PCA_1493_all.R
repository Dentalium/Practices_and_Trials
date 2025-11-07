library(ggplot2)
library(dplyr)

setwd("~/proj/junks/50807_Gs_wild_rice/01_structure/pca/rplot/")

myvec <- read.table("../1493_all_pca20.eigenvec")
myvec <- myvec[,2:5]
colnames(myvec) <- c("ID_vcf","PC1","PC2","PC3")

myinfo <- read.table("all_info_mapped.txt", header = T, sep = "\t") %>% 
  select(ID,ID_vcf,Group)

myvec_infoed <- myvec %>% 
  left_join(myinfo, by = "ID_vcf")

colorscale <- c("Niv1"="#dca432","Niv2"="#cd5144",
                "Ruf1a"="#33b033","Ruf1b"="#33b075","Ruf2"="#1467b2",
                "indica"="#c92559", "aus"="#8a2882", "rayada"="#89ba42", "aromatic"="#cbd23f",
                "temperate japonica"="#369fce", "tropical japonica"="#2a9e96")

ggplot(data = myvec_infoed) +
  geom_point(aes(x=PC1, y=PC2, color=Group), alpha=0.8) +
  scale_color_manual(values = colorscale) +
  scale_y_reverse() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = c(0.75,0.78),
    legend.title = element_blank(),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(
    ncol = 2,
    byrow = F
  ))


