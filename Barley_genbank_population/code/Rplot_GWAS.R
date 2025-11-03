library(ggplot2)
library(dplyr)
library(scales)

setwd("F://projects/19_Barley_population/GWAS/Plot/")

gwas_YMV <- read.csv("../YMV/GAPIT.Association.GWAS_Results.MLM.V2.csv", header = T) |> 
  select(Chr, Pos, P.value) |> 
  mutate(type="YMV")

gwas_MMV <- read.csv("../MMV/GAPIT.Association.GWAS_Results.MLM.V2.csv", header = T) |> 
  select(Chr, Pos, P.value) |> 
  mutate(type="MMV")

gwas_bind <- rbind(gwas_YMV, gwas_MMV)

gwas_bind$Chr <- factor(gwas_bind$Chr)

threshold <- gwas_bind |> 
  group_by(type) |> 
  summarise(threshold=0.05/n())

ggplot() +
  geom_point(data = gwas_bind, aes(x=Pos, y=P.value, color = Chr)) +
  geom_hline(data = threshold, aes(yintercept = threshold), color="red") +
  facet_grid(type~Chr, scales = "free", space="free", switch = "both") +
  scale_y_continuous(trans = transform_compose("log10", "reverse"),
                     breaks = 10^(seq(0,-14,-2)), labels = function(x) -log10(x),
                     name = expression(-lg~italic(p))) +
  scale_x_continuous(breaks = seq(0,8e8,2e8), labels = function(x) x/1e6,
                     name = "Pos/Mb") +
  scale_color_manual(values=rep(c("black", "grey"), length.out=7)) +
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA, color = NA),
        strip.text.y = element_text(angle = 180),
        strip.placement = "outside")

