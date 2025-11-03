library(ggplot2)
library(dplyr)

setwd("~/proj/junks/41110_barley_population/02_plink_pca/Rplot/")

mydf <- read.table("../pca_10.eigenvec")
mydf <- mydf[,-1]
colnames(mydf) <- c("ENA_ID",paste0("PC",seq(1,10)))

myinfo <- read.table("180906_sample_information.tsv", header = T, sep = "\t", quote = "") %>% 
  filter(status_PCA=="domesticated") %>% 
  select(ENA_accession_id, country_of_origin) %>% 
  mutate(group=case_match(
    country_of_origin,
    c("FRA","DEU","GER","GBR","NLD","AUT","BEL","CHE","DNK","IRL","SWE","FIN","NOR",
      "EST","LVA","LTU","POL","CZE","SVK","HUN","BLR","RUS","SUN","DDR","CSK",
      "ALB","BGR","ROU","HRV","SRB","MKD","YUG", "UKR","MDA")~"NE",
    c("ESP","PRT","ITA", "GRC")~"WE",
    c("CHN","MNG","JPN","KOR","PRK")~"FE",
    c("AFG","PAK","IND","NPL","BTN","KAZ","UZB","KGZ","TJK","TKM")~"CA",
    c("TUR","GEO","ARM","AZE","IRN","IRQ","ISR","JOR","LBN","SYR", "CYP")~"ME",
    c("SAU","YEM","OMN","ARE")~"AP",
    c("EGY","LBY","TUN","DZA","MAR","SDN","TCD")~"NA",
    c("ERI","ETH")~"ETH",
    .default = "others"
  )) %>% 
  dplyr::rename(ENA_ID=ENA_accession_id)

while (F) {
  # 输出分组后的文件
  write.table(myinfo, "180906_sample_information_grouped.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
}

color_scale <- c("NE"="#e3181e", "WE"="#9728b2", "FE"="#3f27e3", "CA"="#060404", 
                 "ME"="#31ea47", "AP"="#f4df41", "NA"="#f38231", "ETH"="#21fcf6", "others"="grey")

mydf_infoed <- left_join(mydf,myinfo, by = "ENA_ID")

for (i in list(c("PC1","PC2"),c("PC3","PC4"),c("PC5","PC6"))) {
  p <- ggplot(data = mydf_infoed) +
    geom_point(aes(x=!!as.name(i[2]), y=!!as.name(i[1]), color=group), alpha=0.6) +
    scale_color_manual(values = color_scale)+
    labs(x=i[2], y=i[1], title = paste(i[2], "-", i[1])) +
#    coord_fixed(ratio = 1) +
    theme_classic()
  print(p)
}

# 绘制过滤后的材料

list_filter <- scan("../../04_admixture/Rplot/K12_filtered.list.txt", what = "character")

mydf_infoed_filter <- mydf_infoed %>% 
  mutate(filter=case_match(ENA_ID,
                           list_filter~"yes",
                           .default = "no"))

for (i in list(c("PC1","PC2"),c("PC3","PC4"),c("PC5","PC6"))) {
  p <- ggplot() +
    geom_point(data = mydf_infoed_filter %>% filter(filter=="no"),
               aes(x=!!as.name(i[2]), y=!!as.name(i[1])), color="grey") +
    geom_point(data = mydf_infoed_filter %>% filter(filter=="yes"),
               aes(x=!!as.name(i[2]), y=!!as.name(i[1]), color=group), alpha=0.6) +
    scale_color_manual(values = color_scale[1:8])+
    labs(x=i[2], y=i[1], title = paste(i[2], "-", i[1])) +
#    coord_fixed(ratio = 1) +
    theme_classic()
  print(p)
}

######
ggplot() +
  geom_point(data = mydf, aes(x=PC1, y=PC2, fill=PC3), size=1, shape=21, stroke=0) +
  geom_point(data = mydf_pg,aes(x=PC1, y=PC2, fill=PC3, color=name), shape=21, stroke=1) +
  geom_point(data = mydf[mydf$ENA_ID=="ERX2275088" | mydf$ENA_ID=="ERX2385673",],
             aes(x=PC1, y=PC2, fill=PC3), color="red", shape=21) +
  scale_fill_viridis_c() +
  scale_color_brewer(palette = "Dark2")
