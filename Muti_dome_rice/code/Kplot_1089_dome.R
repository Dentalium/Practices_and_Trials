library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/proj/junks/50807_Gs_wild_rice/01_structure/admixture/Rplot/")

myinfo <- read.table("../../../00_rawdata/info/all_info_mapped.txt", sep = "\t", header = T)
#id_list <- myinfo %>% 
#  select(ID_vcf)

# 自定义组合因子水平用于排序样品
comb_levels <- c("indica","aus","rayada","aromatic","temperate japonica","tropical japonica" )

# 创建组合变量
#myinfo <- myinfo %>% 
#  mutate(comb_lv=paste(Group, Country, sep = "-"),
#         comb_lv=factor(comb_lv, levels = comb_levels))

# 提取出栽培稻样本并赋予因子水平
myinfo <- myinfo %>% 
  filter(Group %in% comb_levels) %>% 
  mutate(Group=factor(Group, levels = comb_levels))

#id_list <- myinfo %>% 
#  select(ID_vcf)

id_list <- read.table("../../../00_rawdata/1089_dome.txt", col.names = "ID_vcf")

# 样品水平
sample_lv <- myinfo %>% 
  select(ID_vcf, Group) %>% 
  arrange(Group) %>% 
  pull(ID_vcf)

myk <- data.frame()

for (i in 2:8) {
  qmatrix <- read.table(paste0("../1089_dome.pruned.",i,".Q"),
                        header = F, col.names = as.character(seq(1,i)), check.names = F) %>% 
    bind_cols(id_list) %>%
    mutate(K=i) %>% 
    pivot_longer(cols = 1:last_col(offset = 2), names_to = "K_group", values_to = "Q")
  
  myk <- myk %>% 
    bind_rows(left_join(myinfo, qmatrix, by = "ID_vcf"))
}

myk$ID_vcf <- factor(myk$ID_vcf, levels = sample_lv)

ggplot(data = myk) +
  geom_col(aes(x=ID_vcf, y=Q, group=K_group, fill=K_group), position = "stack") +
  scale_y_continuous(labels=NULL, breaks = NULL, name = NULL) +
  scale_x_discrete(labels=NULL, name=NULL) +
  facet_grid(K~., switch="y") +
  labs(subtitle = "K:") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, size = 2),
    legend.position = "none"
        )


