library(ggplot2)
library(tidyr)
library(dplyr)

setwd("~/proj/junks/41110_barley_population/04_admixture/Rplot/")

# K值
Kvalue <- 7

df_reorder <- function(df, max_score, ancestors) {
  # 调整顺序，先排序最高得分祖先，然后依序号排序
  ancestors <- c(max_score, setdiff(ancestors, max_score))
  
  ancestor_args <- lapply(ancestors, function(i) df[[i]])
  
  myarg <- c(ancestor_args, decreasing = T)
  
#  print(myarg)
#  print(do.call(order, args = myarg))
  
  df <- df[do.call(order, args = myarg),]
  
  return(df)
}

# 样本名（行名）
samplenames <- scan("../samples.txt", what = character())

# Q矩阵
mydf <- read.table(paste0("../domesticated_barley_19778_samples_76102_SNPs.",Kvalue,".Q"))
mydf <- cbind(samplenames, mydf)

# 所有祖先名
ancestors <- paste0("A",1:Kvalue)

colnames(mydf) <- c("ID", ancestors)

# 读取分组文件
group_info <- read.table("../../02_plink_pca/Rplot/180906_sample_information_grouped.tsv",
                         header = T, sep = "\t", na.strings = "") %>% 
  rename(ID=ENA_ID)

# 将分组信息合并到主数据框
mydf <- merge(mydf, group_info, by = "ID", all.x = TRUE)

mydf$group <- factor(mydf$group, levels = c("others","NE","WE","AP","ME","CA","FE","ETH","NA"))

# 按照最高得分分群
max_score <- colnames(mydf)[2:(Kvalue+1)][apply(mydf[2:(Kvalue+1)], 1, which.max)]
mydf$max_score <- max_score

# 先按分组排序，然后在每个组内按祖先成分排序
mydf_split_by_group <- split(mydf, mydf$group)

mydf_reordered <- data.frame()
for (group_name in names(mydf_split_by_group)) {
  group_df <- mydf_split_by_group[[group_name]]
  
  # 在组内按照最高得分分群
  group_split <- split(group_df, group_df$max_score)
  
  group_reordered <- data.frame()
  for (i in names(group_split)) {
    # 排序
    group_split[[i]] <- df_reorder(group_split[[i]], i, ancestors)
    group_reordered <- rbind(group_reordered, group_split[[i]])
  }
  
  mydf_reordered <- rbind(mydf_reordered, group_reordered)
}

mydf_reordered <- pivot_longer(mydf_reordered, cols = starts_with("A"),
                               names_to = "type", values_to = "score")

# 设置因子水平保持排序
mydf_reordered$ID <- factor(mydf_reordered$ID, levels = unique(mydf_reordered$ID))

ggplot() +
  geom_col(data = mydf_reordered, aes(x=ID, y=score, fill=type), position = "stack", width = 1) +
  scale_fill_manual(values = c("A3"="#e3181e", "A2"="#060404", "A1"="#21fcf6",
                               "A6"="#3f27e3","A4"="#9728b2","A7"="#31ea47","A5"="yellow"),
                    name="Inffered\nancestor\nportion") +
  scale_y_continuous(name = NULL, breaks = NULL) +
  scale_x_discrete(name = NULL, breaks = NULL) +
#  scale_x_discrete(breaks = na.omit(pg_sample$V2)) +
  coord_cartesian(expand = F) +
  labs(title = paste("K =",Kvalue)) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())

while (F) {
  # 输出样本顺序
  write.table(data.frame(levels(mydf_reordered$ID)), "ordered_list.txt", col.names = F, row.names = F, quote = F)
}

