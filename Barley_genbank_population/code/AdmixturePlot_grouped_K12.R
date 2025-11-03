library(ggplot2)
library(tidyr)
library(dplyr)

setwd("~/proj/junks/41110_barley_population/04_admixture/Rplot/")

# K值
Kvalue <- 12

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

mydf_reordered <- pivot_longer(mydf, cols = starts_with("A"),
                               names_to = "type", values_to = "score")
# 样本顺序取自K7
sample_order <- scan("ordered_list.txt", what="character")
# 设置因子水平保持排序
mydf_reordered$ID <- factor(mydf_reordered$ID, levels = sample_order)

ggplot() +
  geom_col(data = mydf_reordered, aes(x=ID, y=score, fill=type), position = "stack", width = 1) +
  scale_fill_manual(values = c("A10"="#ff4444", "A11"="#e3181e", "A12"="red4",
                               "A5"="#9728b2","A4"="#060404", "A1"="#21fcf6",
                               "A3"="#3f27e3","A8"="#31ea47","A2"="yellow",
                               "A6"="darkgreen","A7"="#bb0044","A9"="lightblue"),
                    name="Inffered\nancestor\nportion") +
  scale_y_continuous(name = NULL, breaks = NULL) +
  scale_x_discrete(name = NULL, breaks = NULL) +
#  scale_x_discrete(breaks = na.omit(pg_sample$V2)) +
  coord_cartesian(expand = F) +
  labs(title = paste("K =",Kvalue)) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())

while (F) {
  # 输出最高得分大于0.7的样本列表
  # others组不输出
  mydf_fiter <- mydf %>% 
    left_join(group_info, by = "ID") %>% 
    filter(group!="others") %>% 
    rowwise() %>% 
    mutate(max_score_number=max(c_across(starts_with("A")))) %>% 
    filter(max_score_number>=0.7) %>% 
    select(ID)
  
  write.table(mydf_fiter, "K12_filtered.list.txt", col.names = F, row.names = F, quote = F)
}
