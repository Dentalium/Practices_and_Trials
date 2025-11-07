library(dplyr)
library(tidyr)
library(ggplot2)
library(IRanges)
library(ggbreak)

# 根据Z值计算非标准正态变量对应值
Z_cutoff <- function(data, Z_value){
  return(Z_value*sd(data)+mean(data))
}

setwd("~/proj/junks/50807_Gs_wild_rice/03_ssweep/s_region/perl/")

mydf <- read.table("1355_wild+dome_maj.calculated.txt", header = T, sep = "\t")

mydf_filter0 <- mydf %>% select(chr:snp_num, dome.wild.fst, dome.wild.rod) %>% 
  filter(snp_num>=10)

mydf_filter0$chr <- factor(as.character(mydf_filter0$chr),
                          levels = as.character(seq(1,12)))

while (F) {
  ggplot(data = mydf_filter0, aes(sample=dome.wild.fst),
         dparams=list(mean(dome.wild.fst), sd(dome.wild.fst))) +
    geom_qq() +
    geom_qq_line() +
    labs(x=NULL,y=NULL, title = expression(F[ST])) +
    theme_classic()
  
  ggplot(data = mydf_filter0, aes(sample=dome.wild.rod),
         dparams=list(mean(dome.wild.rod), sd(dome.wild.rod))) +
    geom_qq() +
    geom_qq_line() +
    labs(x=NULL,y=NULL, title = expression(ROD)) +
    theme_classic()
}

# 绘图
while (F) {
  ggplot(data = mydf_filter0) +
    geom_col(aes(x=start_pos, y=dome.wild.fst, fill=chr)) +
    geom_hline(yintercept = 0.3123, color="grey", linetype="dashed") +
    scale_fill_manual(values = rep(c("black","darkgrey"),length.out=12)) +
    scale_y_continuous(expand = expansion(), name=expression("F"["ST"])) +
    facet_grid(.~chr, scales="free_x", space="free_x") +
    theme_classic() +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          strip.text = element_blank(), panel.spacing = unit(0,units = "cm"),
          legend.position = "none")
  
  ggplot(data = mydf_filter0) +
    geom_col(aes(x=start_pos, y=dome.wild.rod, fill=chr)) +
    geom_hline(yintercept = 2.753, color="grey", linetype="dashed") +
    scale_fill_manual(values = rep(c("black","darkgrey"),length.out=12)) +
#    scale_y_continuous(expand = expansion(), name = "ROD") +
    scale_y_break(breaks = c(20,60), expand = expansion()) +
    labs(y="ROD") +
    facet_grid(.~chr, scales="free_x", space="free_x") +
    theme_classic() +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          strip.text = element_blank(), panel.spacing = unit(0,units = "cm"),
          legend.position = "none")
    # 最大值约70，需要节选
#    coord_cartesian(ylim = c(0,25))
}


# ROD偏移过大的位点不参与计算
mydf_filter0 <- mydf_filter0 %>%
  mutate(type=case_when(
    dome.wild.rod<=mean(dome.wild.rod)+2*sd(dome.wild.rod) &
           dome.wild.rod>=mean(dome.wild.rod)-2*sd(dome.wild.rod)~"remain",
    .default = "omitted"
    ))

# z变换
Z_thred <- 2

print(paste("FST cutoff:", Z_cutoff(data = mydf_filter$dome.wild.fst, Z_value = Z_thred)))
print(paste("ROD cutoff:", Z_cutoff(data = mydf_filter$dome.wild.rod, Z_value = Z_thred)))

mydf_filter <- mydf_filter0 %>% 
  filter(type=="remain") %>% 
  mutate(dome.wild.fst_z=(dome.wild.fst-mean(dome.wild.fst))/sd(dome.wild.fst),
                                      dome.wild.rod_z=(dome.wild.rod-mean(dome.wild.rod))/sd(dome.wild.rod),) %>% 
  filter(dome.wild.fst_z>Z_thred & dome.wild.rod_z>Z_thred) %>% 
  select(!ends_with("_z")) %>% 
  bind_rows(mydf_filter0 %>% filter(type=="omitted"))

# 比较与正态分布的差异
while (F) {
  mean_dome.wild.fst <- mydf_filter0 %>%
    filter(type=="remain") %>% 
    summarise(mean=mean(dome.wild.fst)) %>% 
    pull(mean)
  sd_dome.wild.fst<- mydf_filter0 %>%
    filter(type=="remain") %>% 
    summarise(sd=sd(dome.wild.fst)) %>% 
    pull(sd)
  
  ggplot() +
    geom_histogram(data = mydf_filter0 %>% filter(type=="remain"),
                   aes(x=dome.wild.fst), bins = 200) +
    geom_freqpoly(data = data.frame(mock=rnorm(n = 37342,
                                              mean = mean_dome.wild.fst, sd = sd_dome.wild.fst)),
                  aes(x=mock),
                  bins = 200) +
    labs(x=NULL, title = expression(F[ST])) +
    coord_cartesian(expand = F) +
    theme_classic()
  
  mean_dome.wild.rod <- mydf_filter0 %>%
    filter(type=="remain") %>% 
    summarise(mean=mean(dome.wild.rod)) %>% 
    pull(mean)
  sd_dome.wild.rod<- mydf_filter0 %>%
    filter(type=="remain") %>% 
    summarise(sd=sd(dome.wild.rod)) %>% 
    pull(sd)
  
  ggplot() +
    geom_histogram(data = mydf_filter0 %>% filter(type=="remain"),
                   aes(x=dome.wild.rod), bins = 200) +
    geom_freqpoly(data = data.frame(mock=rnorm(n = 37342,
                                               mean = mean_dome.wild.rod, sd = sd_dome.wild.rod)),
                  aes(x=mock),
                  bins = 200) +
    labs(x=NULL, title = expression(ROD)) +
    coord_cartesian(expand = F) +
    theme_classic()
}

# 手动合并
while (F) {
  # 合并相邻PSR
  PSR_combined <- data.frame(chr=character(),
                             no=integer(),
                             start=integer(),
                             end=integer())
  
  current_chr <- 0
  current_start <- 1
  current_end <- 0
  counter=0
  
  for (i in 1:dim(mydf_filter)[1]) {
    this_chr <- mydf_filter[i,1]
    this_start <- mydf_filter[i,2]
    this_end <- mydf_filter[i,3]
    
    if (this_chr != current_chr) {
      PSR_combined <- rbind(PSR_combined, data.frame(
        chr=current_chr, no=counter, start=current_start, end=current_end
      ))
      
      current_chr <- this_chr
      current_end <- 0
      current_start <- this_start
      counter <- counter+1
    }
    else if (this_start > current_end) {
      PSR_combined <- rbind(PSR_combined, data.frame(
        chr=current_chr, no=counter, start=current_start, end=current_end
      ))
      
      current_end <- this_end
      current_start <- this_start
      counter <- counter+1
    }
    else {
      current_end <- this_end
    }
  }
}

# 使用IRanges合并
PSR_combined <- data.frame(chr=character(),
                           no=integer(),
                           start=integer(),
                           end=integer())
counter <- 0

for (i in as.character(1:12)) {
  this_chr <- filter(mydf_filter, chr==i)
  if (nrow(this_chr)!=0) {
    this_chr_PSR <- IRanges(start = this_chr$start_pos, end = this_chr$end_pos)
    this_chr_PSR <- reduce(this_chr_PSR)
    
    num_PSRs <- length(this_chr_PSR)
    counter <- counter+num_PSRs
    
    PSR_combined <- rbind(PSR_combined,
                          data.frame(chr=rep(i, num_PSRs), no=seq(counter-num_PSRs+1, counter), 
                                     start=IRanges::start(this_chr_PSR), end=IRanges::end(this_chr_PSR)))
  }
}

# 输出
while (F) {
  write.table(PSR_combined, "PSR_2.txt", row.names = F, quote = F, col.names = F, sep = "\t")
}

# 绘制联合图
# 提取数据
mydf_jtplt <- mydf %>% 
  filter(snp_num>=10) %>% 
  mutate(pos=start_pos+50000-1) %>% 
  select(chr, pos, dome.pi.1kb., wild.pi.1kb., dome.wild.fst, dome.wild.rod) %>% 
  # 初始化PSR状态
  mutate(PSR=F)

mydf_jtplt$chr <- factor(mydf_jtplt$chr, levels = as.character(seq(1,12)))  

# 选取PSR
for(i in as.character(seq(1,12))) {
  this_chr_PSR <- PSR_combined %>% 
    filter(chr == i)
  this_chr_window <- mydf_jtplt %>%
    filter(chr == i)
  
  # 创建IRanges对象
  PSR_ir <- IRanges(start = this_chr_PSR$start, end = this_chr_PSR$end)
  window_ir <- IRanges(start = this_chr_window$pos, end = this_chr_window$pos)
  
  # 查找重叠
  overlaps <- findOverlaps(PSR_ir, window_ir)
  
  # 写入PSR信息
  overlap_indices <- S4Vectors::subjectHits(overlaps)
  
  mydf_jtplt$PSR[mydf_jtplt$chr==i][overlap_indices] <- T
}

# 按变量类型分组
# type1：统计量类型，用于分面
# type2：颜色，包括以下水平：普通、PSR、两种材料的pi
mydf_jtplt <- mydf_jtplt %>% 
  pivot_longer(cols = dome.pi.1kb.:dome.wild.rod, values_to = "value", names_to = "type0") %>% 
  mutate(type1=case_match(
    type0,
    c("dome.pi.1kb.","wild.pi.1kb.")~"pi",
    "dome.wild.fst"~"fst",
    "dome.wild.rod"~"rod"
  )) %>% 
  mutate(type1=factor(type1, levels = c("pi","rod","fst"))) %>% 
  mutate(type2=case_when(
    PSR & type1 %in% c("fst","rod")~"PSR",
    (!PSR) & type1 %in% c("fst","rod")~"common",
    type0=="dome.pi.1kb."~"dome",
    type0=="wild.pi.1kb."~"wild"
  ))

# 分染色体绘图
for (i in as.character(1:12)) {
  p <- ggplot(data = mydf_jtplt %>% filter(chr==i)) +
    geom_point(aes(x=pos, y=value, color=type2), size=0.5) +
    facet_grid(type1~chr, space="fixed", scales="free", switch="both") +
    scale_color_manual(values = c("PSR"="#f50000", "common"="black", "dome"="#327e06", "wild"="#fc9a1c")) +
    scale_x_continuous(breaks = seq(0,40e6, 10e6), labels = function(x) paste0(as.character(x/1e6),"M")) +
    scale_y_continuous(limits = c(0, NA)) +
    theme_classic() +
    labs(x=NULL, y=NULL) +
    theme(
      legend.position = "none",
      strip.placement = "outside",
      strip.background = element_rect(color = NA)
    )
  
  ggsave(filename = paste0("plot/chr", i, ".joint.pdf"), plot = p,
         width = 12, height = 3, device = "pdf")
}

  

