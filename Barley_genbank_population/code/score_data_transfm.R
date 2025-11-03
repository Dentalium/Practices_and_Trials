library(dplyr)
library(tidyr)

setwd("F://projects/19_Barley_population/rawdata/")

myinfo <- read.table("180906_sample_information.tsv", header = T, sep = "\t", fill=T, quote = "")

mydf_YMV <- read.table("YMV_raw.txt", sep = "\t", header = T) 

mydf_YMV_filter <- mydf_YMV |> select(accession_number,average_score) |> 
  group_by(accession_number) |> 
  summarise(score=mean(average_score)) |> 
  inner_join(myinfo, by = c("accession_number"="accession")) |> 
  select(ENA_accession_id, score)

write.table(mydf_YMV_filter, "YMV_average.txt", sep = "\t", quote = F, col.names = F, row.names = F)

mydf_MMV <- read.table("MMV_raw.txt", sep = "\t", header = T) 

mydf_MMV_filter <- mydf_MMV |> select(accession_number,average_score) |> 
  group_by(accession_number) |> 
  summarise(score=mean(average_score)) |> 
  inner_join(myinfo, by = c("accession_number"="accession")) |> 
  select(ENA_accession_id, score)

write.table(mydf_MMV_filter, "MMV_average.txt", sep = "\t", quote = F, col.names = F, row.names = F)

