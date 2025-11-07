library(dplyr)

setwd("F://projects/42_Gs_wildrice/00_rawdata/")
my_id <- read.table("wild_id_map.txt", header = T, sep = "\t") |> 
  filter(remove=="") |> 
  select(!remove)
my_info <- read.table("wild_info.txt", header = T, sep = "\t") |> 
  filter(Note=="") |> 
  select(!Note)

my_info_mapped <- my_info |> 
  left_join(my_id, by = c("ID", "species")) |> 
  tidyr::drop_na()

write.table(my_info_mapped, "wild_info_mapped.txt", quote = F, row.names = F, sep = "\t")

