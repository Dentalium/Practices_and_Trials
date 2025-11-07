library(ggplot2)
library(sf)
library(maps)

setwd("~/proj/junks/50807_Gs_wild_rice/00_rawdata/map/")

myinfo <- read.table("../info/wild_info_mapped.txt", sep = "\t", header = T)

colorscale <- c("Niv1"="#dca432","Niv2"="#cd5144",
                "Ruf1a"="#33b033","Ruf1b"="#33b075","Ruf2"="#1467b2")

world_map <- st_as_sf(map("world", plot = F, fill = T))

ggplot() +
  geom_sf(data = world_map, fill="white") +
  geom_point(data = myinfo, aes(y=Latitude, x=Longitude, colour = Group)) +
  scale_color_manual(values = colorscale, name=NULL) +
  coord_sf(xlim = c(72,145), ylim = c(-17,32)) +
  theme(
    panel.border = element_rect(color = "black", fill=NA),
    axis.title = element_blank(),
    legend.position = c(0.93,0.85),
    legend.key = element_rect(fill = NA)
  )


