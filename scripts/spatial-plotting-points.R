library(tidyverse)
library(sp)
library(rgdal)
library(raster)
library(sf)
library(ggthemes)
library(RColorBrewer)
source("scripts/spatial-averaging.R")

soilpts = st_read("data/soil-shp/", "wrbfu_pts")
RSG = readxl::read_xlsx("data/WRB_RSGs.xlsx", sheet = "RSG") %>%
  dplyr::select(RSG, Code, Class, Description, Class_06, Description_06)
sp_sp = readOGR("data/soil-shp", "wrbfu_pts")
soildf = create_soil_df(sp_sp, RSG)
crs = 3035
soilpts2 = st_transform(soilpts, crs)
crs2 = st_crs(soilpts2)

get_similarity_treshold = function(decay, min, sp) {
  #sim50 = decay %>% filter(sim >= 0.5)
  #comp50 = max(sim50$dist.km)
  
  test = sp %>% dplyr::select(similarity, dist)
  test$similarity = ifelse(test$dist <= min$dist.km[1], min$sim[1], 0)
  return(test)
}

plot_similarity_zones = function(p1) {
  set.seed(12345)
  p50 = sample.int(nrow(soilpts), nrow(soilpts)*0.5)
  p10 = sample.int(nrow(soilpts), nrow(soilpts)*0.1)
  p5 = sample.int(nrow(soilpts), nrow(soilpts)*0.05)
  point = soilpts2[p1, ]
  sp50 = soilpts2[p50, ]
  sp10 = soilpts2[p10, ]
  sp5 = soilpts2[p5, ]
  
  dd100 = get_decay(sp_sp[p1,], sp_sp, soildf[p1,], soildf)
  mins100 = dd100[findValleys(dd100$sim),]
  dd50 = get_decay(sp_sp[p1,], sp_sp[p50,], soildf[p1,], soildf[p50,])
  mins50 = dd50[findValleys(dd50$sim),]
  dd10 = get_decay(sp_sp[p1,], sp_sp[p10,], soildf[p1,], soildf[p10,])
  mins10 = dd10[findValleys(dd10$sim),]
  dd5 = get_decay(sp_sp[p1,], sp_sp[p5,], soildf[p1,], soildf[p5,])
  mins5 = dd5[findValleys(dd5$sim),]
  
  soilpts2$dist = (get_spatial_distance(sp_sp[p1,], sp_sp)$mat*111.32)
  sp50$dist = (get_spatial_distance(sp_sp[p1,], sp_sp[p50,])$mat*111.32)
  sp10$dist = (get_spatial_distance(sp_sp[p1,], sp_sp[p10,])$mat*111.32)
  sp5$dist = (get_spatial_distance(sp_sp[p1,], sp_sp[p5,])$mat*111.32)
  
  soilpts2$similarity = get_simple_similarity(soildf[p1,], soildf)$sim
  sp50$similarity = get_simple_similarity(soildf[p1,], soildf[p50,])$sim
  sp10$similarity = get_simple_similarity(soildf[p1,], soildf[p10,])$sim
  sp5$similarity = get_simple_similarity(soildf[p1,], soildf[p5,])$sim
  
  soilpts2$similarity2 = get_similarity_treshold(dd100, mins100, soilpts2)$similarity
  sp50$similarity2 = get_similarity_treshold(dd50, mins50, sp50)$similarity
  sp10$similarity2 = get_similarity_treshold(dd10, mins10, sp10)$similarity
  sp5$similarity2 = get_similarity_treshold(dd5, mins5, sp5)$similarity
  
  pp100 = ggplot(data = soilpts2) +
    geom_sf(aes(color = similarity2), size = 0.5) +
    geom_sf(data = point, fill = "orange", color = "orange", shape = 23, size = 1.5) +
    theme_void() +
    scale_color_gradient(low = "grey", high = "blue") +
    labs(title = "100% of points", color = "similarity") +
    theme(axis.text = element_blank(), legend.text = element_text(size = 4), title = element_text(size = 6))
  #plot(pp100)
  
  pp50 = ggplot(data = sp50) +
    geom_sf(aes(color = similarity2), size = 0.5) +
    geom_sf(data = point, fill = "orange", color = "orange", shape = 23, size = 1.5) +
    theme_void() +
    scale_color_gradient(low = "grey", high = "blue") +
    labs(title = "50% of points", color = "similarity") +
    theme(axis.text = element_blank(), legend.text = element_text(size = 4), title = element_text(size = 6))
  
  pp10 = ggplot(data = sp10) +
    geom_sf(aes(color = similarity2), size = 0.75) +
    geom_sf(data = point, fill = "orange", color = "orange", shape = 23, size = 1.5) +
    theme_void() +
    scale_color_gradient(low = "grey", high = "blue") +
    labs(title = "10% of points", color = "similarity") +
    theme(axis.text = element_blank(), legend.text = element_text(size = 4), title = element_text(size = 6))
  
  pp5 = ggplot(data = sp5) +
    geom_sf(aes(color = similarity2), size = 0.75) +
    geom_sf(data = point, fill = "orange", color = "orange", shape = 23, size = 1.5) +
    theme_void() +
    scale_color_gradient(low = "grey", high = "blue") +
    labs(title = "5% of points", color = "similarity") +
    theme(axis.text = element_blank(), legend.text = element_text(size = 4), title = element_text(size = 6))
  
  return(cowplot::plot_grid(pp100, pp50, pp10, pp5, nrow = 2))
}

p1 = 20031
#p1 = 20365
plotnum1 = plot_similarity_zones(p1)
#plot(plotnum1)

p2 = 12127 
#p2 = 12268
#p2 = 12297
plotnum2 = plot_similarity_zones(p2)
#plot(plotnum2)


#cowplot::plot_grid(plotnum2, plotnum1, nrow = 1, labels = c("a)", "b)"), scale = 0.95)

ggsave("output/figures/Fig4.png", cowplot::plot_grid(plotnum1, plotnum2, nrow = 1, labels = c("a)", "b)"), scale = 0.9), 
       dpi = 600, width = 9, height = 4.5)
