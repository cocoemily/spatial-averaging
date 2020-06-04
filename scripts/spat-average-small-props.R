#full spatial averaging script for HPC
#load libraries
library(tidyverse)
library(sp)
library(rgdal)
library(raster)
library(readxl)
library(mgcv)
library(quantmod)

#spatial distances function
get_spatial_distance <- function(point, points) {
  mat = pointDistance(point, points, allpairs = T, lonlat = F)
  return(as.data.frame(mat))
}

#soil dataframe function
create_soil_df <- function(shp, soil) {
  df <- shp@data %>%
    dplyr::select(WRBFU, SLOPEDO, TXSUBDO, DR, ATC, PD_SUB, STR_SUB, MIN_SUB, 
                  TXSRFDO, PD_TOP, STR_TOP, MIN_TOP, PMH)
  df$RSG <- substr(as.character(df$WRBFU), 1, 2)
  df$LL <- substr(as.character(df$WRBFU), 3, length(df$WRBFU))
  df$Class <- 0
  for(i in 1:nrow(df)) {
    df$Class[i] <- RSG$Class_06[RSG$Code == df$RSG[i]]
  }
  return(df)
}

#similarity distance function
get_simple_similarity <- function(med.pt, other.pts) {
  df <- data.frame(id = 1:nrow(other.pts), 
                   sim = NA)
  med = med.pt[1,]
  for(i in 1:nrow(other.pts)) {
    op = other.pts[i,]
    match = op == med
    match = match[!is.na(match)]
    df$sim[i] = sum(match)/length(match)
  }
  return(df)
}

#distance decay dataframe function
get_dist_decay_line = function(spatdist, catdist) {
  ## read in spatial distance matrix
  spat = spatdist %>%
    rename(dist.deg = mat) %>%
    mutate(X = as.numeric(rownames(spatdist)),
           dist.m = dist.deg*111139, #approximate number of meters in a degree
           dist.km = dist.m/1000)
  #hist(spat$dist.km)
  ## read in categorical similarity matrix
  cat <- catdist %>%
    rename(X = id)
  #hist(cat$sim)
  
  dist.dec <- cat %>%
    left_join(spat, by = "X")
  
  #fit smoothing curve
  fit1 <- gam(sim ~ s(dist.km, k=10), data = dist.dec)
  # par(mfrow = c(2,2))
  # gam.check(fit1)
  # summary(fit1)
  # par(mfrow = c(1,1))
  # plot(fit1)
  
  #with(dist.dec, cor(dist.km, sim))
  
  ppdf <- data.frame(dist.km = seq(min(range(dist.dec$dist.km)), max(range(dist.dec$dist.km))))
  ppdf$sim <- predict.gam(fit1, newdata = ppdf)
  
  return(ppdf)
}

#function to put workflow together
get_decay = function(point, soilpts, pointdf, alldf) {
  #get spatial distances
  spatdist = get_spatial_distance(point, soilpts)
  
  #get similarity matrix
  catdist = get_simple_similarity(pointdf, alldf)
  
  #then do something with distance decay function
  dd = get_dist_decay_line(spatdist, catdist)
  return(dd)
}

#get local minimum data for different proportions of points
get_all_small_mins = function(soilpts, soildf) {
  n = nrow(soilpts@data)
  point = sample.int(n, 1)
  soilpts9 = sample.int(n, (n*0.09))
  soilpts8 = sample.int(n, (n*0.08))
  soilpts7 = sample.int(n, (n*0.07))
  soilpts6 = sample.int(n, (n*0.06))
  soilpts5 = sample.int(n, (n*0.05))
  soilpts4 = sample.int(n, (n*0.04))
  soilpts3 = sample.int(n, (n*0.03))
  soilpts2 = sample.int(n, (n*0.02))
  soilpts1 = sample.int(n, (n*0.01))
  point.sp = soilpts[point,]
  point.df = soildf[point,]
  
  dd100 = get_decay(point.sp, soilpts, soildf[point,], soildf)
  dd9 = get_decay(point.sp, soilpts[soilpts9,], point.df, soildf[soilpts9,])
  dd8 = get_decay(point.sp, soilpts[soilpts8,], point.df, soildf[soilpts8,])
  dd7 = get_decay(point.sp, soilpts[soilpts7,], point.df, soildf[soilpts7,])
  dd6 = get_decay(point.sp, soilpts[soilpts6,], point.df, soildf[soilpts6,])
  dd5 = get_decay(point.sp, soilpts[soilpts5,], point.df, soildf[soilpts5,])
  dd4 = get_decay(point.sp, soilpts[soilpts4,], point.df, soildf[soilpts4,])
  dd3 = get_decay(point.sp, soilpts[soilpts3,], point.df, soildf[soilpts3,])
  dd2 = get_decay(point.sp, soilpts[soilpts2,], point.df, soildf[soilpts2,])
  #dd1 = get_decay(point.sp, soilpts[soilpts1,], point.df, soildf[soilpts1,])
  
  mins100 = dd100[findValleys(dd100$sim),]
  mins9 = dd9[findValleys(dd9$sim),]
  mins8 = dd8[findValleys(dd8$sim),]
  mins7 = dd7[findValleys(dd7$sim),]
  mins6 = dd6[findValleys(dd6$sim),]
  mins5 = dd5[findValleys(dd5$sim),]
  mins4 = dd4[findValleys(dd4$sim),]
  mins3 = dd3[findValleys(dd3$sim),]
  mins2 = dd2[findValleys(dd2$sim),]
  #mins1 = dd1[findValleys(dd1$sim),]
  
  allmins1 = data.frame(#prop = c(1,2,3,4,5,6,7,8,9,100),
                        prop = c(2,3,4,5,6,7,8,9,100),
                       dist.km = c(#ifelse(nrow(mins1) == 0, NA, mins1$dist.km[1]), 
                                   ifelse(nrow(mins2) == 0, NA, mins2$dist.km[1]), 
                                   ifelse(nrow(mins3) == 0, NA, mins3$dist.km[1]), 
                                   ifelse(nrow(mins4) == 0, NA, mins4$dist.km[1]), 
                                   ifelse(nrow(mins5) == 0, NA, mins5$dist.km[1]), 
                                   ifelse(nrow(mins6) == 0, NA, mins6$dist.km[1]),
                                   ifelse(nrow(mins7) == 0, NA, mins7$dist.km[1]), 
                                   ifelse(nrow(mins8) == 0, NA, mins8$dist.km[1]), 
                                   ifelse(nrow(mins9) == 0, NA, mins9$dist.km[1]),
                                   ifelse(nrow(mins100) == 0, NA, mins100$dist.km[1])), 
                       sim = c(#ifelse(nrow(mins1) == 0, NA, mins1$sim[1]), 
                               ifelse(nrow(mins2) == 0, NA, mins2$sim[1]), 
                               ifelse(nrow(mins3) == 0, NA, mins3$sim[1]), 
                               ifelse(nrow(mins4) == 0, NA, mins4$sim[1]), 
                               ifelse(nrow(mins5) == 0, NA, mins5$sim[1]), 
                               ifelse(nrow(mins6) == 0, NA, mins6$sim[1]),
                               ifelse(nrow(mins7) == 0, NA, mins7$sim[1]), 
                               ifelse(nrow(mins8) == 0, NA, mins8$sim[1]), 
                               ifelse(nrow(mins9) == 0, NA, mins9$sim[1]),
                               ifelse(nrow(mins100) == 0, NA, mins100$sim[1])))
  #allmins1$rel.dif = allmins1$dist.km - allmins1$dist.km[10]
  
  return(allmins1)
}

#read in files
soilpts = readOGR("data/soil-shp/", "wrbfu_pts")
RSG = readxl::read_xlsx("data/WRB_RSGs.xlsx", sheet = "RSG") %>%
  dplyr::select(RSG, Code, Class, Description, Class_06, Description_06)
soildf = create_soil_df(soilpts, RSG)

#run through sampling multiple times
minimums = get_all_small_mins(soilpts, soildf)
for(i in 2:100) { #can change the second number to run a certain number of times
  minimums = rbind(minimums, get_all_small_mins(langpts, langdf))
}

#write.csv(minimums, "output/results/small_minimums.csv")

