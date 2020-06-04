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
get_dist_decay_line = function(spatdist, catdist, type = "soil") {
  spat = NULL
  if(type == "soil") {
    ## read in spatial distance matrix
    spat = spatdist %>%
      rename(dist.deg = mat) %>%
      mutate(X = as.numeric(rownames(spatdist)),
             dist.m = dist.deg*111139, #approximate number of meters in a degree
             dist.km = dist.m/1000)
    #hist(spat$dist.km)
  }else {
    spat = spatdist %>%
      rename(dist.m = mat) %>%
      mutate(X = as.numeric(rownames(spatdist)),
             dist.km = dist.m/1000)
  }
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
get_decay = function(point, soilpts, pointdf, alldf, type = "soil") {
  #get spatial distances
  spatdist = get_spatial_distance(point, soilpts)
  
  catdist = NULL
  #get similarity matrix
  if(type == "soil") { catdist = get_simple_similarity(pointdf, alldf) }
  else { catdist = get_euclid_similarity(pointdf, alldf) }
  
  #then do something with distance decay function
  dd = get_dist_decay_line(spatdist, catdist)
  return(dd)
}

#get local minimum data for different proportions of points
get_all_mins = function(soilpts, soildf) {
  n = nrow(soilpts@data)
  point = sample.int(n, 1)
  soilpts90 = sample.int(n, (n*0.9))
  soilpts80 = sample.int(n, (n*0.8))
  soilpts70 = sample.int(n, (n*0.7))
  soilpts60 = sample.int(n, (n*0.6))
  soilpts50 = sample.int(n, (n*0.5))
  soilpts40 = sample.int(n, (n*0.4))
  soilpts30 = sample.int(n, (n*0.3))
  soilpts20 = sample.int(n, (n*0.2))
  soilpts10 = sample.int(n, (n*0.1))
  point.sp = soilpts[point,]
  point.df = soildf[point,]
  
  dd100 = get_decay(point.sp, soilpts, soildf[point,], soildf)
  dd90 = get_decay(point.sp, soilpts[soilpts90,], point.df, soildf[soilpts90,])
  dd80 = get_decay(point.sp, soilpts[soilpts80,], point.df, soildf[soilpts80,])
  dd70 = get_decay(point.sp, soilpts[soilpts70,], point.df, soildf[soilpts70,])
  dd60 = get_decay(point.sp, soilpts[soilpts60,], point.df, soildf[soilpts60,])
  dd50 = get_decay(point.sp, soilpts[soilpts50,], point.df, soildf[soilpts50,])
  dd40 = get_decay(point.sp, soilpts[soilpts40,], point.df, soildf[soilpts40,])
  dd30 = get_decay(point.sp, soilpts[soilpts30,], point.df, soildf[soilpts30,])
  dd20 = get_decay(point.sp, soilpts[soilpts20,], point.df, soildf[soilpts20,])
  dd10 = get_decay(point.sp, soilpts[soilpts10,], point.df, soildf[soilpts10,])
  
  mins100 = dd100[findValleys(dd100$sim),]
  mins90 = dd90[findValleys(dd90$sim),]
  mins80 = dd80[findValleys(dd80$sim),]
  mins70 = dd70[findValleys(dd70$sim),]
  mins60 = dd60[findValleys(dd60$sim),]
  mins50 = dd50[findValleys(dd50$sim),]
  mins40 = dd40[findValleys(dd40$sim),]
  mins30 = dd30[findValleys(dd30$sim),]
  mins20 = dd20[findValleys(dd20$sim),]
  mins10 = dd10[findValleys(dd10$sim),]
  
  allmins = data.frame(prop = c(10,20,30,40,50,60,70,80,90,100), 
                       dist.km = c(ifelse(nrow(mins10) == 0, NA, mins10$dist.km[1]), 
                                   ifelse(nrow(mins20) == 0, NA, mins20$dist.km[1]), 
                                   ifelse(nrow(mins30) == 0, NA, mins30$dist.km[1]), 
                                   ifelse(nrow(mins40) == 0, NA, mins40$dist.km[1]), 
                                   ifelse(nrow(mins50) == 0, NA, mins50$dist.km[1]), 
                                   ifelse(nrow(mins60) == 0, NA, mins60$dist.km[1]),
                                   ifelse(nrow(mins70) == 0, NA, mins70$dist.km[1]), 
                                   ifelse(nrow(mins80) == 0, NA, mins80$dist.km[1]), 
                                   ifelse(nrow(mins90) == 0, NA, mins90$dist.km[1]),
                                   ifelse(nrow(mins100) == 0, NA, mins100$dist.km[1])), 
                       sim = c(ifelse(nrow(mins10) == 0, NA, mins10$sim[1]), 
                               ifelse(nrow(mins20) == 0, NA, mins20$sim[1]), 
                               ifelse(nrow(mins30) == 0, NA, mins30$sim[1]), 
                               ifelse(nrow(mins40) == 0, NA, mins40$sim[1]), 
                               ifelse(nrow(mins50) == 0, NA, mins50$sim[1]), 
                               ifelse(nrow(mins60) == 0, NA, mins60$sim[1]),
                               ifelse(nrow(mins70) == 0, NA, mins70$sim[1]), 
                               ifelse(nrow(mins80) == 0, NA, mins80$sim[1]), 
                               ifelse(nrow(mins90) == 0, NA, mins90$sim[1]),
                               ifelse(nrow(mins100) == 0, NA, mins100$sim[1])))
  #allmins$rel.dif = allmins$dist.km - allmins$dist.km[10]
  
  return(allmins)
}

# #read in files
soilpts = readOGR("data/soil-shp/", "wrbfu_pts")
RSG = readxl::read_xlsx("data/WRB_RSGs.xlsx", sheet = "RSG") %>%
  dplyr::select(RSG, Code, Class, Description, Class_06, Description_06)
soildf = create_soil_df(soilpts, RSG)

#run through sampling multiple times
minimums = get_all_mins(soilpts, soildf)
for(i in 2:100) { #can change the second number to run a certain number of times
  minimums = rbind(minimums, get_all_mins(soilpts, soildf))
}

#write.csv(minimums, "output/results/minimums_lang2.csv")

