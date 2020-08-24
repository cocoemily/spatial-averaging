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
create_lang_df = function(shp) {
  df = shp@data %>%
    dplyr::select(category, family_pk, father_pk)
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
get_mins = function(soilpts, soildf, age) {
  numsites = num_sites(age, soilpts)
  n = nrow(soilpts@data)
  point = sample.int(n, 1)
  pts = sample.int(n, numsites)
  point.sp = soilpts[point,]
  point.df = soildf[point,]
  
  dd100 = get_decay(point.sp, soilpts, point.df, soildf)
  dd_prop = get_decay(point.sp, soilpts[pts,], point.df, soildf[pts,])
  
  mins100 = dd100[findValleys(dd100$sim),]
  mins_prop = dd_prop[findValleys(dd_prop$sim),]
  
  allmins = data.frame(prop = c(numsites, "all"), 
                       dist.km = c(ifelse(nrow(mins_prop) == 0, NA, mins_prop$dist.km[1]),
                                   ifelse(nrow(mins100) == 0, NA, mins100$dist.km[1])), 
                       sim = c(ifelse(nrow(mins_prop) == 0, NA, mins_prop$sim[1]),
                               ifelse(nrow(mins100) == 0, NA, mins100$sim[1])))
  allmins$rel.dif = allmins$dist.km - allmins$dist.km[2]
  
  return(allmins)
}

#read in files
# soilpts = readOGR("data/soil-shp/", "wrbfu_pts")
# RSG = readxl::read_xlsx("data/WRB_RSGs.xlsx", sheet = "RSG") %>%
#   dplyr::select(RSG, Code, Class, Description, Class_06, Description_06)
# soildf = create_soil_df(soilpts, RSG)
langpts = readOGR("data/lang-shp/", "samerica_lang")
langdf = create_lang_df(langpts)

#specify number of sites used for different ages based on exponential model fit to European dates database
num_sites = function(age, soilpts) {
  #nt = 3.611*10^2*exp((-6.225*10^-5)*age) + 2.79
  nt = 412.5*exp((-4.073*10^-5)*age)
  prop = nt/14384 #proportion of dates in original dataset
  return(floor(nrow(soilpts)*prop))
}

#run through sampling multiple times
minimums10k = get_mins(langpts, langdf, 10000)
minimums50k = get_mins(langpts, langdf, 50000)
minimums100k = get_mins(langpts, langdf, 100000)

for(i in 2:100) { #can change the second number to run a certain number of times
  minimums10k = rbind(minimums, get_mins(soilpts, soildf, 10000))
  minimums50k = rbind(minimums, get_mins(soilpts, soildf, 50000))
  minimums100k = rbind(minimums, get_mins(soilpts, soildf, 100000))
}

write.csv(minimums10k, "output/results/minimums10k.csv")
write.csv(minimums50k, "output/results/minimums50k.csv")
write.csv(minimums100k, "output/results/minimums100k.csv")
