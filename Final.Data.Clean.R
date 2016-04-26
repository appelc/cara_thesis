
###############  *** GPS DATA CLEANING ALGORITHM ***  ###############

#            v------CSV File Goes Here
data.file <- porc
#           v-------Desired Minimum Threshold For Outliers (meters)
minimum <- 20
#    v-------Hours to remove after first point/before last point
f <- 10 # First
l <- 1  # Last
#   Enter desired OUTLIER csv and shape file name *IN QUOTES  *WITHOUT .CSV/.SHP
out.file <- "16.17_outlier_20.csv"
#   Enter desired NO-OUTLIER csv and shape file name *IN QUOTES  *WITHOUT .CSV/.SHP
no.out.file <- "16.17_clean_20.csv"

#   *** Files will save to working directory ***

library(adehabitatHR)
library(rgdal)
library(sp)

# VVV***Preliminary Setup***VVV
imp.gps <- data.file
imp.sub.gps <- imp.gps[with(imp.gps, Latitude > 41.8333 &  Latitude < 41.95 & Longitude > -124.2166 & Longitude < -124.1333),]
spdf.gps <- SpatialPointsDataFrame(data.frame(imp.sub.gps$Longitude, imp.sub.gps$Latitude), 
                                   proj4string = CRS("+proj=longlat +datum=WGS84"), data=imp.sub.gps)
utm.spdf.gps <- spTransform(spdf.gps, CRS="+proj=utm +zone=10 +datum=NAD83")
utm.gps <- as.data.frame(utm.spdf.gps)
names(utm.gps)[(ncol(utm.gps)-1)] <- "UTM.E"
names(utm.gps)[(ncol(utm.gps))] <- "UTM.N"
length <- nrow(utm.gps)
head.utm.gps <- rbind(utm.gps[-(1:length),])
# ^^^***Preliminary Setup***^^^

# VVV***Remove First and Last points***VVV  
utm.gps.nofl <- NULL
for(i in levels(gps$Session)){
  nofl.sub <- subset(utm.gps, Session==i)
  nofl.sub <- nofl.sub[nofl.sub$posix > (nofl.sub$posix[1] + (f*(60^2))),]
  nofl.sub <- nofl.sub[nofl.sub$posix < (nofl.sub$posix[length(nofl.sub$posix)] - (l*60^2)),]
  utm.gps.nofl <- rbind(utm.gps.nofl, nofl.sub)
}
utm.gps <- utm.gps.nofl
# ^^^***Remove First and Last points***^^^

# VVV***Remove Duplicates***VVV
utm.nd.gps <- rbind(head.utm.gps, utm.gps[1,])
i <- 1
j <- 0
while(i < (nrow(utm.gps)-1)){
  first.pt <- utm.gps[i,]
  second.pt <- utm.gps[i+1,]
  if(!(first.pt$Latitude == second.pt$Latitude && first.pt$Longitude == second.pt$Longitude)){
    no.dup.pt <- data.frame(second.pt)
    utm.nd.gps <-rbind(utm.nd.gps, no.dup.pt)
    }else{
    j <- j + 1
  }
  i <- i + 1
}
# ^^^***Remove Duplicates***^^^

# VVV***Removes Outliers***VVV
clean_gps_data <- function(utm.nd.gps, minimum){
  i <- 1
  while(i <= (nrow(utm.nd.gps)-2)){
    first.pt <- utm.nd.gps[i,]
    second.pt <- utm.nd.gps[i+1,]
    third.pt <- utm.nd.gps[i+2,]
    dist.1.2 <- sqrt(((first.pt$UTM.E-second.pt$UTM.E)^2) + ((first.pt$UTM.N-second.pt$UTM.N)^2))
    dist.1.3 <- sqrt(((first.pt$UTM.E-third.pt$UTM.E)^2) + ((first.pt$UTM.N-third.pt$UTM.N)^2))
    dist.2.3 <- sqrt(((second.pt$UTM.E-third.pt$UTM.E)^2) + ((second.pt$UTM.N-third.pt$UTM.N)^2))
    if(dist.1.3 < dist.1.2 && dist.1.3 < dist.2.3 && dist.1.2 > minimum && dist.2.3 > minimum){
      utm.nd.gps <- utm.nd.gps[-(i+1),]
      i <- i
      } else {
      i <- i + 1
    }
  }
  return(utm.nd.gps)
} 

clean.data <- clean_gps_data(utm.nd.gps, minimum)
# ^^^***Removes Outliers***^^^

# VVV***Count Outliers***VVV
count_outliers <- function(utm.nd.gps, minimum){
  i <- 1
  out.gps <- NULL
  while(i < (nrow(utm.nd.gps)-2)){
    first.pt <- utm.nd.gps[i,]
    second.pt <- utm.nd.gps[i+1,]
    third.pt <- utm.nd.gps[i+2,]
    dist.1.2 <- sqrt(((first.pt$UTM.E-second.pt$UTM.E)^2) + ((first.pt$UTM.N-second.pt$UTM.N)^2))
    dist.1.3 <- sqrt(((first.pt$UTM.E-third.pt$UTM.E)^2) + ((first.pt$UTM.N-third.pt$UTM.N)^2))
    dist.2.3 <- sqrt(((second.pt$UTM.E-third.pt$UTM.E)^2) + ((second.pt$UTM.N-third.pt$UTM.N)^2))
    if(dist.1.3 < dist.1.2 && dist.1.3 < dist.2.3 && dist.1.2 > minimum && dist.2.3 > minimum){
      out.pt <- data.frame(second.pt)
      out.gps <-rbind(out.gps, out.pt)
      utm.nd.gps <- utm.nd.gps[-(i+1),]
      i <- i
    }
    i <- i + 1
  }
  return(nrow(out.gps))
} 


num.outliers <- NULL
minima <- seq(from=1, to=20, by=1)
for(i in 1:length(minima)){
  num.outliers[i] <- count_outliers(utm.nd.gps, minima[i])
}
plot(num.outliers~minima, type="l")
# ^^^***Count Outliers***^^^

write.csv(clean.data, no.out.file)
