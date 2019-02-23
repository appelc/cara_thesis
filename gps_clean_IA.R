library(googlesheets)
library(sp)
library(chron)

#### First, read data from googlesheets 
## (SKIP TO LINE 34 IF USING A CSV, BUT CHECK CORRECT CLASSES FOR EACH COLUMN)

gs_ls()                 
gps.data <- gs_title("Porc GPS data (2)")
## "gs_read_csv" is MUCH faster than "gs_read" but I can't get range=cell_cols(1:7) to work
gps <- data.frame(gs_read_csv(ss=gps.data, ws="PorcGPS", is.na(TRUE), stringsAsFactors=FALSE, range = cell_cols(1:7)))
head(gps)

#gps <- gps[!is.na(gps$Animal.ID),]
gps$Animal.ID[gps$Animal.ID == 16.2] <- '16.20' ## why do I have to do this?
gps$Animal.ID <- as.factor(gps$Animal.ID)

gps$Date <- as.Date(gps$Date, "%m/%d/%Y") ## we were having trouble with date formatting so may need to try both
#gps$Date <- as.Date("1900-01-01") + gps$Date

test.date <- as.character(gps$Date)
test.time <- as.character(gps$Time)
#gps$Time <- chron(times = gps$Time, format = c(times = "%I:%M:%S %p"))

posix.test <- as.POSIXct(strptime(paste(test.date, test.time), "%Y-%m-%d %I:%M:%S %p"), tz="America/Los_Angeles")
posix.test.pdt <- as.POSIXct(format(posix.test, tz="America/Los_Angeles", usetz=TRUE))
gps$posix <- posix.test.pdt

gps <- gps[,-c(8:13)] ## get rid of extra columns
head(gps)

###############  *** GPS DATA CLEANING ALGORITHM ***  ###############

## can read a CSV here instead of using Google Drive import above:
#gps <- read.csv('D:/GIS DATA/Porc_GPS/16.20/16.20_gt120_090916.csv') 

#            v------CSV File Name Goes Here
data.file <- gps
#           v-------Desired Minimum Threshold For Outliers (meters)
minimum <- 20
#    v-------Hours to remove after first point/before last point
f <- 8 # First (I changed this from 10, could probably go lower than 8 even -CA)
l <- 0  # Last (change to 1 if battery not dead when recovered)
#   Enter desired OUTLIER csv and shape file name *IN QUOTES  *WITHOUT .CSV/.SHP
out.file <- "PB_outliers_012418.csv"
#   Enter desired NO-OUTLIER csv and shape file name *IN QUOTES  *WITHOUT .CSV/.SHP
no.out.file <- "PB_cleanGPS_012418.csv"

#   *** Files will save to working directory so double-check***
getwd() 
#setwd(/Users/...)  ## if necessary, input file path here

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

## Added by CA if using from Google Drive:
utm.gps$ID_Session <- paste(utm.gps$Animal.ID, utm.gps$Session, sep = '_')
utm.gps$ID_Session <- as.factor(utm.gps$ID_Session)

# VVV***Remove First and Last points***VVV  
utm.gps.nofl <- NULL
for(i in levels(utm.gps$ID_Session)){
  nofl.sub <- subset(utm.gps, ID_Session==i)
  nofl.sub$ID_Session <- droplevels(nofl.sub$ID_Session) ## added by CA
  nofl.sub <- nofl.sub[nofl.sub$posix > (nofl.sub$posix[1] + (f*(60^2))),]
  nofl.sub <- nofl.sub[nofl.sub$posix < (nofl.sub$posix[length(nofl.sub$posix)] - (l*60^2)),]
  utm.gps.nofl <- rbind(utm.gps.nofl, nofl.sub)
}
utm.gps <- utm.gps.nofl
# ^^^***Remove First and Last points***^^^

# VVV***Remove Duplicates***VVV
utm.nd.gps <- rbind(head.utm.gps, utm.gps[1,]) ## necessary? -CA
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

## From CA: what does remove duplicates do? Mostly important for stationary units?
## I'll just skip it and continue from here
utm.nd.gps <- utm.gps

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

################
## now "clean.data" is final data.frame
################

## PLAYING AROUND WITH SOME VISUALIZATIONS TO SEE HOW IT WORKED - CA:

## try plotting points from just one animal/session (use 'unique(clean.data$ID_Session') to see options)
copper1 <- clean.data[clean.data$ID_Session == "16.20_1",]
  copper1$ID_Session <- droplevels(copper1$ID_Session)
  copper1.spdf <- SpatialPointsDataFrame((data.frame(copper1$UTM.E, copper1$UTM.N)), data = copper1,
                                    proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

dodger3 <- clean.data[clean.data$ID_Session == "16.19_3",]
  dodger3$ID_Session <- droplevels(dodger3$ID_Session)
  dodger3.spdf <- SpatialPointsDataFrame((data.frame(dodger3$UTM.E, dodger3$UTM.N)), data = dodger3,
                                         proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

  
## how does it compare to before removing outliers?
copper1.outliers <- utm.gps[utm.gps$ID_Session == "16.20_1",]
  copper1.outliers$ID_Session <- droplevels(copper1.outliers$ID_Session)
  copper1.outliers.spdf <- SpatialPointsDataFrame((data.frame(copper1.outliers$UTM.E, copper1.outliers$UTM.N)),
                                             data = copper1.outliers,
                                             proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

dodger3.outliers <- utm.gps[utm.gps$ID_Session == "16.19_3",]
  dodger3.outliers$ID_Session <- droplevels(dodger3.outliers$ID_Session)
  dodger3.outliers.spdf <- SpatialPointsDataFrame((data.frame(dodger3.outliers$UTM.E, dodger3.outliers$UTM.N)),
                                                data = dodger3.outliers,
                                                proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

plot(copper1.outliers.spdf, add=TRUE, col="black")
  plot(copper1.spdf, col="red")

plot(dodger3.outliers.spdf, add=TRUE, col="black")
  plot(dodger3.spdf, col="red")

## write shapefiles if you want
writeOGR(copper1.outliers.spdf, dsn = 'Shapefiles/1620_gps_outliers', layer = '1620_gps_outliers', driver = 'ESRI Shapefile')
writeOGR(dodger2.spdf, dsn = 'Shapefiles/1619_2_gps', layer = '1619_gps_2', driver = 'ESRI Shapefile')
writeOGR(dodger3.spdf, dsn = 'Shapefiles/1619_3_gps', layer = '1619_gps_3', driver = 'ESRI Shapefile')


##################
## NOW SAMPLE 1 POINT PER 24 HOURS FOR EACH ANIMAL/SESSION:
##################

#1. subset for ID_Session
#2. get unique date names
#3. sample 1 point in each unique date name

daily.gps <- NULL
gps.samples <- NULL

for (i in levels(clean.data$ID_Session)){
        gps.i <- clean.data[clean.data$ID_Session == i,]
        gps.i$ID_Session <- droplevels(gps.i$ID_Session)
        gps.samples <- NULL
      for (j in unique(gps.i$Date)){
          date.j <- gps.i[gps.i$Date == j,]
          sample <- sample(1:nrow(date.j), 1)
          sample.j <- date.j[sample,]
          gps.samples <- rbind(gps.samples, sample.j)
    }
      daily.gps <- rbind(daily.gps, gps.samples)
}

##################
### EXPORT CSVS AND SHAPEFILES
##################

## one point per day (change filenames):

  write.csv(daily.gps, 'csvs/daily_gps_PB_012418.csv') 
  
  daily.gps.spdf <- SpatialPointsDataFrame(data.frame(daily.gps$UTM.E, daily.gps$UTM.N),
                                           data = daily.gps,
                                           proj4string = CRS('+proj=utm +zone=10 +datum=NAD83'))
  
  writeOGR(daily.gps.spdf, dsn = 'Shapefiles', layer = 'daily_gps_PB_012418', driver="ESRI Shapefile")


## all cleaned GPS points, if desired (change filenames):
  
  write.csv(clean.data, 'csvs/all_gps_PB_012418.csv') 
  
  all.gps.spdf <- SpatialPointsDataFrame(data.frame(clean.data$UTM.E, clean.data$UTM.N),
                                           data = clean.data,
                                           proj4string = CRS('+proj=utm +zone=10 +datum=NAD83'))
  all.gps.spdf <- all.gps.spdf[veg,]
  writeOGR(all.gps.spdf, dsn = 'Shapefiles', layer = 'all_gps_PB_012418', driver = 'ESRI Shapefile')


###########################
## To create KML from specific ID_Session (for viewing in Google Earth):
###########################

library(tlocoh)

my.data <- clean.data[clean.data$ID_Session == '16.19_3',]
test.date <- as.character(my.data$Date)
test.time <- as.character(my.data$Time)
posix.test <- as.POSIXct(strptime(paste(test.date, test.time), "%m/%d/%Y %H:%M:%S"), tz="GMT")
posix.test.pdt <- as.POSIXct(format(posix.test, tz="America/Los_Angeles", usetz=TRUE))
my.data$real.date <- posix.test.pdt

new.gps <- my.data
new.sp <- SpatialPointsDataFrame(data.frame(new.gps$Longitude, new.gps$Latitude), 
                                 proj4string = CRS("+proj=longlat +datum=WGS84"), data=new.gps)
new.sp.utm <- spTransform(new.sp, CRS("+proj=utm +zone=10 +datum=NAD83"))

porc.lxy <- xyt.lxy(xy=new.sp.utm@coords, dt=new.sp.utm$posix,
                    proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"),
                    id="porc", show.bad.timestamps = TRUE)

lxy.exp.kml(porc.lxy, "D:/GIS DATA/Porc_GPS/16.19/16.19_072617.kml")

############
