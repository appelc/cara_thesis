library(googlesheets)
library(sp)
library(chron)

#### First, read data from googlesheets
gs_ls()                 
gps.data <- gs_title("Porc GPS data")
## "gs_read_csv" is MUCH faster than "gs_read" but I can't get range=cell_cols(1:7) to work
gps <- data.frame(gs_read_csv(ss=gps.data, ws="porcGPS", is.na(TRUE), stringsAsFactors=FALSE, range = cell_cols(1:7)))
head(gps)
#gps <- gps[!is.na(gps$Animal.ID),]
gps$Animal.ID[gps$Animal.ID == 16.2] <- '16.20' ## why do I have to do this?
gps$Animal.ID <- as.factor(gps$Animal.ID)
gps$Date <- as.Date(gps$Date, "%m/%d/%Y")
#gps$Date <- as.Date("1900-01-01") + gps$Date
test.date <- as.character(gps$Date)
test.time <- as.character(gps$Time)
#gps$Time <- chron(times = gps$Time, format = c(times = "%I:%M:%S %p"))
posix.test <- as.POSIXct(strptime(paste(test.date, test.time), "%Y-%m-%d %I:%M:%S %p"), tz="America/Los_Angeles")
posix.test.pdt <- as.POSIXct(format(posix.test, tz="America/Los_Angeles", usetz=TRUE))
gps$posix <- posix.test.pdt
gps <- gps[,-c(8:13)] ## get rid of extra columns
head(gps)
tail(gps)

###############  *** GPS DATA CLEANING ALGORITHM ***  ###############
#gps <- read.csv('D:/GIS DATA/Porc_GPS/16.20/16.20_gt120_090916.csv') ## or read here instead of Google Drive

#            v------CSV File Goes Here
data.file <- gps
#           v-------Desired Minimum Threshold For Outliers (meters)
minimum <- 20
#    v-------Hours to remove after first point/before last point
f <- 8 # First (I changed this from 10, could probably go lower than 8 even -CA)
l <- 0  # Last (change to 1 if battery not dead when recovered)
#   Enter desired OUTLIER csv and shape file name *IN QUOTES  *WITHOUT .CSV/.SHP
out.file <- "test_all_outlier_20.csv"
#   Enter desired NO-OUTLIER csv and shape file name *IN QUOTES  *WITHOUT .CSV/.SHP
no.out.file <- "test_all_clean_20.csv"

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

## Added by Cara if using from Google Drive:
utm.gps$ID_Session <- paste(utm.gps$Animal.ID, utm.gps$Session, sep = '_')
utm.gps$ID_Session <- as.factor(utm.gps$ID_Session)

# VVV***Remove First and Last points***VVV  
utm.gps.nofl <- NULL
for(i in levels(utm.gps$ID_Session)){
  nofl.sub <- subset(utm.gps, ID_Session==i)
  nofl.sub$ID_Session <- droplevels(nofl.sub$ID_Session) ## added by Cara
  nofl.sub <- nofl.sub[nofl.sub$posix > (nofl.sub$posix[1] + (f*(60^2))),]
  nofl.sub <- nofl.sub[nofl.sub$posix < (nofl.sub$posix[length(nofl.sub$posix)] - (l*60^2)),]
  utm.gps.nofl <- rbind(utm.gps.nofl, nofl.sub)
}
utm.gps <- utm.gps.nofl
# ^^^***Remove First and Last points***^^^

# VVV***Remove Duplicates***VVV
utm.nd.gps <- rbind(head.utm.gps, utm.gps[1,]) ## necessary?
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

## From Cara: what does remove duplicates do? Mostly important for stationary units?
## I'll just skip it and continue from here
#utm.nd.gps <- utm.gps

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

### 
## now what do we have? "clean.data" is final data.frame
unique(clean.data$ID_Session)

## try with just one animal/session
copper1 <- clean.data[clean.data$ID_Session == "16.20_1",]
copper1$ID_Session <- droplevels(copper1$ID_Session)
copper1.spdf <- SpatialPointsDataFrame((data.frame(copper1$UTM.N, copper1$UTM.E)), data = copper1,
                                    proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

dodger1 <- clean.data[clean.data$ID_Session == '16.19_1',]
dodger1$ID_Session <- droplevels(dodger1$ID_Session)
dodger1.spdf <- SpatialPointsDataFrame((data.frame(dodger1$UTM.N, dodger1$UTM.E)), data = dodger1,
                                       proj4string = CRS('+proj=utm +zone=10 +datum=NAD83'))

## how does it compare to before removing outliers?
copper1.outliers <- utm.gps[utm.gps$ID_Session == "16.20_1",]
copper1.outliers$ID_Session <- droplevels(copper1.outliers$ID_Session)
copper1.outliers.spdf <- SpatialPointsDataFrame((data.frame(copper1.outliers$UTM.N, copper1.outliers$UTM.E)),
                                             data = copper1.outliers,
                                             proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
plot(copper1.outliers.spdf, col="red")
plot(copper1.spdf, add=TRUE, col="black")

dodger1.outliers <- utm.gps[utm.gps$ID_Session == '16.19_1',]
dodger1.outliers$ID_Session <- droplevels(dodger1.outliers$ID_Session)
dodger1.outliers.spdf <- SpatialPointsDataFrame((data.frame(dodger1.outliers$UTM.N, dodger1.outliers$UTM.E)),
                                                data = dodger1.outliers,
                                                proj4string = CRS('+proj=utm +zone=10 +datum=NAD83'))
plot(dodger1.outliers.spdf, col = 'red')
plot(dodger1.spdf, add = TRUE, col = 'black')

## now sample 1 per 24 hrs:
#1. subset for ID_Session
#2. get unique date names
#3. sample 1 in each unique date name

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

write.csv(daily.gps, 'csvs/daily_gps_090416.csv') 
daily.gps.spdf <- SpatialPointsDataFrame(data.frame(daily.gps$UTM.E, daily.gps$UTM.N),
                                         data = daily.gps,
                                         proj4string = CRS('+proj=utm +zone=10 +datum=NAD83'))
writeOGR(daily.gps.spdf, dsn = 'Shapefiles', layer = 'daily_gps_090416', driver="ESRI Shapefile")

## Also create shapefile with all cleaned GPS points, not just one per day
write.csv(clean.data, 'csvs/all_gps_090416.csv') 
all.gps.spdf <- SpatialPointsDataFrame(data.frame(clean.data$UTM.E, clean.data$UTM.N),
                                         data = clean.data,
                                         proj4string = CRS('+proj=utm +zone=10 +datum=NAD83'))
all.gps.spdf <- all.gps.spdf[veg,]
writeOGR(all.gps.spdf, dsn = 'Shapefiles', layer = 'all_gps_090416', driver = 'ESRI Shapefile')

## To create KML from specific ID_Session:
library(tlocoh)

my.data <- clean.data[clean.data$ID_Session == '16.19_1',]
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
lxy.exp.kml(porc.lxy, "D:/GIS DATA/Porc_GPS/16.19/16.19_090316.kml")
