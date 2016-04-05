################################
### Get GPS data from Google Drive
### & do basic data exploration
###############################

library(googlesheets)
library(chron)
library(sp)

setwd("C:/Users/Cara/Documents/__RESEARCH/ANALYSIS/040416")

my_sheets <- gs_ls()                 
my_sheets                            
gps.pts <- gs_title("Porc GPS data")
gps.pts

# this next step will take a few minutes
gps <- data.frame(gs_read(ss=gps.pts, ws="porcGPS", is.na(TRUE), range=cell_cols(1:26)))
head(gps)

# check class() for the important columns and change accordingly
gps$Animal.ID <- as.factor(gps$Animal.ID)

#gps$Date <- as.Date(gps$Date, "%m/%d/%Y")
gps$Date <- as.Date("1900-01-01") + gps$Date
gps$Time <- chron(times = gps$Time, format = c(times = "h:m:s")) #why is time NA?
gps$posix <- paste(gps$Date, gps$Time)
gps$posix <- as.POSIXct(strptime(gps$posix, "%Y-%m-%d %H:%M:%S"), tz="America/Los_Angeles")


# example to get data for one animal
hen <- subset(gps, Animal.ID == "15.01")
hen <- droplevels(hen)
unique(hen$Animal.ID)

# convert to spatial points data frame
hen.sp <- SpatialPointsDataFrame(data.frame(hen$Longitude, hen$Latitude),
                                 data=data.frame(hen),
                                 proj4string=CRS("+proj=longlat +datum=WGS84"))
# is this coordinate system correct?
plot(hen)

# spdf for all animals
gps.sp <- SpatialPointsDataFrame(data.frame(gps$Longitude, gps$Latitude),
                                    data=data.frame(gps),
                                    proj4string=CRS("+proj=longlat +datum=WGS84"))
plot(gps.sp)

