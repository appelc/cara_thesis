###########################
### Get GPS data from Google Drive
### & do basic data exploration
###########################

library(googlesheets)
library(chron)

setwd("C:/Users/Cara/Documents/__RESEARCH/ANALYSIS/033116")

my_sheets <- gs_ls()
my_sheets                               ## lists all your sheets
gps.pts <- gs_title("Porc GPS data")
gps.pts                                 ## lists all worksheets within the file

# this next step will take a few minutes
gps <- data.frame(gs_read(ss=gps.pts, ws="porcGPS", is.na(TRUE), range=cell_cols(1:26)))
colnames(gps)

gps$Animal.ID <- as.factor(gps$Animal.ID)
gps$Date <- as.Date(gps$Date, "%m/%d/%Y")
gps$Time <- chron(times = gps$Time, format = c(times = "h:m:s")) #why is time NA?

hen <- subset(gps, Animal.ID == "15.01")
hen <- droplevels(hen)
unique(hen$Animal.ID)

gps.sp <- SpatialPointsDataFrame(data.frame(gps$Longitude, gps$Latitude),
                                    data=data.frame(gps),
                                    proj4string=CRS())