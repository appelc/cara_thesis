## generate KML file from GPS points

library(tlocoh)

my.data <- read.csv("D:/GIS DATA/Porc_GPS/15.14/15.14_igotu120_070216/15.14_gt120_070216_cleaned.csv")

test.date <- as.character(my.data$Date)
test.time <- as.character(my.data$Time)
posix.test <- as.POSIXct(strptime(paste(test.date, test.time), "%m/%d/%Y %H:%M:%S"), tz="GMT")
posix.test.pdt <- as.POSIXct(format(posix.test, tz="America/Los_Angeles", usetz=TRUE))
my.data$real.date <- posix.test.pdt

new.gps <- my.data
new.sp <- SpatialPointsDataFrame(data.frame(new.gps$Longitude, new.gps$Latitude), 
                                 proj4string = CRS("+proj=longlat +datum=WGS84"), data=new.gps)
new.sp.utm <- spTransform(new.sp, CRS("+proj=utm +zone=10 +datum=NAD83"))
porc.lxy <- xyt.lxy(xy=new.sp.utm@coords, dt=new.sp.utm$real.date,
                    proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"),
                    id="porc")
lxy.exp.kml(porc.lxy, "D:/GIS DATA/Porc_GPS/15.14/15.14_070214.kml")
