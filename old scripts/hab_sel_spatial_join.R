##########################
## Porcupine habitat selection
## Try spatial join in R to get used/available PER ANIMAL
##########################

library(rgdal)
library(sp)
library(googlesheets)
library(maps)
library(httr)
library(readr)
library(devtools)

devtools::install_github("jennybc/googlesheets")

## load location points
## pull out only collared animals, keep only visual/patch/LOAS, and format correctly
my_sheets <- gs_ls()
locs <- gs_title("Porc relocation data")
my_sheets
porc.locs.all <- data.frame(gs_read(ss=locs, ws="Relocations", is.na(TRUE)))
## can add "range=cell_cols(1:12)" if I don't want everything

colnames(porc.locs.all)
colnames(porc.locs.all) <- c("date", "id", "sess", "type", "time", "az", "utm_e", "utm_n", 
                             "obs", "loc", "pos", "notes", "xvar", "yvar", "cov", "error")
## c <- grepl('^1', porc.locs.all$id) ## ask which ones start with "1" (TRUE/FALSE)
## then: porc.locs <- porc.loc.all[c,] ## so only keep the rows where it's true
## OR: change "type" for new porc to something like "N" ... go with this for now
porc.locs <- subset(porc.locs.all, type %in% c("V","V*","P","P*","L"))
unique(porc.locs$id)
## gs_read_csv correctly reads numeric and characters, but still need to format date
porc.locs$date <- as.Date(porc.locs$date, "%m/%d/%Y")

## now make it spatial points data frame
locs <- SpatialPointsDataFrame(data.frame(porc.locs$utm_e, porc.locs$utm_n),
                                    data=data.frame(porc.locs),
                                    proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
locs@data ## can see the attributes
locs@coords ## here are the coordinates

## load veg polygons shapefile
?readOGR
veg <- readOGR(dsn="D:/GIS DATA/Veg map", layer="Veg categories CA", verbose=TRUE)
proj4string(veg) <- proj4string(locs) ## it's almost the same but not exactly; needs to be for "over"
is(veg) # it's a spatial polygons data frame
plot(veg) # cool!

## do spatial join using package "sp"
?over

locs@data$class <- over(locs, veg)$Class
locs@data$class2 <- over(locs, veg)$Class_2
locs@data$area <- over(locs, veg)$Area_1

head(locs)
View(locs)

plot(coordinates(locs))
map("world", region="usa", add=TRUE)
plot(veg, border="green", add=TRUE)
legend("topright", cex=0.85,
       c("Bear in park", "Bear not in park", "Park boundary"),
       pch=c(16, 1, NA), lty=c(NA, NA, 1),
       col=c("red", "grey", "green"), bty="n")
title(expression(paste(italic("Ursus arctos"),
                       " sightings with respect to national parks")))

# now plot bear points with separate colors inside and outside of parks
points(bears[!inside.park, ], pch=1, col="gray")
points(bears[inside.park, ], pch=16, col="red")

## or, to keep everything in addition to Class (ID, Area_1, Class_2)
overlay <- cbind(locs, over(locs,veg)) # this didn't work...
head(overlay)

## great! need to figure out what to do with NAs
## (points not in a veg polygon, like Puck, Henrietta's capture, ?love, etc.)

## now can do a table to sum locations in each veg type per animal
## for design ii and iii analysis in adehabitatHS
