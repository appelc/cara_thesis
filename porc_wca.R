## Try weighted compositional analysis
## ala Millspaugh et al. 2006

#install.packages("ruf",repos="http://www.stat.ucla.edu/~handcock")
library(adehabitatHR)
library(googlesheets)
library(raster)
library(rgdal)
library(rgeos)
library(ruf)

######################
## 1. First, load porcupine location data & veg data
######################
gs_ls()
locs <- gs_title("Porc relocation data")
porc.locs <- data.frame(gs_read(ss=locs, ws="Relocations", is.na(TRUE), range=cell_cols(1:8)))
colnames(porc.locs) <- c("date", "id", "sess", "type", "time", "az", "utm_e", "utm_n")
porc.locs <- subset(porc.locs, type %in% c("V","V*","P","P*","L"))
porc.locs$utm_e <- as.numeric(porc.locs$utm_e)
porc.locs$utm_n <- as.numeric(porc.locs$utm_n)
## check date format before running line 51 or 52
#porc.locs$date <- as.Date(porc.locs$date, "%m/%d/%Y") 
#porc.locs$date <- as.Date(porc.locs$date, origin = as.Date("1899-12-30"))

## OPTIONAL: only keep summer locations (before Nov 1)
## could even incorporate this into the for-loop below
sum.locs <- porc.locs[porc.locs$date < "2015-11-01",]

## Turn these into a Spatial Points Data Frame
## Delete... I never actually use these! Other than to assign projection for "veg"
porc.sp <- SpatialPointsDataFrame(data.frame(porc.locs$utm_e, porc.locs$utm_n),
                                  data=data.frame(porc.locs$id),
                                  proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

sum.sp <- SpatialPointsDataFrame(data.frame(sum.locs$utm_e, sum.locs$utm_n),
                                 data=data.frame(sum.locs$id),
                                 proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

## Load veg data
veg <- readOGR(dsn="shapefiles", layer="Veg categories CA", verbose=TRUE)
proj4string(veg) <- proj4string(porc.sp)
veg.ext <- readOGR(dsn="shapefiles", layer="Veg extent new", verbose=TRUE)
proj4string(veg.ext) <- proj4string(veg)

######################
## 2. Then, extract the UD from "adehabitatHR" package
###################### 

## Calculate grid & extent based on desired cell size (# meters on each side)
## For for each animal separately 

## Also calculate KUD based on summer points ONLY, but within grid of the extent
## for all of the points. Then clip to the 99% contour for all the points, as well
## as the veg layer extent.

ids <- unique(sum.locs$id)
ud.list <- list()
ud.summer.list <- list()
ud.clipped.list <- list()
contour.list <- list()
contour.summer.list <- list()
kde.areas <- list()

for (i in ids){
      locs.i <- porc.locs[porc.locs$id == i,]
      locs.i$id_season <- rep(paste(i, "_all", sep = ""), nrow(locs.i))
      locs.sum.i <- sum.locs[sum.locs$id == i,]
      locs.sum.i$id_season <- rep(paste(i, "_sum", sep = ""), nrow(locs.sum.i))
      locs.all.i <- rbind(locs.i, locs.sum.i)
      sp.i <- SpatialPointsDataFrame(data.frame(locs.all.i$utm_e, locs.all.i$utm_n),
                                 data=data.frame(locs.all.i$id_season),
                                 proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
      c = 10   ## desired cell size (meters)
      fake.kern <- kernelUD(xy = sp.i, extent = 1)
      spdf <- raster(as(fake.kern[[1]], "SpatialPixelsDataFrame"))
      eas <- diff(range(spdf@extent[1:2]))
      nor <- diff(range(spdf@extent[3:4]))
      if(eas > nor){
        g <- (eas/c)
      } else {
        g <- (nor/c)
      }
      # calculate UD on both IDs ("all" and "summer") with same4all = TRUE
      kern.i <- kernelUD(xy = sp.i, h = 60, grid = g, extent = 1, same4all = TRUE)
      kde.i <- kernel.area(kern.i, percent = c(50, 90, 95, 99), unin = "m", unout = "km2", standardize = FALSE)
      data.frame(kde.i, row.names = c("50", "90", "95", "99"))
      kde.areas[[i]] <- kde.i
  
      # make 99% contours (full and summer)
      cont99.all.i <- getverticeshr.estUD(kern.i[[1]], percent = 99, unin = "m", unout = "km2", standardize = FALSE)
      cont99.sum.i <- getverticeshr.estUD(kern.i[[2]], percent = 99, unin = "m", unout = "km2", standardize = FALSE)
  
      # clip summer UD to 99% contour from ALL points (not just summer), and veg extent
      sum.ud.i <- (kern.i[[2]])[cont99.all.i,]
      sum.ud.i <- sum.ud.i[veg,]
  
      # save full UD, summer UD, and clipped UD:
      ud.list[[i]] <- kern.i[[1]]
      ud.summer.list[[i]] <- kern.i[[2]]
      ud.clipped.list[[i]] <- sum.ud.i ##it's now a "SpatialPixelsDataFrame"
  
      # and save the contours:
      contour.list[[i]] <- cont99.all.i
      contour.summer.list[[i]] <- cont99.sum.i # don't actually use this for anything later
}

## it's cool to look at a few here:
image(ud.clipped.list[[12]])
plot(veg, add=TRUE)
plot(contour.list[[12]], add=TRUE, border="blue", lwd=2)
plot(contour.summer.list[[12]], add=TRUE, border="green", lwd=2)

## output KDE areas
#write.csv(kde.areas, "csvs/kde_areas_050316.csv")

######################
## 3. Then, create a list of tables with id, coord, and UD height for each porc
##    a. For UD height at each pixel
######################

ids <- unique(sum.locs$id)
height.list <- list()

for(i in ids){
  ud.i <- ud.clipped.list[[i]]
  ud.height.i <- ud.i$ud
  coords.i <- ud.i@coords
  ht.coords.i <- data.frame((rep(i, length(ud.height.i))), ud.height.i, coords.i)
  colnames(ht.coords.i) <- c("id", "height", "x", "y")
  height.list[[i]] <- data.frame(ht.coords.i) 
}

## wireframe plots! better function to get lat/lon or put it on a map?
library(lattice)
wireframe(height ~ x * y, data=height.list[[12]], drape=TRUE, main="15.12 summer UD height")

######################
## 4. Assign values of covariates (veg class, canopy height) to cells
## and include a column for normalizing the UD height: (x - min) / (max - min)
## (and/or log of UD height)
##    a. For UD height at each pixel
######################

## this loop creates a SPDF for each animal, does 'overlay' with veg class, 
## gets rid of cells where veg=NA (there shouldn't be many but they may mess up ruf.fit),
## then turns it back into a data frame for calculating RUF in next step

ids <- unique(sum.locs$id)
final.list <- list()

for (i in ids){
  ht.i <- height.list[[i]]
  spdf.i <- SpatialPointsDataFrame(data.frame(ht.i$x, ht.i$y),
                                   data=data.frame(ht.i$id, ht.i$height),
                                   proj4string = CRS(proj4string(veg)))
  spdf.i@data$veg <- over(spdf.i, veg)$Class_2
  df.i <- data.frame(i, spdf.i@data$ht.i.height, spdf.i@coords, spdf.i@data$veg)
  colnames(df.i) <- c("id", "ud", "x", "y", "veg")
  df.i <- df.i[!is.na(df.i$veg),]
  min <- min(df.i$ud)
  max <- max(df.i$ud)
  df.i$height_norm <- ((df.i$ud) - min) / (max - min)
  df.i$height_log <- log(df.i$ud)      
  final.list[[i]] <- df.i
}

## another cool figure:
plot(spdf.i)
plot(veg, add=TRUE)
plot(contour.list[[14]], add=TRUE, border="blue", lwd=2)
points(utm_n ~ utm_e, data=porc.locs[porc.locs$id == "15.14",], col="red", pch=16)
points(utm_n ~ utm_e, data=sum.locs[sum.locs$id == "15.14",], col="green", pch=16)

######################
## 4. For weighted compositional analysis: 
## sum raw UD values by veg type and divide the summed UD values by the 
## total UD value of all patches to obtain a UD-weighted estimate of use 
## for each habitat type for each individual animal 
## - (Millspaugh et al. 2004 p. 391)
######################

ids <- unique(sum.locs$id)
tables <- list()
for (i in ids){
        ud.i <- final.list[[i]]
        table.i <- aggregate(ud ~ veg, data=ud.i, FUN = sum)
        table.i$ud_weight <- table.i$ud / sum(table.i$ud)
        table.i$log_ud_weight <- log(table.i$ud_weight)
        tables[[i]] <- table.i
}

## now, need to calcluate "log-transformed availability data"
veg.99kdes <- list()
veg.areas <- list()
for (i in ids){
        cont99.i <- contour.list[[i]]
        veg.i <- intersect(veg, cont99.i)
        veg.i <- veg.i[!is.na(veg.i@data$Class_2),] #get rid of NAs
        area.all <- gArea(veg.i, byid = TRUE) #units should be m^2
        veg.df.i <- data.frame(veg.i$Class_2, area.all)
        colnames(veg.df.i) <- c("veg", "area")
        veg.areas.i <- aggregate(area ~ veg, data=veg.df.i, FUN = sum)
        veg.areas.i$prop_area <- veg.areas.i$area / sum(veg.areas.i$area)
        veg.99kdes[[i]] <- veg.i
        veg.areas[[i]] <- veg.areas.i
}

## get 7 warnings: "In RGEOSUnaryPredFunc(spgeom, byid, "rgeos_isvalid") :
##  Ring Self-intersection at or near point x, y 
##  (I think just because of fragments/holes between polygons I digitized)

## There's also a problem in the "intersect" step where I think it doesn't 
## incldue some polygons instead of truncating them within the 99% contour
## (see figures below for an example)

## cool figures but I know it's really messy (fix later):
ids <- unique(sum.locs$id)
for (i in ids){
        veg.i <- veg.99kdes[[i]]
        mypath <- file.path("figures", "kdes_with_veg", paste(i, "_veg_99kde", ".png", sep = ""))
        png(file=mypath)
        mytitle = paste("99% KDE ", i, sep = "")
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(veg.i, main = mytitle)
        leg.txt <- sort(unique(veg$Class_2))
        leg.col <- c("khaki1", "khaki3", "aquamarine", "khaki4", "darkolivegreen4", "cadetblue1",
                     "coral1", "yellow3", "darkolivegreen3", "darkseagreen3", "aquamarine4")
        legend("topright", inset=c(-0.3,0), legend = leg.txt, pch = 15, col = leg.col, cex=0.9)
        plot(veg.i[veg.i$Class_2 == "Beach",], add=TRUE, col="khaki1")
        plot(veg.i[veg.i$Class_2 == "Beachgrass dune",], add=TRUE, col="khaki3")
        plot(veg.i[veg.i$Class_2 == "Brackish marsh",], add=TRUE, col="aquamarine")
        plot(veg.i[veg.i$Class_2 == "Coastal scrub",], add=TRUE, col="khaki4")
        plot(veg.i[veg.i$Class_2 == "Conifer forest",], add=TRUE, col="darkolivegreen4")
        plot(veg.i[veg.i$Class_2 == "Freshwater marsh",], add=TRUE, col="cadetblue1")
        plot(veg.i[veg.i$Class_2 == "Fruit tree",], add=TRUE, col="coral1")
        plot(veg.i[veg.i$Class_2 == "Meadow",], add=TRUE, col="yellow3")
        plot(veg.i[veg.i$Class_2 == "Pasture",], add=TRUE, col="darkolivegreen3")
        plot(veg.i[veg.i$Class_2 == "Shrub swale",], add=TRUE, col="darkseagreen3")
        plot(veg.i[veg.i$Class_2 == "Wooded swale",], add=TRUE, col="aquamarine4")
        plot(sum.sp[sum.sp$sum.locs.id == i,], add=TRUE, pch=16, cex=1, col="red")
        dev.off() 
}

## may need to run this again to be able to plot again:
#dev.off()

######################
## 5. Subtract differences in log-transformed availability data from the
##    log-transformed use data for each animal and then test for overall
##    selection using Wilks' lambda
######################

## combine used and avail in the same table
## we have "tables" (a list) and "veg.areas" (a list)

ids <- unique(sum.locs$id)
full.table <- NULL
final.table <- NULL
for (i in ids){
        tables.i <- tables[[i]]
        veg.areas.i <- veg.areas[[i]]       
        tables.i$area <- veg.areas.i$area
        tables.i$prop_area <- veg.areas.i$prop_area
        tables.i$id <- rep(i, nrow(tables.i))
        tables.i$sel <- tables.i$ud_weight - tables.i$prop_area
        full.table <- rbind(full.table, tables.i)
        final.df <- data.frame(tables.i$id, tables.i$veg, tables.i$log_ud_weight,
                               tables.i$prop_area, tables.i$sel)
        colnames(final.df) <- c("id", "veg", "log_ud_wt", "prop_area", "sel")
        final.table <- rbind(final.table, final.df)
}

write.csv(final.table, "csvs/wt_comp_analysis_050316.csv")


