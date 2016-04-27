### Try KUD with all locations to get a study area boundary

library(googlesheets)

gs_ls()
locs <- gs_title("Porc relocation data")
porc.locs <- data.frame(gs_read(ss=locs, ws="Relocations", is.na(TRUE), range=cell_cols(1:8)))
colnames(porc.locs) <- c("date", "id", "sess", "type", "time", "az", "utm_e", "utm_n")
porc.locs <- subset(porc.locs, type %in% c("V","V*","P","P*","L"))
porc.locs$utm_e <- as.numeric(porc.locs$utm_e)
porc.locs$utm_n <- as.numeric(porc.locs$utm_n)
porc.locs$date <- as.Date(porc.locs$date, "%m/%d/%Y")
porc.locs$name <- rep("collared_porcs", nrow(porc.locs))

## make SPDF for all points on collared animals (not distinguished by individual)
porc.sp <- SpatialPointsDataFrame(data.frame(porc.locs$utm_e, porc.locs$utm_n),
                                  data=data.frame(porc.locs$name),
                                  proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

## calculate grid size and estimate kernelUD
c = 50   ## desired cell size (meters)
fake.kern <- kernelUD(xy = porc.sp, extent = 1)
spdf <- raster(as(fake.kern[[1]], "SpatialPixelsDataFrame"))
eas <- diff(range(spdf@extent[1:2]))
nor <- diff(range(spdf@extent[3:4]))
if(eas > nor){
      g <- (eas/c)
    } else {
      g <- (nor/c)
    }

kernel.all <- kernelUD(xy = porc.sp, h = 60, grid = g, extent = 1)
cont99.all <- getverticeshr.estUD(kernel.all[[1]], percent = 99, unin = "m", unout = "km2", standardize = FALSE)

cont99.clipped <- cont99.all[veg,]
