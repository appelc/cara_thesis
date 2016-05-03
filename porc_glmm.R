#############################################################
### Taking a stab at fitting a GLMM in WinBugs with porc data
### based on notes from STAT 510 lab (Mar 28)

## the lab does have code for simulating binomial data
## (would need to do multinomial, actually)
#############################################################

library(googlesheets)
library(rgdal)
library(raster)
library(rgeos)
library(data.table)
library(R2WinBUGS)
bugs.dir <- "C:/Programs/WinBUGS14"

## import data
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

## make a column to say that these are all collared porc locations
## (we're not going to use individual ID)
porc.locs$name <- rep("collared_porcs", nrow(porc.locs))

## make SPDF for all points on collared animals (not distinguished by individual)
all.loc.sp <- SpatialPointsDataFrame(data.frame(porc.locs$utm_e, porc.locs$utm_n),
                                  data=data.frame(porc.locs$name),
                                  proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

## Now make a subset of individual porcupine locations during summer only (before Nov 1)
sum.locs <- subset(porc.locs, date < "2015-11-01")
sum.sp <- SpatialPointsDataFrame(data.frame(sum.locs$utm_e, sum.locs$utm_n),
                                 data=data.frame(sum.locs$id),
                                 proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

## Load veg data
veg <- readOGR(dsn="shapefiles", layer="Veg categories CA", verbose=TRUE)
proj4string(veg) <- proj4string(all.loc.sp)

## I also have a dissolved version of the veg shapefile to use for clipping
veg_ext <- readOGR(dsn="shapefiles", layer="Veg extent", verbose=TRUE)
proj4string(veg_ext) <- proj4string(all.loc.sp)

## Use 99% KDE for all animals as the extent for selecting "available" points
## first, calculate grid size and estimate kernelUD
c = 10   ## desired cell size (meters)
fake.kern <- kernelUD(xy = all.loc.sp, extent = 1)
spdf <- raster(as(fake.kern[[1]], "SpatialPixelsDataFrame"))
eas <- diff(range(spdf@extent[1:2]))
nor <- diff(range(spdf@extent[3:4]))
if(eas > nor){
  g <- (eas/c)
} else {
  g <- (nor/c)
}

## Do 99% KDE on ALL points, not just summer (is this right?)
kernel.all <- kernelUD(xy = all.loc.sp, h = 60, grid = g, extent = 1)
image(kernel.all) #how do they look?
cont99.all <- getverticeshr.estUD(kernel.all[[1]], percent = 99, unin = "m", unout = "km2", standardize = FALSE)

## clip contour to "veg" extent
cont99.clipped <- gIntersection(veg_ext, cont99.all, byid = FALSE, drop_lower_td = TRUE)

## so "cont99.clipped" is the extent within which we will sample "available" points

## also clip "sum.sp" to veg_extent (there are a few points outside)
sum.sp.clipped <- sum.sp[veg_ext,]

## Create random points within thestudy area (clipped 99% contour for all animals)
## use "nrow()" to make sure we have the same no. of avail. pts as used pts (during the summer)
sum.avail.pts <- spsample(cont99.clipped, n=(nrow(sum.sp.clipped)), type="random")

## plot everything:
plot(veg_ext) # need to fix annoying holes between polygons...
plot(cont99.clipped, add=TRUE, col="green")
plot(sum.sp.clipped, add=TRUE, pch=16, cex=0.4)
plot(sum.avail.pts, add=TRUE, col="red", pch=16, cex=0.4)

## Now add 0/1 column for pres/avail
sum.pts <- data.frame(sum.sp.clipped$sum.locs.id, sum.sp.clipped@coords) # convert to data.frame
sum.pts$pres <- rep(1, nrow(sum.pts)) # add a 1 for all presence (used) pts
sum.pts$pres_name <- rep("present", nrow(sum.pts)) # add a column for description
colnames(sum.pts) <- c("id", "x", "y", "pres", "pres_name")

## Do the same for the "available" points we created
avail.pts <- data.frame(sum.avail.pts@coords)
avail.pts$pres <- rep(0, nrow(avail.pts)) # add a 0 for all available pts
avail.pts$pres_name <- rep("available", nrow(avail.pts))
avail.pts$id <- rep("avail", nrow(avail.pts))
avail.pts <- avail.pts[,c("id", "x", "y", "pres", "pres_name")] # rearrange to match "sum.pts"

## combine used and available pts
sum.used.avail <- rbind(sum.pts, avail.pts)
nrow(sum.used.avail) # good, it should be the number of summer points (clipped) * 2

## now overlay "veg"
## need to make SPDF again... should have done this before?
points.sp <- SpatialPointsDataFrame(data.frame(sum.used.avail$x, sum.used.avail$y),
                                            data = data.frame(sum.used.avail$id, sum.used.avail$pres),
                                            proj4string = CRS(proj4string(all.loc.sp)))

points.sp$veg <- over(points.sp, veg)$Class_2
points.sp <- points.sp[!is.na(points.sp$veg),] # get rid of veg=NA

points.df <- data.frame(points.sp@data)
colnames(points.df) <- c("id", "pres", "veg")

data.table(points.df)



sumpts.test <- aggregate.data.frame(points.sp, list("id", "veg"), FUN=sum)
test <- table(points.df)
test <- aggregate(cbind(pres, veg) ~ id, data=points.df, FUN = sum)

for (i in veg){
        i
}

hen.scrub <- points.df[points.df$id == "15.01" & points.df$veg == "Coastal scrub",]
hen.scrub.sum <- sum(hen.scrub$pres)

hen.scrub <- sum(hen.scrub$pres)
hen.scrub

###############

sink("BinomGLM1.txt")
cat("
    model {
    
    # Priors
    beta0 ~ dnorm(0, .00001) # vague prior
    beta1 ~ dnorm(0, .00001) # Prior for slope 
     
    # Likelihood: note key components in one line each
    for (i in 1:n){
    Success[i] ~ dbin(p[i], N[i])          # r.v. (draw from binom w/prob success specific to that location and number of times we tried it at that location)
    logit(p[i]) <- beta0 +  beta1*x1[i]    # linear predictor w/link function
    }
    
    ## calculate deviance 'by-hand' to compare with WinBUGS' deviance
    for (i in 1:n){
    LL[i] <- logfact(N[i]) - logfact(Success[i]) -
    logfact(N[i] - Success[i]) +
    Success[i]*log(p[i]) +
    (N[i] - Success[i])*log(1-p[i])
    }
    
    dev <- -2*sum(LL[])
    
    
    # Now, create deviance for 'ideal' datasets
    for (i in 1:n){
    Pred[i] ~ dbin(p[i], N[i])
    LLP[i] <- logfact(N[i]) - logfact(Pred[i]) -
    logfact(N[i] - Pred[i]) +
    Pred[i]*log(p[i]) +
    (N[i] - Pred[i])*log(1-p[i])
    }
    
    devP <- -2*sum(LLP[])
    test <- step(dev - devP)
    bpvalue <- mean(test)
    
    
    }
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(Success = as.numeric(Sim1),
                 n = as.numeric(length(Sim1)),
                 N = as.numeric(N.trials),
                 x1 = as.numeric(x1))


# Initial values
inits <- function() list(beta0 = runif(1, -2, 2), beta1 = runif(1, -2, 2))

# Parameters monitored
params <- c("beta0", "beta1", "dev","devP", "bpvalue")

# MCMC settings
ni <- 50000 # need 50000 iterations total
nt <- 15 # thin (keep 1 out of every 15)
nb <- 20000 # burn in for 20000
nc <- 3 # 3 chains

# Call WinBUGS from R (BRT <1 min)
out <- bugs(win.data, inits, params, "BinomGLM1.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
