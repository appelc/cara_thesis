##########
## summary statistics

library(googlesheets)

gs_ls()
locs <- gs_title("Porc relocation data")
porc.vhf <- data.frame(gs_read(ss=locs, ws="Relocations", is.na(TRUE), range=cell_cols(1:8)))
colnames(porc.vhf) <- c("date", "id", "sess", "type", "time", "az", "utm_e", "utm_n")
porc.vhf <- subset(porc.vhf, type %in% c("V","V*","P","P*","L"))
porc.vhf$utm_e <- as.numeric(porc.vhf$utm_e)
porc.vhf$utm_n <- as.numeric(porc.vhf$utm_n)
porc.vhf$date <- as.Date(porc.vhf$date, "%m/%d/%Y") 

sum1.vhf <- porc.vhf[porc.vhf$date < "2015-11-01",]
win.vhf <- porc.vhf[(porc.vhf$date > "2015-10-31") & (porc.vhf$date < "2016-03-01"), ]
sum2.vhf <- porc.vhf[porc.vhf$date > "2016-02-29",]

t1 <- table(sum1.vhf$id)
m1 <- mean(t1)
se1 <- (sd(t1)) / (sqrt(length(t1)))

t2 <- table(win.vhf$id)
m2 <- mean(t2)
se2 <- (sd(t2)) / (sqrt(length(t2)))

t3 <- table(sum2.vhf$id)
m3 <- mean(t3)
se3 <- (sd(t3)) / (sqrt(length(t3)))


## Incorporate GPS data (1 random point per day)
gps.pts <- read.csv("daily.gps.csv")
gps.pts$type <- rep("gps", nrow(gps.pts))
gps.pts$az <- rep(NA, nrow(gps.pts))
gps.pts.df <- data.frame(gps.pts$Date, gps.pts$Animal.ID, gps.pts$Session, gps.pts$type,
                         gps.pts$Time, gps.pts$az, gps.pts$UTM.E, gps.pts$UTM.N)
colnames(gps.pts.df) <- colnames(porc.vhf)
gps.pts.df$date <- as.Date(gps.pts.df$date, "%Y-%m-%d")
gps.pts.df$id <- as.factor(gps.pts.df$id)

## combine VHF with GPS points
porc.locs <- rbind(porc.vhf, gps.pts.df)

## exclude 15.04
porc.locs <- porc.locs[porc.locs$id != '15.04',]

## subset summer locations (before Nov 1 or after March 1) and winter locations
sum.locs <- porc.locs[(porc.locs$date < "2015-11-01") | (porc.locs$date > "2016-02-29"), ]
win.locs <- porc.locs[(porc.locs$date > "2015-10-31") & (porc.locs$date < "2016-03-01"), ]

## Keep only animals with >= 5 locations
n <- table(porc.locs$id)
porc.locs <- subset(porc.locs, id %in% names(n[n >= 5]), drop=TRUE)
porc.locs <- droplevels(porc.locs)

n <- table(sum.locs$id)
sum.locs <- subset(sum.locs, id %in% names(n[n >= 5]), drop=TRUE)
sum.locs <- droplevels(sum.locs)

n <- table(win.locs$id)
win.locs <- subset(win.locs, id %in% names(n[n >= 5]), drop=TRUE)
win.locs <- droplevels(win.locs)
