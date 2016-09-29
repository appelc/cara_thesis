library(adehabitatHS)

data(puech)
names(puech)

maps <- puech$maps
locs <- puech$relocations
image(maps)
points(locs, col = as.numeric(slot(locs, 'data')[,1]), pch = 16)

cp <- count.points(locs, maps)
mimage(cp)

X <- slot(maps, 'data')
U <- slot(cp, 'data')

head(X)
head(U)

## create discrete slope categories
slope <- maps[,8]
sl <- slot(slope, 'data')[,1]
av <- factor(cut(sl, c(-0.1, 2, 5, 12, 50)),
             labels = c('Low', 'Medium', 'High', 'Very High'))
(tav <- table(av))

slot(slope, 'data')[,1] <- as.numeric(av)
image(slope)

## how many locations in each category, per individual?
us <- join(locs, slope)
tus <- table(slot(locs, 'data')[,1], us)
class(tus) <- NULL
tus <- as.data.frame(tus)
colnames(tus) <- names(tav)

### SECOND ORDER
## compositional analysis
## 'tav' has the # of resource units in each category
## 'tus' has the # of locations of each animal in each category
tav2 <- matrix(rep(tav, nrow(tus)), nrow = nrow(tus), byrow = TRUE) ## matrix of avail. data
colnames(tav2) <- names(tav)

compana(tus, tav2, test = 'randomisation', rnv = 0.01, nrep = 500, alpha = 0.1)
## not sig. diff. from random
## how are they calculating Lambda? I looked at the function compana and still can't tell

## eigen visual analysis
(eis <- eisera(tus, tav2, scannf = FALSE))
barplot(eis$eig)
scatter(eis)

### THIRD ORDER
pcc <- mcp(locs)
image(maps)
plot(pcc, col = rainbow(6), add = TRUE)

hr <- do.call('data.frame', lapply(1:nrow(pcc), function(i) {
              over(maps, geometry(pcc[i,]))
}))
names(hr) <- slot(pcc, 'data')[,1]
coordinates(hr) <- coordinates(maps)
gridded(hr) <- TRUE
mimage(hr)
