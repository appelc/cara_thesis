## GEOMETRIC MEAN AND CONFIDENCE INTERVAL CALCULATION

x <- c(2580,1279,490,1734,511)

# 1. take log of data
logx <- log(x)

# 2. take arithmetic mean of the logged data
meanlogx <- mean(logx)

# 3. exponentiate that mean
xgm <- exp(meanlogx) # this is the geometric mean
geometric.mean(x) # should give the same answer

# two other ways of doing that:
xgm2 <- prod(x)^(1/length(x))
exp(mean(log(x)))

# 4. now take 95% confidence intervals of the logged data ('logx')
# first find SE
selogx <- sd(logx) / sqrt(length(logx))
# then add/subtract 1.95*SE from the arithmetic mean of the logged data
lower <- meanlogx - 1.96*selogx
upper <- meanlogx + 1.96*selogx

# 5. exponentiate the upper/lower CI bounds to get the 95% CI of the geometric mean
exp(lower)
exp(upper)

ci.gm(x) #should give the same answer

# -----------------
# Write a function to do this
# (*should* be the same as 'ci.gm')

geometric.ci <- function(x) {
    gm = (prod(x))^(1/(length(x)))
    logx = log(x)
    logx.se = sd(logx) / sqrt(length(logx))
    cil = exp(mean(logx) - 1.96*logx.se)
    ciu = exp(mean(logx) + 1.96*logx.se)
    vec = c(round(cil,4), round(ciu,4))
    return (vec)
}
