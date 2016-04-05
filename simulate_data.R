willows <- 0.9
#grass <- 0.05
#others <- 1-(willows+grass)
#sum(observations==1)/length(observations)

#my.obs <- rmultinom(100, 1, prob=all.probs)

all.probs <- c(willows, grass, others)
observations <- rbinom(100, 1, willows)
pres <- rep(1, length(observations))

abs.observations <- c(rep(1, 50), rep(0, 50))
abs <- rep(0, length(abs.observations))

veg.classes <- c(observations, abs.observations)
observations <- c(pres, abs)

my.table <- data.frame(observations, veg.classes)

my.glm <- glm(observations ~ veg.classes, data=my.table, family="binomial")
summary(my.glm)
