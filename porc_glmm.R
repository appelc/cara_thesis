#############################################################
### Taking a stab at fitting a GLMM in WinBugs with porc data
### based on notes from STAT 510 lab (Mar 28)
#############################################################

library(R2WinBUGS)
bugs.dir <- "C:/Programs/WinBugs/winbugs14"

## the lab does have code for simulating binomial data
## (would need to do multinomial, actually)






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
