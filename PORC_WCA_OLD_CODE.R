####################################################################################
####################################################################################
## OLD CODE - removed from 'PORC_WCA.R'
####################################################################################

############################################
### 7. step III (old: only one reference category)
## Create a matrix like Erickson et al. 2001 by subtracting the reference category log-ratio 
## from each each habitat/individual log-ratio (Erickson et al. call this 'd')
############################################
d_matrix_2 <- list()
d_means_2 <- list()
ref <- 4      ## column of 'conifer forest' (can change here for desired reference category)
for (j in 1:3){
  matrix_j <- log_ratios_2[[j]]
  ref_matrix <- matrix(rep(matrix_j[,ref], ncol(matrix_j)), nrow = nrow(matrix_j), byrow = FALSE)
  colnames(ref_matrix) <- colnames(matrix_j)
  d_matrix_2[[j]] <- matrix_j - ref_matrix ## store for doing t-tests
  
  d_means_j <- colMeans(d_matrix_2[[j]]) ## shouldn't be any NAs in 2nd order (no need for rm.na = TRUE)
  d_means_j <- stack(d_means_j)
  colnames(d_means_j) <- c('d', 'veg')
  d_means_j$rank <- rank(-d_means_j$d) ## negative sign so it ranks largest -> smallest
  d_means_2[[j]] <- d_means_j ##store
}

##2nd: II) Test for overall selection (use different from random) using Wilk's lambda / MANOVA
## ** Haven't been able to get these to work
wilks_results_2 <- list() ## store Wilk's lambda results by season
for (j in 1:3){
  use_avail_j <- use_avail_2[[j]]
  groups <- as.factor(use_avail_j$veg)
  x <- as.matrix(use_avail_j[,3:4])
  wilks_j <- Wilks.test(x, grouping = groups, method = 'c') ## which method?
  wilks_results_2[[j]] <- wilks_j
}
## lambda is 0 for overall and summer... that doesn't seem right

## try (M)ANOVA for comparison
anova_2 <- aov(log_ratio ~ veg*id, data = use_avail_2[[1]])
summary(anova_2)
plot(anova_2)

## how is it a MANOVA? what are the 2 dependent variables? use & avail?
## ** revisit this...
Y <- data.frame(use_avail_2[[1]][,3:4])
A <- use_avail_2[[1]][,2] #veg
B <- use_avail_2[[1]][,1] #id
manova_2 <- manova(Y ~ A*B)
manova_2 <- manova(cbind(used_prop, avail_prop) ~ veg*id, data = use_avail_2[[1]])
summary(manova_2, test = 'Wilks')

## OR ...

## 2nd: III) Test for overall selection using randomization simulations, recommended by Pendleton et al. 1998 and Thomas and Taylor 2006
##    Help here: http://media.pearsoncmg.com/aw/aw_kuiper_online_resources/R_Manual.pdf, and https://www.youtube.com/watch?v=ds8nbvHVu0s
##    Also check out code used in 'compana' for randomisation
reps <- 10000
results <- numeric(reps) # establish a vector of the right length, with '0’s initially
x <- c() # combines use & avail vectors
for (i in 1:reps) {
  temp <- sample(x)
  results[i] <- mean(temp[1:5])-mean(temp[6:10])
}
## *** REVISIT THIS. 'compana' does randomization automatically and generates a Lambda and p-value. ***


##########################################
## IV) OLD: t-tests on the differences in log-ratio between each veg category and the reference category,
##      then pairwise between all veg categories (see Erickson et al. 2001, pg 228)
##  - I removed 1-sample t-tests against the reference category using t.test(season[,k], mu = 0) because
##    this is the exact same as doing a paired t-test with the reference category, because it is just a 
##    column of 0s. Maybe double-check that this is theoretically correct but I get the exact same values
##    of t and p from the above test and t.tes(season[,k], season[,ref], paired = TRUE)
ttests_2 <- list()
ref <- 4
for (j in 1:3){
  season <- d_matrix_2[[j]]
  pairwise_k <- NULL
  for (k in 1:length(season)){
    pairwise_l  <- NULL
    for (l in 1:length(season)){
      ttest_l <- t.test(season[,k], season[,l], paired = TRUE)
      ttest_l_df <- data.frame(k, l, ttest_l$estimate, ttest_l$conf.int[1], ttest_l$conf.int[2], ttest_l$p.value)
      colnames(ttest_l_df) <- c('v1', 'v2', 'mean_diff', 'lci_95', 'uci_95', 'p')
      pairwise_l <- bind_rows(pairwise_l, ttest_l_df) ## store
    }
    pairwise_k <- bind_rows(pairwise_k, pairwise_l)
  }
  veg_key <- data.frame('veg' = (names(season)), 'veg_id' = 1:length(season))
  ref_key <- data.frame('veg' = 'conifer forest', 'veg_id' = 0) ## can modify reference category label
  veg_key <- rbind(veg_key, ref_key)
  pairwise_k$veg1 <- veg_key[match(pairwise_k$v1, veg_key$veg_id), 'veg'] 
  pairwise_k$veg2 <- veg_key[match(pairwise_k$v2, veg_key$veg_id), 'veg']
  ttests_2[[j]] <- pairwise_k
}

##3rd: I) Compute the difference between log-transformed 'use_3' and log-transormed 'avail_2' data
use_avail_3 <- list() ## combine use and availability data
for (j in 1:3){
  use.j <- data.frame(use_3[[j]])
  use.j$id <- rownames(use.j)
  use <- reshape(use.j, varying = 1:9, direction = 'long', v.names = 'used_prop', timevar = 'veg',
                 idvar = 'id', times = colnames(use.j[,1:9]))
  avail.j <- data.frame(avail_3[[j]])
  avail.j$id <- rownames(avail.j) ## we have some 0s. leave for now and replace at log-ratio step
  avail <- reshape(avail.j, varying = 1:9, direction = 'long', v.names = 'avail_prop', timevar = 'veg',
                   idvar = 'id', times = colnames(avail.j[,1:9]))
  use_avail <- use
  use_avail$avail_prop <- avail$avail_prop
  rownames(use_avail) <- NULL
  ## What to do with 0 use values?
  use_avail$used_prop[use_avail$avail_prop != 0 & use_avail$used_prop == 0] <- 1e-10 #it really IS no use, but 0 will throw off log-ratios
  use_avail$used_prop[use_avail$avail_prop == 0 & use_avail$used_prop == 0] <- NA #missing data; will replace with mean down below
  ## Compute log-ratios: either log(use)-log(avail) OR log(use/avail)
  use_avail$log_ratio <- log(use_avail$used_prop) - log(use_avail$avail_prop)
  use_avail_3[[j]] <- use_avail ## STORE     
}

log_ratios_3 <- list() ## rearrange log-ratios to matrix shape for steps below
for (j in 1:3){
  use_avail_3[[j]]$log_ratio <- log(use_avail_3[[j]]$used_prop) - log(use_avail_3[[j]]$avail_prop)
  log_ratios_j <- use_avail_3[[j]][,c(1:2, 5)]
  log_ratios_j <- reshape(log_ratios_j, timevar = 'veg', idvar = 'id', direction = 'wide')
  names(log_ratios_j) <- gsub('log_ratio.', '', names(log_ratios_j)) ## get rid of 'sel' in column names
  names(log_ratios_j) <- gsub('[.]', ' ', names(log_ratios_j))
  rownames(log_ratios_j) <- log_ratios_j[,1]
  log_ratios_j <- log_ratios_j[,-1]
  for (k in 1:ncol(log_ratios_j)){
    ## replace missing values (avail & use = NA) with the mean of log-ratios for each veg type (column) based on all non-NA values (ala Aebischer et al. 1993, Appendix 2)
    log_ratios_j[is.na(log_ratios_j[,k]), k] <- mean(log_ratios_j[,k], na.rm = TRUE)
  }
  log_ratios_3[[j]] <- log_ratios_j
}
## This keeps the column means (mean log-ratio per veg type) the same as they were with just the non-NA values
## But see caveats in Aebischer et al. 1993 (Appendix 2) re: independence and suggestion for computing mean lambda

##3rd: II) Test for overall selection (use different from random) using Wilk's lambda / MANOVA
## ** Haven't been able to get these to work
wilks_results_3 <- list() ## store Wilk's lambda results by season
for (j in 1:3){
  use_avail_j <- use_avail_3[[j]]
  groups <- as.factor(use_avail_j$veg)  ## veg classes
  x <- as.matrix(use_avail_j[,3:4])  ## used vs. available proportions
  wilks_j <- Wilks.test(x, grouping = groups, method = 'c')  ## which method?
  wilks_results_3[[j]] <- wilks_j
}
## good, all significantly different from random

##3rd: III) Create a ranking matrix based on pairwise differences between EACH category and each OTHER category, 
##      like Aebischer et al. 1993, then test the mean difference per veg type (across animals) for nonrandom
##      use using a one-sample t-test (H0: mu = 0). This is slightly different from Erickson et al. 2001; there shouldn't be just one reference category, right?

ttests_3 <- list()
ranks_3 <- list()

for (j in 1:3){ 
  matrix_j <- log_ratios_3[[j]]
  ttests_r <- NULL
  ranks_r <- NULL
  for (r in 1:9){ ## cycle thru each veg type (column) as the reference (subtracted from all others)
    ref_matrix <- matrix(rep(matrix_j[,r], ncol(matrix_j)), nrow = nrow(matrix_j), byrow = FALSE)
    diff_r <- matrix_j - ref_matrix ## calculate difference in log-ratio between each column and ref.
    ttests_rs <- NULL
    for (s in 1:9){ ## t-tests for the difference in log-ratio for each veg type (column) against 0 (random use)
      ttest_s <- t.test(diff_r[,s], mu =  0) 
      ttest_s_df <- data.frame(r, s, ttest_s$estimate, ttest_s$conf.int[1], ttest_s$conf.int[2], ttest_s$p.value)
      colnames(ttest_s_df) <- c('r', 's', 'mean_diff', 'lci_95', 'uci_95', 'p')
      ttests_rs <- bind_rows(ttests_rs, ttest_s_df) ## store
    }
    pos_r <- nrow(ttests_rs[ttests_rs$mean_diff > 0,]) ## how many positive means for each 'r' (veg type as reference)?
    ranks_r <- rbind(ranks_r, data.frame(r, pos_r))
    ttests_r <- rbind(ttests_r, ttests_rs)
  }
  veg_key <- data.frame('veg' = (names(matrix_j)), 'veg_id' = 1:length(matrix_j))
  ttests_r$veg1 <- veg_key[match(ttests_r$r, veg_key$veg_id), 'veg'] # match veg names
  ttests_r$veg2 <- veg_key[match(ttests_r$s, veg_key$veg_id), 'veg'] # match veg names
  ttests_3[[j]] <- ttests_r 
  ranks_3[[j]] <- data.frame('veg' = colnames(matrix_j), 'rank' = ranks_r$pos_r) # 0 is the most selected and 8 is the least
}

## Compare 2nd- versus 3rd-order selection ranks:
ranks_2
ranks_3

##3rd: V) Create table of significance codes based on t-tests (see Beasley et al. 2007 for example)
ttest_sig_3 <- list()
for (j in 1:3){
  tests_j <- ttests_3[[j]][,c(3, 6:8)]
  tests_j$sig[tests_j$mean_diff > 0] <- '+'
  tests_j$sig[tests_j$mean_diff > 0 & tests_j$p <= 0.05 & tests_j$p > 0.01] <- '++'
  tests_j$sig[tests_j$mean_diff > 0 & tests_j$p <= 0.01 & tests_j$p > 0.001] <- '+++'
  tests_j$sig[tests_j$mean_diff > 0 & tests_j$p <= 0.001] <- '++++'
  tests_j$sig[tests_j$mean_diff < 0] <- '\u2013' ## en dashes!
  tests_j$sig[tests_j$mean_diff < 0 & tests_j$p <= 0.05 & tests_j$p > 0.01] <- '\u2013 \u2013'
  tests_j$sig[tests_j$mean_diff < 0 & tests_j$p <= 0.01 & tests_j$p > 0.001] <- '\u2013 \u2013 \u2013'
  tests_j$sig[tests_j$mean_diff < 0 & tests_j$p <= 0.001] <- '\u2013 \u2013 \u2013 \u2013'
  matrix_j <- cast(tests_j, veg1 ~ veg2, value = 'sig')
  ttest_sig_3[[j]] <- matrix_j
}

############################################
## 8. Boxplots (difference in log ratio between each veg category and the reference category)
## These are similar to Figure 5 in Millspaugh et al. 2006
## - do I want horizontal or vertical? add '+ coord_flip()' 
## - 'guide = FALSE' in scale_fill_manual turns off the legend (since I have axis labels)
## - plot all on the same scale (i.e., same 'ylim' min and max values)? *
##    * Make sure 'ylim' includes all min and max values! Otherwise it will mess up some of the boxplots 
##      (e.g., it plotted fruit tree mean of -8.7 as positive 3)! Use min(summer_melt_2), etc., to check.
############################################

############################################
##    a. 2nd order
############################################
overall_melt_2 <- melt(d_matrix_2[[1]]) ## message OK
colnames(overall_melt_2) <- c('veg', 'd')
overall_melt_2$rank <- d_means_2[[1]][match(overall_melt_2$veg, d_means_2[[1]]$veg), 'rank']

summer_melt_2 <- melt(d_matrix_2[[2]]) ## message OK
colnames(summer_melt_2) <- c('veg', 'd')
summer_melt_2$rank <- d_means_2[[2]][match(summer_melt_2$veg, d_means_2[[2]]$veg), 'rank'] 

winter_melt_2 <- melt(d_matrix_2[[3]])
colnames(winter_melt_2) <- c('veg', 'd')
winter_melt_2$rank <- d_means_2[[3]][match(winter_melt_2$veg, d_means_2[[3]]$veg), 'rank']

s2 <- ggplot(data = summer_melt_2, aes(x = reorder(veg, rank), y = d, fill = as.factor(veg))) +
  stat_summary(fun.data = min.mean.se.max, geom = 'boxplot') +
  geom_point(position = position_dodge(0.8), size = 2) +
  scale_fill_manual(values = colors$veg_colors, guide = FALSE) +
  geom_hline(yintercept = 0, linetype = 'dashed') + ylim(-24, 10) + ## *see note above
  xlab('Vegetation Type') + ylab('Differences in Log Ratio') +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 35, hjust = 1),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 12, colour = 'black'),
        axis.line.x = element_line(size = 0.5, colour = 'black'),
        axis.line.y = element_line(size = 0.5, colour = 'black'),
        panel.background = element_rect(fill = 'white')) +
  geom_text(label ='*', aes(x = 2, y = -24), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 6, y = -24), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 7, y = -24), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 8, y = -24), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 9, y = -24), size = 8, colour = 'grey50')
s2 

w2 <- ggplot(data = winter_melt_2, aes(x = reorder(veg, rank), y = d, fill = as.factor(veg))) +
  stat_summary(fun.data = min.mean.se.max, geom = 'boxplot') +
  geom_point(position = position_dodge(0.8), size = 2) +
  scale_fill_manual(values = colors$veg_color, guide = FALSE) +
  geom_hline(yintercept = 0, linetype = 'dashed') + ylim(-24, 10) + ## *see note above
  xlab('Vegetation Type') + ylab('Differences in Log Ratio') +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 35, hjust = 1),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 12, colour = 'black'),
        axis.line.x = element_line(size = 0.5, colour = 'black'),
        axis.line.y = element_line(size = 0.5, colour = 'black'),
        panel.background = element_rect(fill = 'white')) +
  geom_text(label ='*', aes(x = 7, y = -24), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 8, y = -24), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 9, y = -24), size = 8, colour = 'grey50')
w2

## overall (at the most basic level, where did they choose their home ranges?)
o2 <- ggplot(data = overall_melt_2, aes(x = reorder(veg, rank), y = d, fill = as.factor(veg))) +
  stat_summary(fun.data = min.mean.se.max, geom = 'boxplot') +
  geom_point(position = position_dodge(0.8), size = 2) +
  scale_fill_manual(values = colors$veg_color, guide = FALSE) +
  geom_hline(yintercept = 0, linetype = 'dashed') + ylim(-24, 10) + ## *see note above
  xlab('Vegetation Type') + ylab('Differences in Log Ratio') +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 35, hjust = 1),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 12, colour = 'black'),
        axis.line.x = element_line(size = 0.5, colour = 'black'),
        axis.line.y = element_line(size = 0.5, colour = 'black'),
        panel.background = element_rect(fill = 'white')) +
  geom_text(label ='*', aes(x = 6, y = -24), size = 8, colour = 'grey50') +      
  geom_text(label ='*', aes(x = 7, y = -24), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 8, y = -24), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 9, y = -24), size = 8, colour = 'grey50')
o2

############################################
##    b. 3rd order
############################################
summer_melt_3 <- melt(d_matrix_3[[2]]) ## message OK
colnames(summer_melt_3) <- c('veg', 'd')
summer_melt_3$rank <- d_means_3[[2]][match(summer_melt_3$veg, d_means_3[[2]]$veg), 'rank'] 

winter_melt_3 <- melt(d_matrix_3[[3]])
colnames(winter_melt_3) <- c('veg', 'd')
winter_melt_3$rank <- d_means_3[[3]][match(winter_melt_3$veg, d_means_3[[3]]$veg), 'rank']

s3 <- ggplot(data = summer_melt_3, aes(x = reorder(veg, rank), y = d, fill = as.factor(veg))) +
  stat_summary(fun.data = min.mean.se.max, geom = 'boxplot') +
  geom_point(position = position_dodge(0.8), size = 2) +
  scale_fill_manual(values = colors$veg_colors, guide = FALSE) +
  geom_hline(yintercept = 0, linetype = 'dashed') + ylim(-5, 5) + ## *see note above
  xlab('Vegetation Type') + ylab('Differences in Log Ratio') +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 35, hjust = 1),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 12, colour = 'black'),
        axis.line.x = element_line(size = 0.5, colour = 'black'),
        axis.line.y = element_line(size = 0.5, colour = 'black'),
        panel.background = element_rect(fill = 'white')) +
  geom_text(label ='*', aes(x = 1, y = -5), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 2, y = -5), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 3, y = -5), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 8, y = -5), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 9, y = -5), size = 8, colour = 'grey50')
s3

w3 <- ggplot(data = winter_melt_3, aes(x = reorder(veg, rank), y = d, fill = as.factor(veg))) +
  stat_summary(fun.data = min.mean.se.max, geom = 'boxplot') +
  geom_point(position = position_dodge(0.8), size = 2) +
  scale_fill_manual(values = colors$veg_color, guide = FALSE) +
  geom_hline(yintercept = 0, linetype = 'dashed') + ylim(-20, 5) + ## *see note above
  xlab('Vegetation Type') + ylab('Differences in Log Ratio') +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 35, hjust = 1),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 12, colour = 'black'),
        axis.line.x = element_line(size = 0.5, colour = 'black'),
        axis.line.y = element_line(size = 0.5, colour = 'black'),
        panel.background = element_rect(fill = 'white')) +
  geom_text(label ='*', aes(x = 1, y = -20), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 8, y = -20), size = 8, colour = 'grey50') +
  geom_text(label ='*', aes(x = 9, y = -20), size = 8, colour = 'grey50')
w3 

# ----------------------------------------------------------------
# ** Try just # points / total for use for fun and to compare with compana
# ----------------------------------------------------------------
ids <- unique(sum.locs$id)
summer.prop <- list()
summer.matrix <- NULL
for (i in ids){
  sum.i <- sum.locs[sum.locs$id == i,]
  sum.sp <- SpatialPointsDataFrame(data.frame(sum.i$utm_e, sum.i$utm_n),
                                   data = data.frame(sum.i),
                                   proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
  sum.sp@data$veg <- over(sum.sp, veg)$Class_4
  sum.over <- data.frame(sum.sp@data[,c(2, 7:9)], 'season' = rep('sum', nrow(sum.sp@data)))
  sum.over.df <- data.frame('id' = i, 'season' = rep('sum', 9), stack(table(sum.over$veg)))
  sum.over.df$prop <- sum.over.df$values / sum(sum.over.df$values)
  colnames(sum.over.df) <- c('id', 'season', 'points', 'veg', 'prop')
  summer.prop[[i]] <- sum.over.df
  summer.matrix <- rbind(summer.matrix, sum.over.df$prop)
}
ids <- unique(win.locs$id)
winter.prop <- list()
winter.matrix <- NULL
for (i in ids){
  win.i <- win.locs[win.locs$id == i,]    
  win.sp <- SpatialPointsDataFrame(data.frame(win.i$utm_e, win.i$utm_n),
                                   data = data.frame(win.i),
                                   proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
  win.sp@data$veg <- over(win.sp, veg)$Class_4
  win.over <- data.frame(win.sp@data[,c(2, 7:9)], 'season' = rep('win', nrow(win.sp@data)))
  win.over.df <- data.frame('id' = i, 'season' = rep('win', 9), stack(table(win.over$veg)))
  win.over.df$prop <- win.over.df$values / sum(win.over.df$values)
  colnames(win.over.df) <- c('id', 'season', 'points', 'veg', 'prop')
  winter.prop[[i]] <- win.over.df
  winter.matrix <- rbind(winter.matrix, win.over.df$prop)
}

rownames(summer.matrix) <- unique(sum.locs$id)
colnames(summer.matrix) <- unique(summer.prop[[1]]$veg)
rownames(winter.matrix) <- unique(win.locs$id)
colnames(winter.matrix) <- unique(winter.prop[[1]]$veg)

compana_summer <- compana(summer.matrix, avail_3[[2]], test = 'randomisation', rnv = 0.000001, nrep = 10000, alpha = 0.05)
(eis <- eisera(summer.matrix, avail_3[[2]], scannf = FALSE))
barplot(eis$eig) ## what does this tell us?
scatter(eis)

compana_winter <- compana(winter.matrix, avail_3[[3]], test = 'randomisation', rnv = 0.000001, nrep = 10000, alpha = 0.05)
(eis <- eisera(winter.matrix, avail_3[[3]], scannf = FALSE))
barplot(eis$eig) ## what does this tell us?
scatter(eis)

## PLOTS

# ----------------------------------------------------------------
## try just bar plots of sel_means with CI
# ----------------------------------------------------------------
## 2nd order:
sel_means_2[[2]]$rank <- ranks_2[[2]][match(rownames(sel_means_2[[2]]), ranks_2[[2]]$veg), 'pos_r']
sel_means_2[[2]] <- sel_means_2[[2]][order(sel_means_2[[2]]$rank),]
sel_means_2[[3]]$rank <- ranks_2[[3]][match(rownames(sel_means_2[[3]]), ranks_2[[3]]$veg), 'pos_r']
sel_means_2[[3]] <- sel_means_2[[3]][order(sel_means_2[[3]]$rank),]

limits <- aes(ymax = sel_means_2[[2]]$uci, ymin = sel_means_2[[2]]$lci)
sum2 <- ggplot(data = sel_means_2[[2]], aes(x = reorder(rownames(sel_means_2[[2]]), rank), y = geo.mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), aes(fill = as.factor(rownames(sel_means_2[[2]])))) +
  geom_errorbar(limits, position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = colors$veg_color, guide = FALSE) +
  # geom_point(data = sr_2s_no, aes(x = reorder(rownames(sr_2s_no), rank)), y = sel_ratio) +
  xlab('Vegetation type') + ylab('Selection Ratios (95% CI') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.line.x = element_line(size = 1, colour = 'black'),
        axis.line.y = element_line(size = 1, colour = 'black'),
        panel.background = element_rect(fill = 'white'),
        axis.ticks.x = element_line(size = 1, colour = 'black'))
sum2

limits <- aes(ymax = sel_means_2[[3]]$uci, ymin = sel_means_2[[3]]$lci)
win2 <- ggplot(data = sel_means_2[[3]], aes(x = reorder(rownames(sel_means_2[[3]]), rank), y = geo.mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), aes(fill = as.factor(rownames(sel_means_2[[3]])))) +
  geom_errorbar(limits, position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = colors$veg_color, guide = FALSE) +
  xlab('Vegetation type') + ylab('Selection Ratios (95% CI') + #ylim(0, 2) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.line.x = element_line(size = 1, colour = 'black'),
        axis.line.y = element_line(size = 1, colour = 'black'),
        panel.background = element_rect(fill = 'white'),
        axis.ticks.x = element_line(size = 1, colour = 'black'))
win2
# ----------------------------------------------------------------
## 3rd order:
sel_means_3[[2]]$rank <- ranks_3[[2]][match(rownames(sel_means_3[[2]]), ranks_3[[2]]$veg), 'pos_r']
sel_means_3[[2]] <- sel_means_3[[2]][order(sel_means_3[[2]]$rank),]
sel_means_3[[3]]$rank <- ranks_3[[3]][match(rownames(sel_means_3[[3]]), ranks_3[[3]]$veg), 'pos_r']
sel_means_3[[3]] <- sel_means_3[[3]][order(sel_means_3[[3]]$rank),]

limits <- aes(ymax = sel_means_3[[2]]$uci, ymin = sel_means_3[[2]]$lci)
sum3 <- ggplot(data = sel_means_3[[2]], aes(x = reorder(rownames(sel_means_3[[2]]), rank), y = geo.mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), aes(fill = as.factor(rownames(sel_means_3[[2]])))) +
  geom_errorbar(limits, position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = colors$veg_color, guide = FALSE) +
  xlab('Vegetation type') + ylab('Selection Ratios (95% CI') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.line.x = element_line(size = 1, colour = 'black'),
        axis.line.y = element_line(size = 1, colour = 'black'),
        panel.background = element_rect(fill = 'white'),
        axis.ticks.x = element_line(size = 1, colour = 'black'))
sum3

limits <- aes(ymax = sel_means_3[[3]]$uci, ymin = sel_means_3[[3]]$lci)
win3 <- ggplot(data = sel_means_3[[3]], aes(x = reorder(rownames(sel_means_3[[3]]), rank), y = geo.mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), aes(fill = as.factor(rownames(sel_means_3[[3]])))) +
  geom_errorbar(limits, position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = colors$veg_color, guide = FALSE) +
  xlab('Vegetation type') + ylab('Selection Ratios (95% CI') + ylim(0, 2) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.line.x = element_line(size = 1, colour = 'black'),
        axis.line.y = element_line(size = 1, colour = 'black'),
        panel.background = element_rect(fill = 'white'),
        axis.ticks.x = element_line(size = 1, colour = 'black'))
win3

