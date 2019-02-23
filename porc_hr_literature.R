## Porcupine home ranges from the literature
## (for manuscript Table 4)

## converting SD to SE 

## Coltrane and Sinnott 2012
  ## all winters (females): mean = 0.89 km2, SD = 0.26, n = 19
  0.26 / sqrt(19)
  ## all winters (males): mean = 1.11 km2, SD = 0.38, n = 12
  0.38 / sqrt(12)
  ## all winters, combined: mean = 0.98 km2, SD = 0.33, n = 31
  0.33 / sqrt(31)

## Roze 1987
  ## nonwinter (1982-84): mean = 0.649 km2, SD = 0.185, n = 14? *the sample size is animal-summers, not animals (but I guess you'd avg for each animal anyway)
  0.185 / sqrt(14)
    ## winter (1981/82 - 1983/84): mean = 0.074 km2, SD = 0.06, n = 9?
  0.06 / sqrt(9)
  # winter (1984/85): mean = 0.599, SD = 0.94, n = 7?
  0.94 / sqrt(7)
  
win82_84 <- c(0.111, 0.08, 0.034, 0.008, 0.036, 0.017, 0.173, 0.052, 0.152)
  mean(win82_84)
  a <- sd(win82_84) / sqrt(length(win82_84))

win84_85 <- c(0.094, 0.643, 0.387, 0.009, 2.665, 0.024, 0.37)
  mean(win84_85)
  b <- sd(win84_85) / sqrt(length(win84_85))

win_both <- c(win82_84, win84_85)
  mean(win_both)
  sd(win_both) / sqrt(length(win_both))

non_win <- c(0.353, 0.784, 0.638, 0.55, 0.67, 0.606, 0.327, 0.713, 0.606, 0.775, 1.095, 0.6, 0.656, 0.701)
  mean(non_win)
  sd(non_win) / sqrt(length(non_win))
  
## Morin et al. 2005 (summer only)
  ## males: mean = 0.209 km2, SE = 0.058, n = 8
  0.058 * sqrt(8) #SD
  
  ## females: mean = 0.154 km2, SE = 0.056, n = 9
  
## Smith 1979 (Oregon)
smith_hr <- c(3.9, 8.1, 28.8, 80.8, 82.1, 22.6)
mean(smith_hr)

## Craig and Keller 1986 (Idaho) - minimum area method
## SUMMER: (7 females; no males had enough to calculate HR)
craig_hr_s_min <- c(23.9, 19.6, 43.8, 27.9, 28.9, 9.1, 8.3) ## min area method
  mean(craig_hr_s_min) #23.07
  sd(craig_hr_s_min) / sqrt(7) #4.66 (they report SD (12.4))

craig_hr_s_mod <- c(5.4, 9.4, 4.7, 8.3, 14.9, 3.3, 4.3)  ## modified min area method
  mean(craig_hr_s_mod) #7.19
  sd(craig_hr_s_mod) / sqrt(7) #1.53 (SD = 4.1)
  
  ## in km2:
  7.19 / 100 #0.07
  1.53 / 100 #0.02

  ## WINTER: (unknown sex)
craig_hr_w_min <- c(0.02, 0.18) ## min area method
  mean(craig_hr_w_min) #0.10
  sd(craig_hr_w_min) / sqrt(2) #0.08 (SD = 0.1)

craig_hr_w_mod <- c(0.02, 0.12) ## modified min area method
  mean(craig_hr_w_mod) #0.07
  sd(craig_hr_w_mod) / sqrt(2) #0.05 (SD = 0.07)

  ## in km2:
  0.07 / 100 #<0.01
  0.05 / 100 #<0.01
  