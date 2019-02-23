## median number of animals overlapped for males / females in summer / winter
## home ranges ('hr'): 95% KDE
## core areas ('ca'): 50% KDE

hr.s.f <- c(9, 8, 9, 2, 9, 9, 8, 7, 1, 4)
hr.s.m <- c(5, 5, 2, 2, 10, 9, 4, 5, 9)

median(hr.s.f)
mean(hr.s.f)
  sd(hr.s.f) / sqrt(length(hr.s.f))

median(hr.s.m)
mean(hr.s.m)
  sd(hr.s.m) / sqrt(length(hr.s.m))

hr.s <- c(hr.s.f, hr.s.m)
median(hr.s)

# --------------

hr.w.f <- c(2, 0, 0, 2, 0, 1)
hr.w.m <- c(2, 0, 3, 1, 1)

median(hr.w.f)
  mean(hr.w.f)
  sd(hr.w.f) / sqrt(length(hr.w.f))
  
median(hr.w.m)
  median(hr.w.m)
  sd(hr.w.m) / sqrt(length(hr.w.m))
  
# --------------
ca.s.f <- c(7, 3, 7, 0, 7, 4, 0, 3, 0, 2)
ca.s.m <- c(5, 2, 0, 2, 7, 5, 3, 1, 8)
  
median(ca.s.f)
  mean(ca.s.f)
  sd(ca.s.f) / sqrt(length(ca.s.f))
  
median(ca.s.m)
  median(ca.s.m)
  sd(ca.s.m) / sqrt(length(ca.s.m))
  