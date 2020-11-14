nt <- 500
cmt <- rep(1, nt)
evid <- rep(2, nt)
addl <- rep(0, nt)
ss <- rep(0, nt)
amt <- rep(0, nt)
time <- seq(0, 49.9, 0.1)
rate <- rep(0, nt)
ii <- rep(0, nt)
nm <- data.frame(cmt, evid, addl, ss, amt, time, rate, ii)
nm$evid[1] <- 1
nm$amt[1] <- 10000
nm$addl[1] <- 1
nm$ii[1] <- 24
nm$evid[181] <- 4
nm$amt[181] <- 8000
