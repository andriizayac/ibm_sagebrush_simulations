bio.dists <- as.matrix(dist(cbind(test$x, test$y)))
bio.dists.inv <- 1/bio.dists
diag(bio.dists.inv) <- 0

yseq <- seq(0,24, by = 1.5)
bio.dists.y <- 1/as.matrix(dist(cbind(rep(1,17), yseq)))
diag(bio.dists.y) <- 0

bio.dists.x <- 1/as.matrix(dist(cbind(rep(1,17), 1:28)))
diag(bio.dists.x) <- 0

xMoranI <- rep(NA, times = 28)
yMoranI <- rep(NA, times = 17)

for(i in c(1:28)){
  a <- orchard1[orchard1$x == i,]
  xMoranI[i] <- unlist(Moran.I(a$volumelast, bio.dists.y, na.rm = TRUE)[1])
}


for(i in seq(0,24,1.5)){
  a <- orchard1[orchard1$y == 1,]
  yMoranI[i] <- unlist(Moran.I(a$volumelast, bio.dists.x, na.rm = TRUE)[1])
}


