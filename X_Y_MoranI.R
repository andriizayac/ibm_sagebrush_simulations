orchard1<- read.csv("orchard1.csv")

xMoranI <- rep(NA, times = 28)
yMoranI <- rep(NA, times = 17)

yseq <- seq(0,24, by = 1.5)

bio.dists.x <- 1/as.matrix(dist(cbind(rep(1,17), 1:28))) #used for y
bio.dists.y <- 1/as.matrix(dist(cbind(rep(1,28), 1:28))) #used for x


for(i in yseq){
  a <- orchard1[orchard1$y == i,]
 yMoranI[i] <- unlist(Moran.I(a$volumelast, bio.dists.x, na.rm = TRUE)[1])
}
for(i in c(1:28)){
 a <- orchard1[orchard1$x == i,]
   xMoranI[i] <- unlist(Moran.I(a$volumelast, bio.dists.y, na.rm = TRUE)[1])
}

moraxes <- data.frame("x" = xMoranI, "y" = append(yMoranI, rep(NA, times = 4)))

b<- moraxes[which(moraxes$y != "NA"),]
aa<-ggplot(moraxes, aes(1:28, x))+geom_point()+geom_line(col = "red")
aa+ylab("Moran's I")+xlab("Rows on x")+ylim(-0.45, 0.25)
bb<- ggplot(b, aes(1:16, y))+geom_point()+geom_line(col = "blue")
bb+ylab("Moran's I")+xlab("Rows on y")+ylim(-0.45, 0.25)

