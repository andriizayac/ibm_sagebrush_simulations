sage_sim<-function(competition_strength,distance_strength,survrate=1, timesteps=59) {
  sizemat1<-matrix(NA,nrow=nrow(test),ncol=timesteps)
  sizemat1[,1]<-test$height1_cm
  dist.fx<-function(x1,x2,y1,y2) {sqrt((x1-x2)^2+(y1-y2)^2)}
  rm.zero <- function(x){x<- x[x != 0]}
  for(i in 2:timesteps) {
    
    for(j in 1:nrow(test)) {
      
      dist_vec<-dist.fx(test$x[j],test$x,test$y[j],test$y)
      #comp_term<- sum(((exp(-distance_strength*dist_vec)/sizemat1[i-1])*competition_strength),na.rm=TRUE)-1 #size-dependent competition
      #comp_term<- sum((((-(dist_vec*(1-exp(-0.1*dist_vec))^(1-dist_vec*distance_strength)))/sizemat1[i-1])*competition_strength),na.rm=TRUE) #size-dependent competition
      comp_term<- sum(exp(-(rm.zero(dist_vec)/sizemat1[i-1])*competition_strength),na.rm=TRUE) #original
      comp_term <- abs(comp_term)
      if (i <= 18 & sizemat1[j, i-1] < 15 & !is.na(sizemat1[j,i-1])){
        sizemat1[j,i]<- (sizemat1[j,i-1] + rnorm(1,mean=1 + 0.05*sizemat1[j,i-1]+-0.013*comp_term, sd=exp(-0.015*sizemat1[j,i-1]))) * rbinom(1, 1, survrate)   
        if (sizemat1[j,i] < 2  | is.na(sizemat1[j,i])){ #& sizemat1[j,i] != 0
          sizemat1[j,i] <- NA
        }
      } else {
        sizemat1[j,i]<- sizemat1[j,i-1] + rnorm(1,mean=1 + 0.05*sizemat1[j,i-1]+-0.013*comp_term, sd=exp(-0.01*sizemat1[j,i-1]))
        if (sizemat1[j,i] < 2 | is.na(sizemat1[j,i])){
          sizemat1[j,i] <- NA
        }
      }
      sizemat1[j,i] <- sizemat1[j,i]/(1+0.0002*sizemat1[j,i])
    } 
  }
  
  matplot(t(sizemat1),type="l", xlim = c(1,60), ylim = c(1,200), xlab = "time (months)", ylab = "size",  main = "ABM simulated growth curves")
  plot(orchard1$x, orchard1$y, cex = 0.015*sizemat1[,timesteps])
  
  #return(sizemat1)
  bio.dists <- as.matrix(dist(cbind(test$x, test$y)))
  bio.dists.inv <- 1/bio.dists
  diag(bio.dists.inv) <- 0
  #bio.dists.inv[1:15, 1:15]
  moransI <- round(unlist(Moran.I(sizemat1[,timesteps], bio.dists.inv, na.rm = TRUE)[1]), digits = 3)
  
  cover <- round(sum(sizemat1[,timesteps], na.rm = TRUE), digits = 3)
  return(c(moransI, cover))
}

comp<-seq(1.5, 0.5, l=30)
moranIout<-rep(NA,times=length(comp))
#coverout<-rep(NA,times=length(comp))


for(k in 1:length(comp)) {
  moranIout[k]=sage_sim(competition_strength=comp[k], distance_strength = 1, timesteps = 59, survrate=1)
}

plot(moranIout~comp, ylab = "Moran's I", xlab = "Competition gradient")
