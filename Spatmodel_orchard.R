library(ape)
library(raster)
library(ggplot2)
# model test dataset:
###df "otest" is a subset from the main df orchard1 (merged) created for explotative purposes ###
orchard1 <- read.csv("orchard1.csv", header = TRUE)
test <- orchard1[, c("x", "y", "height1_cm")]

# replace NA's in otest with zeros #
test[is.na(test)] <- 0
test<-test[!is.na(test$height1_cm),]
#create vectors x and y for the grid
#x<-c(1:30)
#y<-c(1:30)

# create data frame with a grid x by y.
#data<-expand.grid(x,y)     

#name columnes
#colnames(data)<-c("x","y")

#add third column to the x and y dataframe with the first size measurement (all of equal size)
#data$size<-rep(5,length=30*30)

# create vector with the number of timesteps
timesteps <- 150 # 59 months

# create empty matrix with nrow = length(data), and ncol = the naumber of timesteps.
sizemat1 <- matrix(NA,nrow=nrow(test),ncol=timesteps)

# assign the values of the first size measurement to the first timestep column of the matrix
sizemat1[,1] <- test$height1_cm




for(i in 2:timesteps) {
  
  for(j in 1:nrow(test)) {
  
sizemat1[j,i]<-sizemat1[j,i-1]+rnorm(1,mean=1+0.1*sizemat1[j,i-1],sd=0.05*sizemat1[j,i-1]) #rnorm(1,mean=1,sd=0.5)    
    
  }
}

matplot(t(sizemat1),type="l")

#spatial term

dist.fx<-function(x1,x2,y1,y2) {sqrt((x1-x2)^2+(y1-y2)^2)}
timesteps <- 59
sizemat1<-matrix(NA,nrow=nrow(test),ncol=timesteps)
#sizemat2<-matrix(NA,nrow=nrow(test),ncol=timesteps)
sizemat1[,1]<-test$height1_cm

competition_strength=0.9
survrate <- 0.9

for(i in 2:timesteps) {
  
  for(j in 1:nrow(test)) {
    
    dist_vec<-dist.fx(test$x[j],test$x,test$y[j],test$y)
    comp_term<- sum(exp(-(dist_vec/sizemat1[i-1])*competition_strength),na.rm=TRUE)-1 #size-dependent competition
   
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

matplot(t(sizemat1),type="l", xlab = "time (months)", ylab = "size", main = "ABM simulated growth curves")
matplot(t(sizemat1),type="l", xlim = c(1,200), ylim = c(1,300), xlab = "time (months)", ylab = "size",  main = "ABM simulated growth curves")
timesteps <- 59
sizemat1<-matrix(NA,nrow=nrow(test),ncol=timesteps)
sizemat1[,1]<-test$height1_cm

sage_sim<-function(competition_strength,distance_strength,survrate = 1, timesteps = 59) {
  sizemat1<-matrix(NA,nrow=nrow(test),ncol=timesteps)
  sizemat1[,1]<-test$height1_cm
  dist.fx<-function(x1,x2,y1,y2) {sqrt((x1-x2)^2+(y1-y2)^2)}
  for(i in 2:timesteps) {
    
    for(j in 1:nrow(test)) {
      
      dist_vec<-dist.fx(test$x[j],test$x,test$y[j],test$y)
      comp_term<- sum(distance_strength*exp(-((dist_vec)/sizemat1[i-1])*competition_strength),na.rm=TRUE)-1 #size-dependent competition
      
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
  plot(orchard1$x, orchard1$y, cex = 0.02*sizemat1[,59])
  
  #matplot(t(sizemat1),type="l", xlim = c(1,60), ylim = c(1,200), xlab = "time (months)", ylab = "size",  main = "ABM simulated growth curves")
  
#return(sizemat1)
  bio.dists <- as.matrix(dist(cbind(test$x, test$y)))
  bio.dists.inv <- 1/bio.dists
  diag(bio.dists.inv) <- 0
  #bio.dists.inv[1:15, 1:15]
  moransI <- round(unlist(Moran.I(sizemat1[,59], bio.dists.inv, na.rm = TRUE)[1]), digits = 3)

  cover <- round(sum(sizemat1[,59], na.rm = TRUE), digits = 3)
  return(c(moransI, cover))
}


#comp<-c(0.7,0.8,0.9)
comp<-seq(0.05, 0.2, by = 0.05)
moranIout<-rep(NA,times=length(comp))
coverout<-rep(NA,times=length(comp))

for(k in 1:length(comp)) {
moranIout[k]=sage_sim(competition_strength=comp[k], distance_strength = 0.21)
}

plot(moranIout~comp)

for(z in 1:length(comp)){
  coverout[z] <- sage_sim(competition_strength = comp[z], distance_strength = 0.21)[2]
}

plot(coverout ~ moranIout)

#metrics that can be compared to reality
#hist(sizemat1[,59],breaks=20)
#hist(orchard1$height_cm,breaks=20)
plot(density(orchard1$height_cm,na.rm=T),ylim=c(0,0.05), xlab = "size",
     main = "Real vs Simulated size distributions")
lines(density(sizemat1[,59],na.rm=T),col="red",add=T)
#lines(density(rnorm(1000, mean = 90, sd = 30)), col = "green", add = TRUE)



#proportional circles

#for next time:
#1. making the proportion plots to look at data in space
#2. think about how to put survival in here as well

df <- data.frame(orchard1$x[1:468], orchard1$y[1:468], sizemat1[1:468,59])
names(df) <- c("x", "y", "sizemat1")
plot <- ggplot(df, aes(x, y, size = sizemat1, col = "red")) + geom_point()+ ggtitle("AGB common garden simulation")

x <- seq(1, 36, by = 0.5)
plot(x, 4/(exp((x)^2)))

plot(exp(-(x/5)*5))
plot(exp((-x)*3.5))

plot(3*exp(-((x)/sizemat1[i-1])*0.2),na.rm=TRUE)
