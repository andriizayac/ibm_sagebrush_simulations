library(rstan)
library(ggplot2)
library(bayesplot)
library(tidybayes)
dist.fx=function(x1, x2, y1, y2) {sqrt((x2-x1)^2 + (y2-y1)^2)}
#rm.zero<-function(x){x<- x[x != 0]}

# import data
model1=readRDS("fit_aa1_segment_c.rds")
post1=extract(model1)
model2=readRDS("fit_aa2_segment_c.rds")
post2=extract(model2)
plot(model1,pars=c("a3"))
orchard1=read.csv("orchard1.csv")
mf=read.csv("MajorFlats_comp_df.csv")

# parameters
ID <- seq(1,length(x), by = 1)
N=400
t=40


distmat=matrix(NA,sqrt(N),sqrt(N))

# Interaction parameters census1
alpha1=mean(post1$a0)
alphat2n1=mean(post1$Itype[,2])
alphat4n1=mean(post1$Itype[,3])
alphav2n1=mean(post1$Itype[,4])
alphav4n1=mean(post1$Itype[,5])
alphaw4n1=mean(post1$Itype[,6])
beta1=mean(post1$c0)
betat2n1=mean(post1$c[,2])
betat4n1=mean(post1$c[,3])
betav2n1=mean(post1$c[,4])
betav4n1=mean(post1$c[,5])
betaw4n1=mean(post1$c[,6])
a31=mean(post1$a03)
a3t2n1=mean(post1$a3[,2])
a3t4n1=mean(post1$a3[,3])
a3v2n1=mean(post1$a3[,4])
a3v4n1=mean(post1$a3[,5])
a3w4n1=mean(post1$a3[,6])
a21=mean(post1$a2)
a11=mean(post1$a01)
a1t2n1=mean(post1$a1[,2])
a1t4n1=mean(post1$a1[,3])
a1v2n1=mean(post1$a1[,4])
a1v4n1=mean(post1$a1[,5])
a1w4n1=mean(post1$a1[,6])

# Interaction parameters census2
alpha2=mean(post2$a0)
alphat2n2=mean(post2$Itype[,1])
alphat4n2=mean(post2$Itype[,2])
alphav2n2=mean(post2$Itype[,3])
alphav4n2=mean(post2$Itype[,4])
alphaw4n2=mean(post2$Itype[,5])
beta2=mean(post2$c0)
betat2n2=mean(post2$c[,1])
betat4n2=mean(post2$c[,2])
betav2n2=mean(post2$c[,3])
betav4n2=mean(post2$c[,4])
betaw4n2=mean(post2$c[,5])
a32=mean(post2$a03)
a3t2n2=mean(post2$a3[,1])
a3t4n2=mean(post2$a3[,2])
a3v2n2=mean(post2$a3[,3])
a3v4n2=mean(post2$a3[,4])
a3w4n2=mean(post2$a3[,5])
a22=mean(post2$a2)
a12=mean(post2$a01)
a1t2n2=mean(post2$a1[,1])
a1t4n2=mean(post2$a1[,2])
a1v2n2=mean(post2$a1[,3])
a1v4n2=mean(post2$a1[,4])
a1w4n2=mean(post2$a1[,5])

spacing=seq(0.5,4,l=20)
total_cover=rep(NA,length(spacing))
for(i in 1:length(spacing)){
  xspace=spacing[i]
  yspace=spacing[i]
  x=rep(seq(1,l=20,by=xspace),times=20)
  y=sort(rep(seq(1,l=20,by=yspace),times=20))
  
  subspp=sample(rep(1:5,times=N/5),N)
  df<- data.frame(x,y, "vol" = 0.0001, ID = ID,subspp=subspp)
  df$spp=ifelse(df$subspp==1,"T2n",ifelse(df$subspp==2,"T4n",ifelse(df$subspp==3,"V2n",ifelse(df$subspp==4,"V4n",ifelse(df$subspp==5,"W4n",NA)))))
  
  total_cover[i]=sage_sim()[[2]]
}
plot(total_cover/max(total_cover)~spacing)

############################################# include density
sim=function(xspace=5,yspace=5,t=40){
  x=rep(seq(1,l=20,by=xspace),times=20)
  y=sort(rep(seq(1,l=20,by=yspace),times=20))
  
  ID=c(1:N)
  subspp=sample(rep(1:5,times=N/5),N)
  df<- data.frame(x,y, "vol" = 0.0001, ID = ID,subspp=subspp)
  df$spp=ifelse(df$subspp==1,"T2n",ifelse(df$subspp==2,"T4n",ifelse(df$subspp==3,"V2n",ifelse(df$subspp==4,"V4n",ifelse(df$subspp==5,"W4n",NA)))))
  a=sage_sim()
  return(sage_sim()[[2]])
  #print(df)
} # crap!!!
####################################################################### sage_sim
x=rep(seq(1,l=20,by=xspace),times=20)
y=sort(rep(seq(1,l=20,by=yspace),times=20))

subspp=sample(rep(1:5,times=N/5),N)
df<- data.frame(x,y, "vol" = 0.0001, ID = ID,subspp=subspp)
df$spp=ifelse(df$subspp==1,"T2n",ifelse(df$subspp==2,"T4n",ifelse(df$subspp==3,"V2n",ifelse(df$subspp==4,"V4n",ifelse(df$subspp==5,"W4n",NA)))))


sage_sim<-function(survintercept=1,survslope =1, timesteps=t) {
  
  sizemat=matrix(NA,N,timesteps)
  sizemat[,1]=0.0001
  xedge=(max(df$x)/sqrt(N))*3
  yedge=(max(df$y)/sqrt(N))*3
  edge<-df$ID[which(df$x < min(df$x)+xedge | df$x > max(df$x)-xedge | df$y < min(df$y)+yedge | df$y > max(df$y)-yedge)]
  
  for(i in 2:timesteps) {
    for(j in 1:N) {
      dist_vec=dist.fx(df$x[j],df$x,df$y[j],df$y)
      if(df$subspp[j] == 1){
        weight=plogis(2*(i-t/2))
        wa1t2n= weight*a1t2n1+(1-weight)*a1t2n2
        wa2= weight*a21+(1-weight)*a22
        walphat2n=weight*alphat2n1+(1-weight)*alphat2n2
        wbetat2n=weight*betat2n1+(1-weight)*betat2n2
        wa3t2n=weight*a3t2n1+(1-weight)*a3t2n2
      cf<- sum(wa1t2n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
      sizemat[j,i] <- sizemat[j,i-1]+
        rnorm(1, mean=walphat2n + sizemat[j,i-1]*wbetat2n + wa3t2n*cf, sd=0.001)
      
      } else if (df$subspp[j] ==2){
        weight=plogis(2*(i-t/2))
        wa1t4n= weight*a1t4n1+(1-weight)*a1t4n2
        wa2= weight*a21+(1-weight)*a22
        walphat4n=weight*alphat4n1+(1-weight)*alphat4n2
        wbetat4n=weight*betat4n1+(1-weight)*betat4n2
        wa3t4n=weight*a3t4n1+(1-weight)*a3t4n2
        cf<- sum(wa1t4n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
        sizemat[j,i] <- sizemat[j,i-1]+
          rnorm(1, mean=walphat4n + sizemat[j,i-1]*wbetat4n + wa3t4n*cf, sd=0.001)
        
      } else if (df$subspp[j] ==3){
        weight=plogis(2*(i-t/2))
        wa1v2n= weight*a1v2n1+(1-weight)*a1v2n2
        wa2= weight*a21+(1-weight)*a22
        walphav2n=weight*alphav2n1+(1-weight)*alphav2n2
        wbetav2n=weight*betav2n1+(1-weight)*betav2n2
        wa3v2n=weight*a3v2n1+(1-weight)*a3v2n2
        cf<- sum(wa1v2n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
        sizemat[j,i] <- sizemat[j,i-1]+
          rnorm(1, mean=walphav2n + sizemat[j,i-1]*wbetav2n + wa3v2n*cf, sd=0.001)
        
      } else if (df$subspp[j] ==4){
        weight=plogis(2*(i-t/2))
        wa1v4n= weight*a1v4n1+(1-weight)*a1v4n2
        wa2= weight*a21+(1-weight)*a22
        walphav4n=weight*alphav4n1+(1-weight)*alphav4n2
        wbetav4n=weight*betav4n1+(1-weight)*betav4n2
        wa3v4n=weight*a3v4n1+(1-weight)*a3v4n2
        cf<- sum(wa1v4n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
        sizemat[j,i] <- sizemat[j,i-1]+
          rnorm(1, mean=walphav4n + sizemat[j,i-1]*wbetav4n + wa3v4n*cf, sd=0.001)
        
      } else if (df$subspp[j] ==5){
        weight=plogis(2*(i-t/2))
        wa1w4n= weight*a1w4n1+(1-weight)*a1w4n2
        wa2= weight*a21+(1-weight)*a22
        walphaw4n=weight*alphaw4n1+(1-weight)*alphaw4n2
        wbetaw4n=weight*betaw4n1+(1-weight)*betaw4n2
        wa3w4n=weight*a3w4n1+(1-weight)*a3w4n2
        cf<- sum(wa1w4n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
        sizemat[j,i] <- sizemat[j,i-1]+
          rnorm(1, mean=walphaw4n + sizemat[j,i-1]*wbetaw4n + wa3w4n*cf, sd=0.001)
      } else if (sizemat[j,i] <= 0 | is.na(sizemat[j,i])){sizemat[j,i] = NA}
    }
    sizemat[edge,i] <- mean(sizemat[-edge,i-1], na.rm = T)
  }
  sizemat1=sizemat[-edge,]
  cover=((3*sizemat1[,t])/(4*pi))^(1/3)
  return(list(sizemat=sizemat1,
              area=sum(cover),
              area=max(x)*max(y)))
}
###################################################################### graph the results
a=sage_sim()[[2]] 
matplot(t(a),type="l")
#for(i in 1:t){plot(df$x,df$y,cex=(a[,i]),main=i)}
plot(df1$x,df1$y, cex=a[,t]*5, col=df1$subspp, pch=19)
#names(df1)[ncol(df1)-t:ncol(df1)]=c(1:40)
#plot
df1=cbind(df[-edge,], as.data.frame(a))

plot(a[1,], type="l",ylim=c(0,.8), main="Spacing: 1x1", 
     ylab=expression(paste("Crown volume",(m^{3}))),xaxt="n",xlab="Time (years)")
grid()
for(i in 1:nrow(a)){
  lines(a[i,], type="l", col=df1$subspp[i])
}
#xtick=seq(1,8, by=1)
axis(side=1,at=seq(0,40,by=5),labels=c(0:8))
#legend("topleft",legend=c("T2n","T4n","V2n","V4n","W4n"),col=c("green","red","blue","orange","purple"),lty=2,lwd=3,bty="n")
legend("topleft",legend=unique(df1$spp),col=unique(df1$subspp),pch=1,lty=2,lwd=3,bty="n")
#legend(0,1,legend=list("V2n","V4n","T2n","W4n","T4n"),col=unique(df1$subspp),pch=1,lty=2,lwd=3,bty="n")

