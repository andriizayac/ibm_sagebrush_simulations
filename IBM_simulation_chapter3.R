library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(tidybayes)
dist.fx=function(x1, x2, y1, y2) {sqrt((x2-x1)^2 + (y2-y1)^2)}
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
w_param.fx=function(timestep=NULL,param1=NULL,param2=NULL){
  weight=plogis(2*(timestep-t/2))
  wparam=(1-weight)*param1+weight*param2
  return(wparam)
}
scale_surv.fx=function(x){plogis(x)^(1/5*12)}
abun.fx=function(smat=smat){
  r=((3*smat[,t])/(4*pi))^(1/3)
  crown_area=pi*r^2
  total_canopy=sum(crown_area,na.rm = T)
  plot_area=400
  return(list(total_canopy=total_canopy, #square meters?
              plot_area=plot_area,
              relative_cover=total_canopy/plot_area,
              palive=sum(is.na(smat[,t])==F)/nrow(smat),
              height=mean(r*2,na.rm=T)))
}
toroid.dist=function(x1,y1,x2,y2,xmax,ymax){
  xcutoff=xmax/2
  ycutoff=ymax/2
  dx=abs(x2-x1)
  dy=abs(y2-y1)
  if (dx > xcutoff){
    dx=1-dx
  } 
  if(dy > ycutoff){
    dy=1-dy
  }
  return(sqrt(dx^2+dy^2))
}
#rm.zero<-function(x){x<- x[x != 0]}

# import data
model1=readRDS("fit_mf3_segment_c_sage_sim.rds")
post1=extract(model1)
model2=readRDS("fit_mf5_segment_c_sage_sim.rds")
post2=extract(model2)
surv=readRDS("surv_fit_mf_a_sage_sim.rds")
spost=extract(surv)
plot(surv,pars=c("alpha"))
orchard1=read.csv("orchard1.csv")
mf=read.csv("MajorFlats_comp_df.csv")

# parameters
param_mat1=as.data.frame(extract(model1))%>%select(-matches("mu|log_lik|lp__"))
param_mat2=as.data.frame(extract(model2))%>%select(-matches("mu|log_lik|lp__"))
param_mat_surv=as.data.frame(extract(surv))%>%select(-matches("theta|log_lik|lp__"))
names(param_mat1)=gsub("\\.","_",names(param_mat1))
names(param_mat2)=gsub("\\.","_",names(param_mat2))
names(param_mat_surv)=gsub("\\.","_",names(param_mat_surv))
param_mat1=as.data.frame(param_mat1%>%summarize_all(mean))
param_mat2=as.data.frame(param_mat2%>%summarize_all(mean))
param_mat_surv=as.data.frame(param_mat_surv%>%summarize_all(mean))



model1%>%
  spread_draws(a3[subspp])%>%#str()%>%
  mutate(subspecies=as.factor(as.character(subspp)))%>%
  #median_qi()%>%
  ggplot(aes(y=subspecies,x=a3))+
  geom_halfeyeh(fill="gray65",.width=c(0.34,0.68))+#+stat_pointintervalh(.width = c(0.8,0.95))
  #geom_violinh(color=NA,fill="gray65")+stat_pointintervalh(.width = c(0.8,0.95))+
  theme_bw()




############################################################# Density sequence simulations
  n=30
  N=n^2
  ID=1:N
  t=40
  
  spacing=seq(.5,2,l=4)
  relative_cover=rep(NA,length(spacing))
  total_cover=rep(NA,length(spacing))
  palive=rep(NA,length(spacing))
  height=rep(NA,length(spacing))

  for(i in 1:length(spacing)){
    xspace=spacing[i]
    yspace=spacing[i]
    x=rep(seq(1,l=n,by=xspace),times=n)
    y=sort(rep(seq(1,l=n,by=yspace),times=n))
    distmat=matrix(NA,sqrt(N),sqrt(N))
    subspp=sample(rep(1:5,times=N/5),N)
    df<- data.frame(x,y, "vol" = 0.001, ID = ID,subspp=subspp)
    df$spp=ifelse(df$subspp==1,"T2n",ifelse(df$subspp==2,"T4n",
                                            ifelse(df$subspp==3,"V2n",
                                                   ifelse(df$subspp==4,"V4n",
                                                          ifelse(df$subspp==5,"W4n",NA)))))
    a=sage_sim()
    total_cover[i]=abun.fx(smat=a)[[1]]
    relative_cover[i]=abun.fx(smat = a)[[3]]
    palive[i]=abun.fx(smat=a)[[4]]
    height[i]=abun.fx(smat=a)[[5]]
  }

#df=readRDS("sim_df.rds")
plot(relative_cover~spacing, ylab="", xlab="Spacing (m)",col="purple",pch=19)
plot(palive~spacing, col="red",pch=19)
plot(height~spacing, col="blue",pch=19)
legend("right",legend=c("# alive individuals","% cover"),col=c("red","blue"),pch=1,lty=2,lwd=3,bty="n")

apply(outcome,2,function(x){plot(df$x,df$y,cex=x)})
####################################################################### sage_sim
n=50 # n x n gridx
x=rep(seq(1,l=n,by=1),times=n)
y=sort(rep(seq(1,l=n,by=1),times=n))
N=length(x)
ID=1:length(x)
subspp=sample(rep(1:5,times=N/5),N)
df<- data.frame(x,y, "vol" = rnorm(length(x),.001,sd=0.0001), ID = ID,subspp=subspp)
df$spp=ifelse(df$subspp==1,"T2n",ifelse(df$subspp==2,"T4n",ifelse(df$subspp==3,"V2n",ifelse(df$subspp==4,"V4n",ifelse(df$subspp==5,"W4n",NA)))))
edge=df$ID[which(df$x>mean(df$x)+10 | df$x<mean(df$x)-10 | df$y>mean(df$y)+10 | df$y<mean(df$y)-10)]


sage_sim<-function(timesteps=t) {
  dist_vec=rep(NA,times=N)
  sizemat=matrix(NA,N,timesteps)
  sizemat[,1]=rnorm(N,0.001,.0001)
  edge=df$ID[which(df$x>mean(df$x)+10 | df$x<mean(df$x)-10 | df$y>mean(df$y)+10 | df$y<mean(df$y)-10)]
  for(i in 2:timesteps) {
    for(j in 1:N) {
      dist_vec=dist.fx(df$x[j],df$x,df$y[j],df$y)
      weight=plogis(2*(i-timesteps/2))
      ############### growth
      wa1= w_param.fx(i,param_mat1[,paste("a1_",df$subspp[j],sep="")],param_mat2[,paste("a1_",df$subspp[j],sep="")])
      wa2= w_param.fx(i,param_mat1$a2,param_mat2$a2)
      walpha=w_param.fx(i,param_mat1[,paste("Itype_",df$subspp[j],sep="")],param_mat2[,paste("Itype_",df$subspp[j],sep="")])
      wbeta=w_param.fx(i,param_mat1[,paste("c_",df$subspp[j],sep="")],param_mat2[,paste("c_",df$subspp[j],sep="")])
      wa3=w_param.fx(i,param_mat1[,paste("a3_",df$subspp[j],sep="")],param_mat2[,paste("a3_",df$subspp[j],sep="")])
      sd=w_param.fx(i,param_mat1$sigma,param_mat2$sigma)
      ############### survival
      wa1_surv=param_mat_surv[,paste("a1_",df$subspp[j],sep="")]
      wa2_surv= param_mat_surv$a2
      walpha_surv=param_mat_surv[,paste("alpha_",df$subspp[j],sep="")]
      wbeta_surv=param_mat_surv[,paste("beta_",df$subspp[j],sep="")]
      wa3_surv=param_mat_surv[,paste("a3_",df$subspp[j],sep="")]
      cf<- sum((sizemat[,i-1]^wa1)/exp(dist_vec^2*wa2),na.rm=TRUE)
      cf_surv=sum(exp(wa1_surv*log(sizemat[,i-1])-dist_vec^2*wa2_surv))
      if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ sizemat[j,i]=NA
      }else{
        sizemat[j,i] <- (sizemat[j,i-1]+
                           rnorm(1, mean=walpha + sizemat[j,i-1]*wbeta + wa3*cf, sd=sd))
        rbinom(1,1,plogis(walpha_surv+wbeta_surv*sizemat[j,i-1]+wa3_surv*cf_surv)^(1/60))
      }
      sizemat[edge,i] <- mean(sizemat[-edge,i-1], na.rm = T)
    }}
  sizemat=ifelse(sizemat<0,NA,sizemat)
  return(sizemat[-edge,])
}




abun.fx=function(smat=smat){
    r=((3*smat[,t])/(4*pi))^(1/3)
    crown_area=pi*r^2
    total_canopy=sum(crown_area,na.rm = T)
    plot_area=400
    return(list(total_canopy=total_canopy, #square meters?
                  plot_area=plot_area,
                  relative_cover=total_canopy/plot_area,
                  palive=sum(is.na(smat[,t])==F)/nrow(smat),
                  height=mean(r*2,na.rm=T)))
}
###################################################################### graph the results
tic=Sys.time()
a=sage_sim()
Sys.time()-tic


matplot(t(a),type="l")
abun.fx(smat=a)
#for(i in 1:t){plot(df1$x,df1$y,cex=(a[,i]),main=i)}
plot(df1$x,df1$y, cex=a[,19]*5, col=df1$subspp, pch=19,xlab="x",ylab="y") 
#names(df1)[ncol(df1)-t:ncol(df1)]=c(1:40)
#plot
df1=cbind(df[-edge,], as.data.frame(a))
plot(a[1,], type="l",ylim=c(0,1.5),ylab="",xaxt="n",xlab="")
grid()
for(i in 1:nrow(a)){lines(a[i,], type="l", col=df1$subspp[i])}
mtext(expression(paste("Crown volume"~(m^{3}))),side = 2, line=2.1)
mtext("Time (years)",side=1,line=2)
title("Spacing: 1 x 1", line=1.2)
axis(side=1,at=seq(0,40,by=5),labels=c(0:8))
legend("topleft",legend=unique(df1$spp),col=unique(df1$subspp),pch=1,lty=2,lwd=3,bty="n")



a=df[edge,]
plot(a$x,a$y)

############### older versions of the simulation functions
#sage_sim<-function(timesteps=t) {
  dist_vec=rep(NA,times=N)
  sizemat=matrix(NA,N,timesteps)
  sizemat[,1]=0.001
  edge=df$ID[which(df$x>mean(df$x)+10 | df$x<mean(df$x)-10 | df$y>mean(df$y)+10 | df$y<mean(df$y)-10)]
  for(i in 2:timesteps) {
    for(j in 1:N) {
      dist_vec=dist.fx(df$x[j],df$x,df$y[j],df$y)
      if(df$subspp[j] == 1){
        weight=plogis(2*(i-t/2))
        wa1t2n= (1-weight)*a1t2n1+weight*a1t2n2
        wa2= (1-weight)*a21+weight*a22
        walphat2n=(1-weight)*alphat2n1+weight*alphat2n2
        wbetat2n=(1-weight)*betat2n1+weight*betat2n2
        wa3t2n=(1-weight)*a3t2n1+weight*a3t2n2
        sd=(1-weight)*sd1+weight*sd2 #################################################
        cf<- sum(wa1t2n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
        if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ sizemat[j,i]=NA
        }else{ 
          sizemat[j,i] <- sizemat[j,i-1]+
            rnorm(1, mean=walphat2n + sizemat[j,i-1]*wbetat2n + wa3t2n*cf, sd=sd)
        }
      } else if (df$subspp[j] ==2){
        weight=plogis(2*(i-t/2))
        wa1t4n= (1-weight)*a1t4n1+weight*a1t4n2
        wa2= (1-weight)*a21+weight*a22
        walphat4n=(1-weight)*alphat4n1+weight*alphat4n2
        wbetat4n=(1-weight)*betat4n1+weight*betat4n2
        wa3t4n=(1-weight)*a3t4n1+weight*a3t4n2
        sd=(1-weight)*sd1+weight*sd2
        cf<- sum(wa1t4n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
        if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ sizemat[j,i]=NA
        }else{  
          sizemat[j,i] <- sizemat[j,i-1]+
            rnorm(1, mean=walphat4n + sizemat[j,i-1]*wbetat4n + wa3t4n*cf, sd=sd)
        }
      } else if (df$subspp[j] ==3){
        weight=plogis(2*(i-t/2))
        wa1v2n= (1-weight)*a1v2n1+weight*a1v2n2
        wa2= (1-weight)*a21+weight*a22
        walphav2n=(1-weight)*alphav2n1+weight*alphav2n2
        wbetav2n=(1-weight)*betav2n1+weight*betav2n2
        wa3v2n=(1-weight)*a3v2n1+weight*a3v2n2
        sd=(1-weight)*sd1+weight*sd2
        cf<- sum(wa1v2n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
        if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ sizemat[j,i]=NA
        }else{   
          sizemat[j,i] <- sizemat[j,i-1]+
            rnorm(1, mean=walphav2n + sizemat[j,i-1]*wbetav2n + wa3v2n*cf, sd=sd)
        }
      } else if (df$subspp[j] ==4){
        weight=plogis(2*(i-t/2))
        wa1v4n= (1-weight)*a1v4n1+weight*a1v4n2
        wa2= (1-weight)*a21+weight*a22
        walphav4n=(1-weight)*alphav4n1+weight*alphav4n2
        wbetav4n=(1-weight)*betav4n1+weight*betav4n2
        wa3v4n=(1-weight)*a3v4n1+weight*a3v4n2
        sd=(1-weight)*sd1+weight*sd2
        cf<- sum(wa1v4n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
        if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ sizemat[j,i]=NA
        }else{ 
          sizemat[j,i] <- sizemat[j,i-1]+
            rnorm(1, mean=walphav4n + sizemat[j,i-1]*wbetav4n + wa3v4n*cf, sd=sd)
        }
      } else if (df$subspp[j] ==5){
        weight=plogis(2*(i-t/2))
        wa1w4n= (1-weight)*a1w4n1+weight*a1w4n2
        wa2= (1-weight)*a21+weight*a22
        walphaw4n=(1-weight)*alphaw4n1+weight*alphaw4n2
        wbetaw4n=(1-weight)*betaw4n1+weight*betaw4n2
        wa3w4n=(1-weight)*a3w4n1+weight*a3w4n2
        sd=(1-weight)*sd1+weight*sd2
        cf<- sum(wa1w4n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
        if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ sizemat[j,i]=NA
        }else{ 
          sizemat[j,i] <- sizemat[j,i-1]+
            rnorm(1, mean=walphaw4n + sizemat[j,i-1]*wbetaw4n + wa3w4n*cf, sd=sd)
        }
      }
    }
    sizemat[edge,i] <- mean(sizemat[-edge,i-1], na.rm = T)
  }
  sizemat=ifelse(sizemat<0,NA,sizemat)
  return(sizemat[-edge,])
}
#plogis(0)^(1/(5*12)) #SCALE SURVIVAL!!! 
#sage_surv=function(sizemat=sizemat){
  dist_vec=rep(NA,times=nrow(sizemat))
  for(j in 1:nrow(sizemat)){
    if(df$subspp[j] == 1 & any(is.na(sizemat[j,]))==F){
      cf_surv=sum(a1t2n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
      sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphat2n_surv+sizemat[j,t/2]*betat2n_surv+a3t2n_surv*cf_surv))
    } else if (df$subspp[j] == 2 & any(is.na(sizemat[j,]))==F){
      cf_surv=sum(a1t4n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
      sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphat4n_surv+sizemat[j,t/2]*betat4n_surv+a3t4n_surv*cf_surv))
    } else if (df$subspp[j] == 3 & any(is.na(sizemat[j,]))==F){
      cf_surv=sum(a1v2n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
      sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphav2n_surv+sizemat[j,t/2]*betav2n_surv+a3v2n_surv*cf_surv))
    } else if (df$subspp[j] == 4 & any(is.na(sizemat[j,]))==F){
      cf_surv=sum(a1v4n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
      sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphav4n_surv+sizemat[j,t/2]*betav4n_surv+a3v4n_surv*cf_surv))
    } else if (df$subspp[j] == 5 & any(is.na(sizemat[j,]))==F){
      cf_surv=sum(a1w4n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
      sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphaw4n_surv+sizemat[j,t/2]*betaw4n_surv+a3w4n_surv*cf_surv))
    }
  }
  sizemat=ifelse(sizemat==0,NA,sizemat)
  return(sizemat)
}
#sage1backup<-function(timesteps=t) {
dist_vec=rep(NA,times=N)
sizemat=matrix(NA,N,timesteps)
sizemat[,1]=0.001
edge=df$ID[which(df$x>mean(df$x)+10 | df$x<mean(df$x)-10 | df$y>mean(df$y)+10 | df$y<mean(df$y)-10)]
for(i in 2:timesteps) {
  for(j in 1:N) {
    dist_vec=dist.fx(df$x[j],df$x,df$y[j],df$y)
    weight=plogis(2*(i-t/2))
    wa1= (1-weight)*param_mat1[,paste("a1_",df$subspp[j],sep="")]+weight*param_mat2[,paste("a1_",df$subspp[j],sep="")]
    wa2= (1-weight)*param_mat1$a2+weight*param_mat2$a2
    walpha=(1-weight)*param_mat1[,paste("Itype_",df$subspp[j],sep="")]+weight*param_mat2[,paste("Itype_",df$subspp[j],sep="")]
    wbeta=(1-weight)*param_mat1[,paste("c_",df$subspp[j],sep="")]+weight*param_mat2[,paste("c_",df$subspp[j],sep="")]
    wa3=(1-weight)*param_mat1[,paste("a3_",df$subspp[j],sep="")]+weight*param_mat2[,paste("a3_",df$subspp[j],sep="")]
    sd=(1-weight)*param_mat1$sigma+weight*param_mat2$sigma
    cf<- sum(wa1t2n*sizemat[,i-1]/exp(dist_vec^2*wa2),na.rm=TRUE)
    if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ sizemat[j,i]=NA
    }else{ 
      sizemat[j,i] <- sizemat[j,i-1]+
        rnorm(1, mean=walphat2n + sizemat[j,i-1]*wbetat2n + wa3t2n*cf, sd=sd)
    }
    sizemat[edge,i] <- mean(sizemat[-edge,i-1], na.rm = T)
  }}
sizemat=ifelse(sizemat<0,NA,sizemat)
return(sizemat[-edge,])
}
#sage_surv=function(sizemat=sizemat){
dist_vec=rep(NA,times=nrow(sizemat))
for(j in 1:nrow(sizemat)){
  if(df$subspp[j] == 1 & any(is.na(sizemat[j,]))==F){
    cf_surv=sum(a1t2n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
    sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphat2n_surv+sizemat[j,t/2]*betat2n_surv+a3t2n_surv*cf_surv))
  } else if (df$subspp[j] == 2 & any(is.na(sizemat[j,]))==F){
    cf_surv=sum(a1t4n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
    sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphat4n_surv+sizemat[j,t/2]*betat4n_surv+a3t4n_surv*cf_surv))
  } else if (df$subspp[j] == 3 & any(is.na(sizemat[j,]))==F){
    cf_surv=sum(a1v2n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
    sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphav2n_surv+sizemat[j,t/2]*betav2n_surv+a3v2n_surv*cf_surv))
  } else if (df$subspp[j] == 4 & any(is.na(sizemat[j,]))==F){
    cf_surv=sum(a1v4n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
    sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphav4n_surv+sizemat[j,t/2]*betav4n_surv+a3v4n_surv*cf_surv))
  } else if (df$subspp[j] == 5 & any(is.na(sizemat[j,]))==F){
    cf_surv=sum(a1w4n_surv*sizemat[,t/2]/exp(dist_vec^2*a2_surv),na.rm=TRUE)
    sizemat[j,]=sizemat[j,]*rbinom(1,1,plogis(alphaw4n_surv+sizemat[j,t/2]*betaw4n_surv+a3w4n_surv*cf_surv))
  }
}
sizemat=ifelse(sizemat==0,NA,sizemat)
return(sizemat)
}
#sim=function(xspace=5,yspace=5,t=40){
  x=rep(seq(1,l=20,by=xspace),times=20)
  y=sort(rep(seq(1,l=20,by=yspace),times=20))
  
  ID=c(1:N)
  subspp=sample(rep(1:5,times=N/5),N)
  df<- data.frame(x,y, "vol" = 0.0001, ID = ID,subspp=subspp)
  df$spp=ifelse(df$subspp==1,"T2n",ifelse(df$subspp==2,"T4n",ifelse(df$subspp==3,"V2n",ifelse(df$subspp==4,"V4n",ifelse(df$subspp==5,"W4n",NA)))))
  a=sage_sim()
  return(sage_sim()[[2]])
  #print(df)
}  # crap!!!

################################################# parameters extracted but not needed for now
# variance
sd1=mean(post1$sigma)
sd2=mean(post2$sigma)

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


# Interaction parameters SURVIVAL
alpha_surv=mean(spost$a0)
alphat2n_surv=mean(spost$alpha[,1])
alphat4n_surv=mean(spost$alpha[,2])
alphav2n_surv=mean(spost$alpha[,3])
alphav4n_surv=mean(spost$alpha[,4])
alphaw4n_surv=mean(spost$alpha[,5])
beta_surv=mean(spost$beta0)
betat2n_surv=mean(spost$beta[,1])
betat4n_surv=mean(spost$beta[,2])
betav2n_surv=mean(spost$beta[,3])
betav4n_surv=mean(spost$beta[,4])
betaw4n_surv=mean(spost$beta[,5])
a3_surv=mean(spost$a03)
a3t2n_surv=mean(spost$a3[,1])
a3t4n_surv=mean(spost$a3[,2])
a3v2n_surv=mean(spost$a3[,3])
a3v4n_surv=mean(spost$a3[,4])
a3w4n_surv=mean(spost$a3[,5])
a2_surv=mean(spost$a2)
a1_surv=mean(spost$a01)
a1t2n_surv=mean(spost$a1[,1])
a1t4n_surv=mean(spost$a1[,2])
a1v2n_surv=mean(spost$a1[,3])
a1v4n_surv=mean(spost$a1[,4])
a1w4n_surv=mean(spost$a1[,5])
