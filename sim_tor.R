library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(tidybayes)
dist.fx=function(x1, x2, y1, y2) {sqrt((x2-x1)^2 + (y2-y1)^2)}
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
w_param.fx=function(timestep=NULL,param1=NULL,param2=NULL, bp=breakpoint){
  weight=plogis(2*(timestep-t/bp))
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
  # this function generates a distance matrix that is wrapped on edges
  # xmax,ymax define the size of the plot in x and y 
  xcutoff=xmax/2
  ycutoff=ymax/2
  dx=abs(x2-x1)
  dy=abs(y2-y1)
  if (dx > xcutoff){
    dx=xmax-dx
  } 
  if(dy > ycutoff){
    dy=ymax-dy
  }
  return(sqrt(dx^2+dy^2))
}
#rm.zero<-function(x){x<- x[x != 0]}


######################################### IMPORT DATA
model1=readRDS("~/Desktop/Orchard/OrchardProjectR/CaughlinLab_sagebrushABM/fit_mf3_segment_c.rds")
post1=extract(model1)
model2=readRDS("~/Desktop/Orchard/OrchardProjectR/CaughlinLab_sagebrushABM/fit_mf5_segment_c.rds")
post2=extract(model2)
surv=readRDS("~/Desktop/Orchard/OrchardProjectR/CaughlinLab_sagebrushABM/surv_fit_mf_b.rds")
spost=extract(surv)


mf=read.csv("~/Desktop/Orchard/OrchardProjectR/CaughlinLab_sagebrushABM/MajorFlats_comp_df1.csv")
mf$subspp=as.numeric(as.factor(mf$subsppcyt))
mf2=filter(mf,census==2 & subsppcyt!="AR")%>%droplevels()
mf3=filter(mf, census==3 & subsppcyt!="AR")%>%droplevels()
mf5=filter(mf,census==5 & subsppcyt!="AR")%>%droplevels()


# parameters
param_mat1=as.data.frame(extract(model1))%>%select(-matches("mu|log_lik|lp__"))
param_mat2=as.data.frame(extract(model2))%>%select(-matches("mu|log_lik|lp__"))
param_mat_surv=as.data.frame(extract(surv))%>%select(-matches("theta|log_lik|lp__"))
names(param_mat1)=gsub("\\.","_",names(param_mat1))
names(param_mat2)=gsub("\\.","_",names(param_mat2))
names(param_mat_surv)=gsub("\\.","_",names(param_mat_surv))
param_mat1=as.data.frame(param_mat1%>%summarize_all(median))
param_mat2=as.data.frame(param_mat2%>%summarize_all(median))
param_mat_surv=as.data.frame(param_mat_surv%>%summarize_all(median))


############################################ sage_sim
n=20 # n x n grid
x=rep(seq(1,l=n,by=1),times=n) # x and y reproduce orchard designx=orchard1$x[1:476]  
y=sort(rep(seq(1,l=n,by=1.5),times=n)) #y=orchard1$y[1:476]+1  
N=length(x)
ID=1:length(x)
subspp=sample(c(2:6),N, replace=T)
df<- data.frame(x,y, ID = ID,subspp=subspp)
df$spp=ifelse(df$subspp==2,"t2x",ifelse(df$subspp==3,"t4x",ifelse(df$subspp==4,"v2x",ifelse(df$subspp==5,"v4x",ifelse(df$subspp==6,"w4x",NA)))))


sage_sim.tor<-function(timesteps=t,breakpoint=2) {
  #dist_vec=rep(NA,times=N)
  sizemat=matrix(NA,N,timesteps)
  sizemat[,1]=rnorm(N,mean=0.001,sd=.0001)
  tdistmat=matrix(NA,N,N)
  for(i in 1:N){
    for(j in 1:N){
      tdistmat[i,j]=toroid.dist(df$x[i],df$y[i],df$x[j],df$y[j],max(df$x),max(df$y))
    }
  }
  ifelse(tdistmat==0,NA,tdistmat)
  for(i in 2:timesteps) {
    cat(" ",i,"/",timesteps,sep="")
    #cat(" ",format(Sys.time()))
    for(j in 1:N) {
      ############### growth
      wa1= w_param.fx(i,param_mat1[,paste0("a1_",df$subspp[j])],param_mat2[,paste0("a1_",df$subspp[j])], bp=breakpoint)
      wa2= w_param.fx(i,param_mat1$a2,param_mat2$a2, bp=breakpoint)
      walpha=w_param.fx(i,param_mat1[,paste0("Itype_",df$subspp[j])],param_mat2[,paste0("Itype_",df$subspp[j])], bp=breakpoint)
      wbeta=w_param.fx(i,param_mat1[,paste0("c_",df$subspp[j])],param_mat2[,paste0("c_",df$subspp[j])], bp=breakpoint)
      wa3=w_param.fx(i,param_mat1[,paste0("a3_",df$subspp[j])],param_mat2[,paste0("a3_",df$subspp[j])], bp=breakpoint)
      sd=w_param.fx(i,param_mat1$sigma,param_mat2$sigma, bp=breakpoint)
      ############### survival
      wa1_surv=param_mat_surv[,paste0("a1_",df$subspp[j])]
      wa2_surv= param_mat_surv$a2
      walpha_surv=param_mat_surv[,paste0("alpha_",df$subspp[j])]
      wbeta_surv=param_mat_surv[,paste0("beta_",df$subspp[j])]
      wa3_surv=param_mat_surv[,paste0("a3_",df$subspp[j])]
      
      cf<- sum((sizemat[,i-1]^wa1)/exp(tdistmat[j,]^2*wa2),na.rm=T)
      cf_surv=sum(exp(wa1_surv*log(sizemat[,i-1])-tdistmat[j,]^2*wa2_surv),na.rm=T)
      if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ 
        sizemat[j,i-1]=NA
      }else{
        rd=rnorm(1, mean=walpha + sizemat[j,i-1]*wbeta + wa3*cf, sd=sd)
        Pr_surv=plogis(walpha_surv+wbeta_surv*sizemat[j,i-1]+wa3_surv*cf_surv)^(1/timesteps)
        sizemat[j,i] <- (sizemat[j,i-1]+rd)*rbinom(1,1,Pr_surv)
      }
      }}
  #sizemat=ifelse(sizemat<=0,NA,sizemat)
  return(sizemat)
}
#t=20
#breakpoint=3
tic=Sys.time()
a=sage_sim.tor(timesteps = 30,breakpoint=1.4)
Sys.time()-tic
#matplot(t(a),type="l")
#abun.fx(smat=a)
#for(i in 1:t){plot(df1$x,df1$y,cex=(a[,i]),main=i)}
#plot(df1$x,df1$y, cex=a[,19]*5, col=df1$subspp, pch=19,xlab="x",ylab="y") 
#names(df1)[ncol(df1)-t:ncol(df1)]=c(1:40)
#plot
a=sage_sim.tor(timesteps = 30,breakpoint=1.3)
layout(matrix(c(1,2,3,4),nrow=1))
par(mar=c(4,4,3,2))
df1=cbind(df, as.data.frame(a))
plot(a[1,], type="l",ylim=c(0,4),ylab="",xaxt="n",xlab="")
grid()
for(i in 1:nrow(a)){lines(a[i,], type="l", col=df1$subspp[i])}
mtext(expression(paste("Crown volume"~(m^{3}))),side = 2, line=2.1)
mtext("Time (years)",side=1,line=2)
title("Spacing: 1 x 1.5 | Sim", line=1.2)
#axis(side=1,at=seq(0,21,by=3),labels=c(0:7))
legend("topleft",legend=unique(df1$spp),col=unique(df1$subspp),pch=1,lty=2,lwd=3,bty="n")

boxplot(df1[,ncol(df1)]~df1$subspp, col=levels(as.factor(df1$subspp)),xaxt="n",ylim=c(0,4),
        main="Spacing: 1 x 1.5 | Sim")
axis(side=1,at=seq(1,5,by=1),labels=c("t2x","t4x","v2x","v4x","w4x"))
#####
#par(mfrow=c(1,1))
boxplot(mf3$size_t~mf3$subsppcyt, col=levels(as.factor(mf3$subspp)),xaxt="n",ylim=c(0,4),
        main="MajorsFlat2012")
axis(side=1,at=seq(1,5,by=1),labels=c("t2x","t4x","v2x","v4x","w4x"))

boxplot(mf5$size_t~mf5$subsppcyt, col=levels(as.factor(mf5$subspp)),xaxt="n",ylim=c(0,4),
        main="MajorsFlat2018")
axis(side=1,at=seq(1,5,by=1),labels=c("t2x","t4x","v2x","v4x","w4x"))

#####################################################################
df1=df              #CHECK NNUMBER OF NEIGHBORS
df1$y=df1$y+1
tdistmat=matrix(NA, N,N)
for(i in 1:N){
  for(j in 1:N){
    tdistmat[i,j]=toroid.dist(df1$x[i],df1$y[i],df1$x[j],df1$y[j],max(df1$x),max(df1$y))
  }
}
tna=ifelse(tdistmat>2 | tdistmat==0,NA,tdistmat)
tna_na=apply(tna,1,function(x){sum(!is.na(x))})
a=df1[tna_na<8,]
a=df1[!is.na(tna[60,]),]
plot(1,type="n",xlim=c(0,28),ylim=c(0,25))
points(a$x,a$y)
##################################################################
######## TWEAKING BREAKPOINT AND TIME PARAMETERS
########################################################

a=sage_sim.tor(timesteps = 28,breakpoint=1)
layout(matrix(c(1,2,3,4),nrow=1))
par(mar=c(4,4,3,2))
df1=cbind(df, as.data.frame(a))
plot(a[1,], type="l",ylim=c(0,4),ylab="",xaxt="n",xlab="")
grid()
for(i in 1:nrow(a)){lines(a[i,], type="l", col=df1$subspp[i])}
mtext(expression(paste("Crown volume"~(m^{3}))),side = 2, line=2.1)
mtext("Time",side=1,line=2)
title("Spacing: 1 x 1.5 | Sim", line=1.5)
#axis(side=1,at=seq(0,21,by=3),labels=c(0:7))
legend("topleft",legend=unique(df1$spp),col=unique(df1$subspp),pch=1,lty=2,lwd=3,bty="n")

boxplot(df1[,ncol(df1)]~df1$subspp, col=levels(as.factor(df1$subspp)),xaxt="n",ylim=c(0,4),
        main="Spacing: 1 x 1.5 | Sim")
axis(side=1,at=seq(1,5,by=1),labels=c("t2x","t4x","v2x","v4x","w4x"))
#####
#par(mfrow=c(1,1))
boxplot(mf3$size_t~mf3$subsppcyt, col=levels(as.factor(mf3$subspp)),xaxt="n",ylim=c(0,4),
        main="MajorsFlat2012")
axis(side=1,at=seq(1,5,by=1),labels=c("t2x","t4x","v2x","v4x","w4x"))

boxplot(mf5$size_t~mf5$subsppcyt, col=levels(as.factor(mf5$subspp)),xaxt="n",ylim=c(0,4),
        main="MajorsFlat2018")
axis(side=1,at=seq(1,5,by=1),labels=c("t2x","t4x","v2x","v4x","w4x"))


#############################################################
###### POSTERIOR LOOP
#############################################################


#################################################### extract parameters
param_mat1=as.data.frame(extract(model1))%>%select(-matches("mu|log_lik|lp__"))
param_mat2=as.data.frame(extract(model2))%>%select(-matches("mu|log_lik|lp__"))
param_mat_surv=as.data.frame(extract(surv))%>%select(-matches("theta|log_lik|lp__"))
names(param_mat1)=gsub("\\.","_",names(param_mat1))
names(param_mat2)=gsub("\\.","_",names(param_mat2))
names(param_mat_surv)=gsub("\\.","_",names(param_mat_surv))

##### posterior uncertainty loop simulation function: full model
sage_sim.tor2012.v1F<-function(timesteps=timesteps) {
  
  #dist_vec=rep(NA,times=N)
  sizematF=sizemat #define output matrix for a full model
  tdistmat=matrix(NA,N,N)
  for(i in 1:N){
    for(j in 1:N){
      tdistmat[i,j]=toroid.dist(df$x[i],df$y[i],df$x[j],df$y[j],max(df$x),max(df$y))
    }
  }
  ifelse(tdistmat==0,NA,tdistmat)
  for(i in 2:timesteps) {
    #cat(" ",i,"/",timesteps,sep="")
    #cat(" ",format(Sys.time()))
    for(j in 1:N) {
      ############### growth
      wa1=param_mat1[k,paste0("a1_",df$subspp[j])]
      wa2=param_mat1$a2[k]
      walpha=param_mat1[k,paste0("Itype_",df$subspp[j])]
      wbeta=param_mat1[k,paste0("c_",df$subspp[j])]
      wa3=param_mat1[k,paste0("a3_",df$subspp[j])]
      ############### survival: took out as only 3 plants died during that interval
      
      cf<-sum((sizemat[,i-1]^wa1)/exp(tdistmat[j,]^2*wa2),na.rm=T)
      if(sizematF[j,i-1]<= 0 | is.na(sizematF[j,i-1])){ 
        sizematF[j,i-1]=NA
      }else{
        rd=walpha + sizematF[j,i-1]*wbeta + wa3*cf
        sizematF[j,i] <- (sizematF[j,i-1]+rd)
      }
    }}
  #sizemat=ifelse(sizemat<=0,NA,sizemat)
  return(sizematF)
}
#### posterior uncertainty loop simulation function: base model (no crowding)
sage_sim.tor2012.v1B<-function(timesteps=timesteps) {
  sizematB=sizemat #define output matrix for a full model
  
  for(i in 2:timesteps) {
    #cat(" ",i,"/",timesteps,sep="")
    #cat(" ",format(Sys.time()))
    for(j in 1:N) {
      ############### growth
      wa1=param_mat1[k,paste0("a1_",df$subspp[j])]
      wa2=param_mat1$a2[k]
      walpha=param_mat1[k,paste0("Itype_",df$subspp[j])]
      wbeta=param_mat1[k,paste0("c_",df$subspp[j])]
      wa3=param_mat1[k,paste0("a3_",df$subspp[j])]
      ############### survival: took out as only 3 plants died during that interval
      
      #cf<-sum((sizemat[,i-1]^wa1)/exp(tdistmat[j,]^2*wa2),na.rm=T)
      if(sizematB[j,i-1]<= 0 | is.na(sizematB[j,i-1])){ 
        sizematB[j,i-1]=NA
      }else{
        rd=walpha + sizematB[j,i-1]*wbeta
        sizematB[j,i] <- (sizematB[j,i-1]+rd)
      }
    }}
  return(sizematB)
}

## this loop is used to go through the posterior
#####
mlist=list()

 
for(k in 1:10){
  cat(" ",k,"/",timesteps,sep="")
  n=20 # n x n grid
  timesteps=12 # number of timesteps
  x=rep(seq(1,l=n,by=1),times=n) 
  y=sort(rep(seq(1,l=n,by=1.5),times=n))
  N=length(x)
  ID=1:length(x)
  subspp=sample(c(2:6),N, replace=T) # randomly assign subspecies identity
  df<- data.frame(x,y, ID = ID,subspp=subspp) # a df with initial conditions and spatial arrangement, which is different for every posterior combination
  df$spp=ifelse(df$subspp==2,"t2x",ifelse(df$subspp==3,"t4x",ifelse(df$subspp==4,"v2x",ifelse(df$subspp==5,"v4x",ifelse(df$subspp==6,"w4x",NA)))))
  
  sizemat=matrix(NA,N,timesteps)
  sizemat[,1]=rnorm(N,mean=0.001,sd=.0001) #a data matrix with reasonable intial size conditions
  #this introduces random size variation
  
  sizematF=sage_sim.tor2012.v1F(timesteps)
  sizematB=sage_sim.tor2012.v1B(timesteps)
  
  df$initial_size=sizemat[,1]
  df$sizeF=sizematF[,timesteps]
  df$sizeB=sizematB[,timesteps]
  df$diff=df$sizeB-df$sizeF
  m=lm(df$diff~df$spp)
  mlist[[k]]=summary(m)
}
# this function is good but introduces noise via rnorm(), instead the uncertainty should come from the posterior directly
sage_sim.tor2012<-function(timesteps=t) {
    #dist_vec=rep(NA,times=N)
    sizemat=matrix(NA,N,timesteps)
    sizemat[,1]=rnorm(N,mean=0.001,sd=.0001)
    tdistmat=matrix(NA,N,N)
    for(i in 1:N){
      for(j in 1:N){
        tdistmat[i,j]=toroid.dist(df$x[i],df$y[i],df$x[j],df$y[j],max(df$x),max(df$y))
      }
    }
    ifelse(tdistmat==0,NA,tdistmat)
    for(i in 2:timesteps) {
      cat(" ",i,"/",timesteps,sep="")
      #cat(" ",format(Sys.time()))
      for(j in 1:N) {
        ############### growth
        wa1= rnorm(1,param_mat1[,paste0("a1_",df$subspp[j])][k], param_mat1$sigma_a1[k])
        wa2= param_mat1$a2[k]
        walpha=rnorm(1,param_mat1[,paste0("Itype_",df$subspp[j])][k],param_mat1$sigma_Itype[k])
        wbeta=rnorm(1,param_mat1[,paste0("c_",df$subspp[j])][k],param_mat1$sigma_c[k])
        wa3=rnorm(1,param_mat1[,paste0("a3_",df$subspp[j])][k],param_mat1$sigma_a3[k])
        sd=param_mat1$sigma[k]
        ############### survival: took out as only 3 plants died during that interval
        
        cf<-sum((sizemat[,i-1]^wa1)/exp(tdistmat[j,]^2*wa2),na.rm=T)
        if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ 
          sizemat[j,i-1]=NA
        }else{
          rd=rnorm(1, mean=walpha + sizemat[j,i-1]*wbeta + wa3*cf, sd=sd)
          sizemat[j,i] <- (sizemat[j,i-1]+rd)
        }
      }}
    #sizemat=ifelse(sizemat<=0,NA,sizemat)
    return(sizemat)
  }
  


