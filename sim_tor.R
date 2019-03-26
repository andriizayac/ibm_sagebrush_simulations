library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(tidybayes)
dist.fx=function(x1, x2, y1, y2) {sqrt((x2-x1)^2 + (y2-y1)^2)}
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
w_param.fx=function(timestep=NULL,param1=NULL,param2=NULL){
  #weight=plogis(2*(timestep-t/2))
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

mf=read.csv("MajorFlats_comp_df.csv")
mf3=filter(df,census=)

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


############################################ sage_sim
n=20 # n x n grid
x=rep(seq(1,l=n,by=1),times=n) # x and y reproduce orchard designx=orchard1$x[1:476]  
y=sort(rep(seq(1,l=n,by=1.5),times=n)) #y=orchard1$y[1:476]+1  
N=length(x)
ID=1:length(x)
subspp=sample(c(2:6),N, replace=T)
df<- data.frame(x,y, ID = ID,subspp=subspp)
df$spp=ifelse(df$subspp==2,"T2n",ifelse(df$subspp==3,"T4n",ifelse(df$subspp==4,"V2n",ifelse(df$subspp==5,"V4n",ifelse(df$subspp==6,"W4n",NA)))))


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
    weight=plogis(2*(i-timesteps/breakpoint))
    for(j in 1:N) {
      ############### growth
      wa1= w_param.fx(i,param_mat1[,paste0("a1_",df$subspp[j])],param_mat2[,paste0("a1_",df$subspp[j])])
      wa2= w_param.fx(i,param_mat1$a2,param_mat2$a2)
      walpha=w_param.fx(i,param_mat1[,paste0("Itype_",df$subspp[j])],param_mat2[,paste0("Itype_",df$subspp[j])])
      wbeta=w_param.fx(i,param_mat1[,paste0("c_",df$subspp[j])],param_mat2[,paste0("c_",df$subspp[j])])
      wa3=w_param.fx(i,param_mat1[,paste0("a3_",df$subspp[j])],param_mat2[,paste0("a3_",df$subspp[j])])
      sd=w_param.fx(i,param_mat1$sigma,param_mat2$sigma)
      ############### survival
      wa1_surv=param_mat_surv[,paste0("a1_",df$subspp[j])]
      wa2_surv= param_mat_surv$a2
      walpha_surv=param_mat_surv[,paste0("alpha_",df$subspp[j])]
      wbeta_surv=param_mat_surv[,paste0("beta_",df$subspp[j])]
      wa3_surv=param_mat_surv[,paste0("a3_",df$subspp[j])]
      
      cf<- sum((sizemat[,i]^wa1)/exp(tdistmat[j,]^2*wa2),na.rm=T)
      cf_surv=sum(exp(wa1_surv*log(sizemat[,i-1])-tdistmat[j,]^2*wa2_surv),na.rm=T)
      if(sizemat[j,i-1]<= 0 | is.na(sizemat[j,i-1])){ 
        sizemat[j,i-1]=NA
      }else{
        sizemat[j,i] <- (sizemat[j,i-1]+
                           rnorm(1, mean=walpha + sizemat[j,i-1]*wbeta + wa3*cf, sd=sd))*
        rbinom(1,1,plogis(walpha_surv+wbeta_surv*sizemat[j,i-1]+wa3_surv*cf_surv)^(1/timesteps))
      }
      }}
  #sizemat=ifelse(sizemat<=0,NA,sizemat)
  return(sizemat)
}
t=20
breakpoint=2
tic=Sys.time()
a=sage_sim.tor(timesteps = 21,breakpoint=2)
Sys.time()-tic


matplot(t(a),type="l")
abun.fx(smat=a)
#for(i in 1:t){plot(df1$x,df1$y,cex=(a[,i]),main=i)}
plot(df1$x,df1$y, cex=a[,19]*5, col=df1$subspp, pch=19,xlab="x",ylab="y") 
#names(df1)[ncol(df1)-t:ncol(df1)]=c(1:40)
#plot
df1=cbind(df, as.data.frame(a))
plot(a[1,], type="l",ylim=c(0,3),ylab="",xaxt="n",xlab="")
grid()
for(i in 1:nrow(a)){lines(a[i,], type="l", col=df1$subspp[i])}
mtext(expression(paste("Crown volume"~(m^{3}))),side = 2, line=2.1)
mtext("Time (years)",side=1,line=2)
title("Spacing: 1 x 1.5", line=1.2)
axis(side=1,at=seq(0,40,by=5),labels=c(0:8))
legend("topleft",legend=unique(df1$spp),col=unique(df1$subspp),pch=1,lty=2,lwd=3,bty="n")
boxplot(a[,t]~df1$spp)
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


for(i in 2:10){
  for(j in 1:N){
    print(sum(sizemat[,i-1]^wa1/exp(tdistmat[j,]^2*wa2),rm.na=T))#
}}

