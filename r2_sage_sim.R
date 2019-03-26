#lib and functions
library(rstan)
library(ggplot2)
library(dplyr)
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
sage_sim<-function(timesteps=t) {
  dist_vec=rep(NA,times=N)
  sizemat=matrix(NA,N,timesteps)
  sizemat[,1]=0.001
  edge=df$ID[which(df$x>mean(df$x)+10 | df$x<mean(df$x)-10 | df$y>mean(df$y)+10 | df$y<mean(df$y)-10)]
  for(i in 2:timesteps) {
    for(j in 1:N) {
      dist_vec=dist.fx(df$x[j],df$x,df$y[j],df$y)
      weight=plogis(2*(i-t/2))
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
#################################################### import models
model1=readRDS("fit_aa1_segment_c_sage_sim.rds")
post1=extract(model1)
model2=readRDS("fit_aa2_segment_c_sage_sim.rds")
post2=extract(model2)
surv=readRDS("surv_fit_orchard_a_sage_sim.rds")
spost=extract(surv)
#################################################### exptract parameters
param_mat1=as.data.frame(extract(model1))%>%select(-matches("mu|log_lik|lp__"))
param_mat2=as.data.frame(extract(model2))%>%select(-matches("mu|log_lik|lp__"))
param_mat_surv=as.data.frame(extract(surv))%>%select(-matches("theta|log_lik|lp__"))
names(param_mat1)=gsub("\\.","_",names(param_mat1))
names(param_mat2)=gsub("\\.","_",names(param_mat2))
names(param_mat_surv)=gsub("\\.","_",names(param_mat_surv))
param_mat1=as.data.frame(param_mat1%>%summarize_all(mean))
param_mat2=as.data.frame(param_mat2%>%summarize_all(mean))
param_mat_surv=as.data.frame(param_mat_surv%>%summarize_all(mean))
############################################################# Density sequence simulations
#set up the environment
n=100
N=n^2
ID=1:N
t=40

spacing=seq(.25,3,l=50)
relative_cover=rep(NA,length(spacing))
total_cover=rep(NA,length(spacing))
palive=rep(NA,length(spacing))
height=rep(NA,length(spacing))
################################# run simulations
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


sim_df=data.frame(rel_cov=relative_cover, tot_cov=total_cover, palive=palive, avg_height=height,spacing=spacing)
saveRDS(df, file="sim_df.rds")



