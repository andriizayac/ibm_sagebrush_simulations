n=30
N=n^2
ID=1:N
t=40

spacing=seq(.5,2,l=4)
relative_cover=rep(NA,length(spacing))
total_cover=rep(NA,length(spacing))
palive=rep(NA,length(spacing))
height=rep(NA,length(spacing))



function(x,y,z){
  param_mat1=x
  param_mat2=y
  param_mat_surv=z
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
}