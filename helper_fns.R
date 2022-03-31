# === helper functions
# euclidean distance function
dist.fx <- function(x1, x2, y1, y2) {sqrt((x2-x1)^2 + (y2-y1)^2)}

# eucledean distance on the toroid wrap
toroid.dist <- function(x1,y1,x2,y2,xmax,ymax){
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

# competition kernel function based reproducing the relationship from the statistical model
cf.fn <- function(model=NULL,data=NULL){
  n = data$N
  nssp = data$k
  post = extract(model)
  iter = dim(post[[1]])[1]
  # extract parameters
  a01 = with(post,a01)
  a2 = with(post,a2)
  # spatial data
  sobs = data$size_observations
  dobs = data$dist_observations^2
  n_nb = data$n_nb
  pos = data$pos
  type = data$type
  cf_mod = matrix(NA,nrow=iter,ncol=n)
  for(i in 1:iter){
    for(j in 1:n){
      smat = sobs[pos[j] : (pos[j] + n_nb[j] - 1)]
      dmat = dobs[pos[j] : (pos[j] + n_nb[j] - 1)]
      cf_mod[i,j]=sum(smat^a01[i] / exp(dmat * a2[i]))
    }
  }
  return(cf_mod)
}

