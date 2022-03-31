# Script to fit a log-log NDD model for MajorsFlat CG


pkgs <- c("dplyr", "lme4")
sapply(pkgs, require, character.only = T)
# === helper functions
scale.fn <- function(x) { (x - mean(x))/(2*sd(x)) }
nb.fn <- function(x, dmat) {
  smat <- matrix(rep(x, length(x)), length(x))
  smat <- ifelse(is.na(dmat), 0, smat)
  rowSums(smat/exp(dmat^2*.1), na.rm = T)
}
toroid.dist <- function(x1,y1,x2,y2,xmax,ymax){
  # this function generates a distance matrix that is wrapped on edges
  # xmax,ymax define the size of the plot in x and y 
  
  xcutoff <- xmax/2
  ycutoff <- ymax/2
  dx <- abs(x2-x1)
  dy <- abs(y2-y1)
  if (dx > xcutoff){
    dx <- xmax-dx
  } 
  if(dy > ycutoff){
    dy <- ymax-dy
  }
  return(sqrt(dx^2 + dy^2))
}


# load data
dat <- read.csv("MajorFlats_comp_df1.csv") %>%
  filter(subsppcyt != "AR")
outs <- dat$Sample.tag.[which(is.na(dat$size_t) | dat$size_t == 0)]
dat <- dat[!c(dat$Sample.tag. %in% outs),]

y <- dat$size_t[dat$census != 1]
x <- dat$size_t[dat$census !=max(dat$census)]
type <- dat$subsppcyt[dat$census != 1]

# --- create distance matrix
dat.nb1 <- dat[dat$census == 1, ]
n = nrow(dat.nb1)
dat.nb <- dat[dat$census != 1, ]
distmat <- as.matrix( dist(cbind(dat.nb1$x, dat.nb1$y), diag = F, upper = T) )
distmat <- ifelse(distmat > 4 | distmat == 0, NA, distmat)

# --- creat a wrapped distance matrix
tdmat <- matrix(NA, n, n)
for(i in 1:n){
  for(j in 1:n){
    tdmat[i,j] <- toroid.dist(dat.nb1$x[i], dat.nb1$y[i],
                              dat.nb1$x[j], dat.nb1$y[j],
                              max(dat.nb1$x), max(dat.nb1$y))
  }
}
tdmat <- ifelse(tdmat > 4 | tdmat == 0, NA, tdmat)
# --- calculate crowding terms standardized for each timestep
lout <- list()
for(i in 1:(max(dat.nb$census) - 1)) {
  hood <- nb.fn( dat$size_t[dat$census==i], distmat )
  lout[[i]] = scale.fn(log(hood))
}
crowd <- do.call(c, lout)
# --- organize
df <- data.frame(y = y*1e6, x = x*1e6, type = as.factor(type), nb = crowd, t = dat.nb$census)
df$x[df$t == 4] = df$x[df$t == 4]*5 # correct for the long sampling interval

# --- fit glm's
m.nonspat <- glmer(y ~ (1|type) + log(x) + offset(log(x)),
           data = df, 
           family = Gamma(link = log))
m <- glmer(y ~ (1|type) + log(x) + (0+nb|type) + offset(log(x)),
           data = df, 
           family = Gamma(link = log)) 

ahat <- exp(predict(m.nonspat))
mhat <- exp(predict(m))
# --- compare prediction error
sqrt(mean((df$y - mhat)^2))*1e-6
sqrt(mean((df$y - ahat)^2))*1e-6

# --- extract non-spatial coefficients
coefs.n <- coefficients(m.nonspat)$type
names(coefs.n) <- c("Int", "slope")
coefs.n$type <- rownames(coefs.n)
# --- extract spatial coefficients
coefs.s <- coefficients(m)$type
names(coefs.s) <- c("nb", "Int.nb", "slope.nb")
coefs.s$type <- rownames(coefs.s)

# === IBM simulation
meta.df <- dat.nb1 %>% 
  left_join(coefs.n, by = c("subsppcyt" = "type")) %>% 
  left_join(coefs.s, by = c("subsppcyt" = "type")) %>%
  rename(n0 = size_t) %>%
  rename(type = subsppcyt) %>% 
  mutate(type = as.factor(type)) %>%
  select(Int, slope, Int.nb, slope.nb, nb, type, n0)

# === naive
n <- nrow(meta.df)
T <- 20
mat <- matrix(NA, n, T)
mat[,1] <- log(meta.df$n0)
for(i in 2:T) {
  mat[,i] <- meta.df$Int + meta.df$slope*mat[,i-1] + mat[,i-1]
}
boxplot(exp(-meta.df$Int/meta.df$slope)*1e-6 ~ meta.df$type)
points(exp(mat[,T])*1e-6 ~ jitter(as.numeric(meta.df$type)), col = rgb(.5,0,0,.1), pch = 19)

# === non-spatial vs spatial
n <- nrow(meta.df)
T <- 20
mat.n <- matrix(NA, n, T)
mat.n[,1] <- log(meta.df$n0)
mat.s <- matrix(NA, n, T)
mat.s[,1] <- log(meta.df$n0)

for(i in 2:T) {
  hood <- ( log( nb.fn(exp(mat.s[,i-1]), (.0010*distmat) ) ) )
  # NDD removed
  mat.n[,i] <- meta.df$Int.nb + 
    meta.df$slope.nb*mat.n[,i-1] + mat.n[,i-1] 
  # NDD present
  mat.s[,i] <- meta.df$Int.nb + 
    meta.df$slope.nb*mat.s[,i-1] + mat.s[,i-1] + meta.df$nb*hood 
}

plot(exp(mat[,T])*1e-6 ~ jitter(as.numeric(meta.df$type)), col = rgb(.5,1,1,.1), pch = 19)

# observed data
boxplot(df$y[df$t==5]*1e-6 ~ df$type[df$t==5], col = rgb(0,0,0,.1), 
        ylab = "Crown volume", xlab = "Subsppcyt")
# predicted NDD removed
t <- 10
points(exp(mat.n[,t])*1e-6 ~ jitter(as.numeric(meta.df$type)),
     col = rgb(.5,0,0,.25), pch = 19)
# predicted NDD present
plot(exp(mat.s[,t])*1e-6 ~ jitter(as.numeric(meta.df$type)), col = rgb(.5,.5,1,.25), pch = 19)



