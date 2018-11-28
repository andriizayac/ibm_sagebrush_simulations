library(plyr)
library(dplyr)
library(tidyr)

# subset data from the main data set (orchard1, generated earlier)
spatial <- orchard1[, c("Sample.tag.", "volumelast", "Subspp.Prov.ID", "subsppcyt", "x", "y")]

# add x and y coordinates for each plant, and generate a grid
#x <- 1:28
#y <- 1:17
#a <- expand.grid(x, y)
#spatial[, c("x", "y")] <- a[, c("Var1", "Var2")]
#rm(a)
spatial$y <- spatial$y/1.5 + 1 # transform from actual grid (1.5m on y) to 1 by 1 grid
# distance function
dist.fx <- function(x1, x2, y1, y2) {sqrt((x2-x1)^2 + (y2-y1)^2)}

# generate two lists to find the ID of neighbors for each plant
# list distmatbyID contans 476 matrices with with a zero at each consecutive focal plant
# list neighborlist subsets the tag.id's of neighbors. limits neighbors to the sed disctance range
#distmat <- matrix(NA, nrow = 17, ncol = 28)
distmatbyID <- list()
neighborlist <- list()
for (j in 1:476){
  distmatbyID[[j]] <- matrix(dist.fx(spatial$x[j], spatial$x, spatial$y[j], spatial$y), nrow = 17, ncol = 28, byrow = TRUE)
  neighborlist[[j]] <- as.numeric(which((distmatbyID[[j]] < 2) & (distmatbyID[[j]] > 0)))
  
}

#data.frame(matrix(unlist(b[1:10]), nrow = length(b[1:10]), ncol = 8, byrow = TRUE))

# drop x and y variables
spattable <- spatial[,!names(spatial) %in% c("x", "y")]
spattable$subsppcyt <- as.factor(spattable$subsppcyt)


# generate/calculate additional columns
spattable$numdiploids <- rep(NA, times = 476)
spattable$numtetraploids <- rep(NA, times = 476)
neighborgen <- list()
spattable$numT2n <- rep(NA, times = 476)
spattable$numT4n <- rep(NA, times = 476)
spattable$numW4n <- rep(NA, times = 476)
spattable$numV4n <- rep(NA, times = 476)
spattable$numV2n <- rep(NA, times = 476)

for (i in 1:476){
  c <- unlist(neighborlist[i])
  spattable$numdiploids[i] <- sum(spattable$Cytotype[c] == "2n", na.rm = TRUE)
  spattable$numtetraploids[i] <- sum(spattable$Cytotype[c] == "4n", na.rm = TRUE)
  spattable$numT2n[i] <- sum(spattable$subsppcyt[c] == "T2n", na.rm = TRUE)
  spattable$numT4n[i] <- sum(spattable$subsppcyt[c] == "T4n", na.rm = TRUE)
  spattable$numW4n[i] <- sum(spattable$subsppcyt[c] == "W4n", na.rm = TRUE)
  spattable$numV4n[i] <- sum(spattable$subsppcyt[c] == "V4n", na.rm = TRUE)
  spattable$numV2n[i] <- sum(spattable$subsppcyt[c] == "V2n", na.rm = TRUE)
  
  
  
  neighborgen[[i]] <- paste(unique(spattable$subsppcyt[c]), collapse = ' ')
  rm(c)
}

rm(i,j,x,y)
# convert newly generated columns into numeric type
spattable[,5:11] <- sapply(spattable[5:11], as.numeric)

#random check-up
#spattable$Cytotype[c(unlist(b[469]))]
#numdiploids[469]
###concatinates character values of a list into a string: 
###paste(neighborgen[[1]], sep = "", collapse = ' ')

#spattable <- data.frame(spattable, numdiploids, numtetraploids, numT2n, numT4n, numW4n, numV4n, numV2n)

reg <- lm(volumelast ~ numdiploids + 
            numtetraploids + numT2n + 
            numT4n + numW4n + numV4n + 
            numV2n, data = spattable)
summary(reg)


 

# interaction matrix on growth
reg <- glm(volumelast ~ subsppcyt*numT2n + 
                  subsppcyt*numT4n + subsppcyt*numW4n + subsppcyt*numV4n + 
                  subsppcyt*numV2n, data = spattable)
summary(reg)

#interaction matrix on survival
# survival variable (0 - dead, 1 - alive)
spattable$surv<-ifelse(is.na(spattable$volumelast)==T,0,1)


binomreg <- glm(surv ~ subsppcyt*numT2n + 
            subsppcyt*numT4n + subsppcyt*numW4n + subsppcyt*numV4n + 
            subsppcyt*numV2n, data = spattable, family="binomial")
summary(binomreg)



summary(reg)
reg1 <- predict(reg)
reg1 <- data.frame(x = c(1:length(reg1)), y = reg1)
ggplot(reg1, aes(x,y)) + geom_point() + geom_smooth(method = "lm")

