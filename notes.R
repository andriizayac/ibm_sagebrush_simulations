library(manipulate)
manipulate({
  lambda=seq(0,5,l=100)
  prior=dgamma(lambda,a,b)
  post=dgamma(lambda,a+sum(y),b+length(y))
  tibble(lambda,prior,post)%>%
    ggplot(aes(lambda,prior))+geom_line()+
    geom_line(aes(post),colour="red")+
    ylab("density")+xlab(expression(lambda))
},
a=slider(1,10,step=1),
b=slider(1,10,step=1),
y=slider(1,10,step=1)
)

manipulate({
  lambda=seq(0,5,l=100)
  prior=dnorm(lambda,a,b)
  post=dnorm(lambda,a+sum(y),b+length(y))
  tibble(lambda,prior,post)%>%
    ggplot(aes(lambda,prior))+geom_line()+
    geom_line(aes(post),colour="red")+
    ylab("density")+xlab(expression(lambda))
},
a=slider(1,10,step=1),
b=slider(1,10,step=1),
y=slider(1,10,step=1)
)


entropy <- function(p) -sum(p * log(p))
n <- 1e3
p <- runif(n)
q <- 1 - p
z <- matrix(c(p, q), ncol = 2)
ent <- apply(z, 1, entropy)
qplot(q, ent, geom = "line") + 
  ylab("Entropy") + xlab("Probability") +
  ggtitle("Entropy For Two Events")
