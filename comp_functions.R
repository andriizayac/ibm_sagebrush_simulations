a <- seq (1, 38, l = 476)
b <- seq(from = 1, to = 100, l = 100)
plot(dist_vec, (dist_vec*exp(-0.2*dist_vec/0.4)))
sort(dist_vec)

plot(dist_vec, dist_vec*(1-exp(-0.2*dist_vec))^(1/1-dist_vec*0.9))


plot(dist_vec, 1/log10(dist_vec^5))
sizemat2 <- seq(5, 150, l = 476)

plot(dist_vec, (exp(dist_vec*0.7)/(1+exp(dist_vec))))
plot(dist_vec, ((exp(-distance_strength*dist_vec)/sizemat2)*competition_strength))
     
plot(dist_vec, ((exp(-3*dist_vec)/sizemat2)*10))
