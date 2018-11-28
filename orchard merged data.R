#import source-population map as a matrix
#matrix1 <- as.matrix(read.csv("orchardcmmap1.csv", header = FALSE, sep = ","))
tags <- as.matrix(read.csv("orchardcmmap2.csv", sep = ",", header = FALSE))

#create x and y variables for the grid
x <- rep(1:28, times = 17)
y <- sort(rep(seq(from = 0, to = 24, by = 1.5), times = 28))

#df <- matrix(x, nrow = 17, ncol = 28, byrow = TRUE)
#df1 <- matrix(sort(rep(seq(from = 0, to = 24, by = 1.5), times = 28)), nrow = 17, ncol = 28, byrow = TRUE)

#convert the matrices by row - "t(matrix)" and combine the vectors into a data frame

# merge x and y into a data frame and add the tag ID vector
map <- data.frame(x, y)
#names(map)[names(map) == "tagID"] <- "Sample.tag."
map$Sample.tag. <- as.vector(t(tags))



# join two dataframes by source ID

# extra rows for the spaces in the vector (t(map)) to match the length: orchard1[nrow(orchard1)+1, ] <- NA
orchard1 <- merge(orchard, map, by.x = "Sample.tag.", by.y = "Sample.tag.")

#remove duplicated rows to return to the original data frame size
orchard1 <- orchard1[-c(477:532),]

######### data frame "orchard1" is used from here on as a merged copy of two files ##############
# random plots
ggplot(orchard1, aes(x = x, y=y, colour = subsppcyt)) + 
 geom_point(aes(size = volumelast)) + ggtitle("Orchard data")
  #geom_point(aes(size = areastart, colour = "start"))

rm(x,y)

write.csv(orchard1, "orchard1.csv")
