##read files
#3gardsurv is downloaded from the paper attributes 
#orchardsurv <- filter(read.csv("3gardsurv_2015.csv", sep = ",", header = TRUE), garden == "Orchard")

#orchardcm.csv created from datafixed.xlsx which is the original data with a couple modifications of Subspp.Prov.ID
orchard <- read.csv("orchardcm.csv", header = TRUE, sep = ",")
#orchardlive <- filter(orchard, survivorship != "Dead")

#create date variables POSIXct is used because full_join() doesn't like PPOSIXlt
orchard$datestart <- as.POSIXct("2010-05-01", format = "%Y-%m-%d")
orchard$datelast <- as.POSIXct("2015-04-01", format = "%Y-%m-%d")

##transform/rename measurements to numeric variables
# height at outplant is reasonably selected as 2 in
orchard$height1_in <- as.numeric(2)
# width at outplnat is reasonably selected as 1 in
orchard$width1_in <- as.numeric(1)
#renaming for conveniences of typing in variables
names(orchard)[names(orchard) == "width..in."] <- "width_in"
names(orchard)[names(orchard) == "height..in."] <- "height_in"

# if measurements are imported as integers - transform to numerical values
orchard$height_in <- as.numeric(orchard$height_in)
orchard$width_in <- as.numeric(orchard$width_in)

#calculate crown area for the first and last measurements
#orchard$areastart <- (orchard$height1/2) * (orchard$width1/2) *pi
#orchard$arealast <- (orchard$height_in/2) * (orchard$width_in/2) *pi

# create test copy of the data frame: 
#test1 <- orchard


# inches to metric
orchard$width_cm <- orchard$width_in *2.54
orchard$height_cm <- orchard$height_in *2.54
#rename measurment variables
#names(test1)[names(test1) == "width_in"] <- "width"
#names(test1)[names(test1) == "width_in"] <- "width"
# inches to metric
orchard$width1_cm <- orchard$width1_in *2.54
orchard$height1_cm <- orchard$height1_in *2.54

# calc crouwn volume a/2 * b2/2 * c/2 *pi
## * 1e-06 transforms to cubic
#names(test1)[names(test1) == "arealast"] <- "volumelast"
orchard$volumestart <- ((orchard$height1_cm/2) * (orchard$width1_cm/2)^2 * pi * (4/3))/1e+06
orchard$volumelast <- ((orchard$height_cm/2) * (orchard$width_cm/2)^2 * pi * (4/3))/1e+06
# calc crown area
orchard$areastart <- (orchard$width1_cm/200)^2 * pi #pi*r^2 = pi*(d/2*100)^2: since measurements are in cm, devide by 100 to get sq m
orchard$arealast <- (orchard$width_cm/200)^2 * pi #pi*r^2 = pi*(d/2*100)^2: since measurements are in cm, devide by 100 to get sq m

# add empty rows to the data frame to match the matrix/map size of tag id's
orchard[nrow(orchard)+8, ] <- NA

# add a column merging subspp with cytotype
#orchard$subsppcyt <- paste(orchard$Subspecies, orchard$Cytotype, sep = "")
#orchard$subsppcyt[orchard$subsppcyt == "NANA"] <- NA

##### data frame "orchard" is produced after running all this code ###########
str(orchard)
write.csv(orchard, "orchard.csv")
