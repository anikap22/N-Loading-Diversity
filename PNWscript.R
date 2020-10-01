#####################################
##
## PNW Soils project -- Anika Petach
## Contents: Mixed Models
##
##
## Q: At what level of nitrogen limitation (Soil N) does productivity shift from multispecies domination to single species domination?
######################################
library(bbmle)
library(rjags)
library(lattice)
library(foreign)
library(ggplot2)
library(plyr)
library(rworldmap)
require(GGally)
require(ggmap)
require(lmerTest)
require(LMERConvenienceFunctions)
require(arm)
require(coefplot)
library(merTools)
library(MCMCglmm)
library(dplyr)
library(lme4)
library(RLRsim)
source("/Users/Anika/Documents/R/StatsClass/PNWproject/glmm_funs.R.R")
library(sjPlot)
library(sjmisc)
library(mcmcplots)

setwd("/Users/Anika/Documents/R/StatsClass/PNWproject") #set working directory
getwd()                           #print the working directory

# Read in the data & format dataframe -------------------------------------
#############################
productivity <- read.dbf("NACP_TERRA_PNW_forest_biomass_productivity.dbf", as.is = FALSE)
leaf <- read.dbf("NACP_TERRA_PNW_leaf_trait.dbf", as.is = FALSE)
soil <- read.dbf("NACP_TERRA_PNW_soil.dbf", as.is = FALSE)

# for each plot there are many different plants in the leaf file
# each plot has one soil and productivity measure

# Look at the types of variables
head(productivity)
head(leaf)
head(soil)
length(unique(productivity$PLOT_ID))   #265 plots
nm <- c("ASA","BULK_DENSI")
plot.all <- as.data.frame(NULL)
plot.all[nm] <- lapply(nm, function(x){productivity[[x]][match(productivity$PLOT_ID, soil$PLOT_ID)]})
  
length(unique(leaf$SPECIES))  #35 leaf species
xtabs(~PLOT_ID + MAP, data=productivity)

# Convert files to csv to look at them more easily
write.csv(productivity,"productivity.csv")
write.csv(leaf,"leaf.csv")
write.csv(soil,"soil.csv")
productivity <- read.csv("productivity.csv", header=TRUE)
leaf <- read.csv("leaf.csv", header=TRUE)
soil <- read.csv("soil.csv", header=TRUE)

summary(productivity)
summary(leaf) #leaf carbon, sla_psa, leaf nitro, leaf CN, leaf life
summary(soil)  #nitrogen c, bulk density, MAT, 

#productivity assign NA to none
productivity[] <- lapply(productivity, function(x){replace(x, x == "none", NA)})  # Recode "none" as Na
soil[] <- lapply(soil, function(x){replace(x, x == -9999, NA)})  # Recode -9999 as Na
leaf[] <- lapply(leaf, function(x){replace(x, x == -9999, NA)})  # Recode -9999 as Na
productivity[] <- lapply(productivity, function(x){replace(x, x == -9999, NA)})  # Recode -9999 as Na

# Merge productivity and soil
length(intersect(productivity$PLOT_ID, soil$PLOT_ID)) #211 plot IDs in common between productivity and soil
vector.list <- list(productivity$PLOT_ID, soil$PLOT_ID)
length(intersect(productivity$PLOT_ID, soil$PLOT_ID))
commonPlots <- Reduce("intersect", vector.list)   #list of all common plots between prod and soil
length(commonPlots)
productivity.sub <- productivity[productivity$PLOT_ID %in% commonPlots,]
# no middle on 86, 201, 203, 204, 208, 213, 222, 226, 234, 235, 236, 237, 242, 243, 244, 245, 246, 247, 248, 253, 255, 256, 260, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 311, 312, 313(use top)
nomiddle <- list(86, 201, 203, 204, 208, 213, 222, 226, 234, 235, 236, 237, 242, 243, 244, 245, 246, 247, 248, 253, 255, 256, 260, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 311, 312, 313)
soil.sub <- soil[(soil$SOIL_LAYER == "middle"|soil$PLOT_ID %in% nomiddle),]
soil.sub <- soil.sub[soil.sub$PLOT_ID %in% commonPlots,]
merged.soil.prod <- merge(productivity.sub, soil.sub, by='PLOT_ID')

# Plot numbers used
plots <- unique(merged.soil.prod$PLOT_ID)

##Find number of dominant species (productivity df)
#make matrix of productivity$SPP_O1_ABB, productivity$SPP_02_ABB, productivity$SPP_03_ABB, productivity$SPP_04_ABB
prod.matrix <- as.matrix(cbind(merged.soil.prod$SPP_O1_ABB, merged.soil.prod$SPP_O2_ABB, merged.soil.prod$SPP_O3_ABB, merged.soil.prod$SPP_O4_ABB))
# RowSums(!is.na(matrix))
NumDom <- rowSums(!is.na(prod.matrix))
NumDom <- cbind(merged.soil.prod$PLOT_ID, NumDom)
merged.soil.prod <- cbind(merged.soil.prod, NumDom)
colnames(NumDom) <- c("PLOT_ID","NumDom")
NumDom <- data.frame(NumDom)
head(NumDom)
table(NumDom$NumDom)

# Subset from the leaf data
leaf.sub <- leaf[leaf$PLOT_ID %in% plots,]
merged.leaf <- join(leaf.sub, merged.soil.prod, by="PLOT_ID", match="all")
write.csv(merged.leaf, "mergeleaf.csv")
merged.leaf.sub <- merged.leaf[,c(2,5:17,20,21,22:27,39:41,43:60,77:85)]
data1 <- merged.leaf.sub
head(data1)
data1[] <- lapply(data1, function(x){replace(x, x == -9999, NA)})  # Recode -9999 as Na

# Extract key variables
MAP <- data1$MAP
MAT <- data1$MAT_C
SoilN <- data1$NITROGEN_C

################################################# End data formatting ################
# Histograms --------------------------------------------------------------

#2. Histograms of distributions
x <- seq(0,25,0.1)
par(mfrow=c(2,3))
hist(data1$MAT_C, prob=TRUE, xlab="MAT", main="Hist MAT")   #normal
mean.mat <- mean(data1$MAT_C, na.rm=TRUE)
sd.mat <- sd(data1$MAT_C, na.rm=TRUE)
curve(dnorm(x, mean=mean.mat, sd=sd.mat), add = TRUE, col = 'blue') 
hist(data1$NumDom, prob=TRUE, xlab="# Dom Species", main="Hist of Dom Species")
hist(data1$ELEVATION.x, prob=TRUE, breaks=20, xlab="Elev", main="Hist of Elev")   #normal (gamma?)
mean.elev <- mean(data1$ELEVATION.x, na.rm=TRUE)
sd.elev <- sd(data1$ELEVATION.x, na.rm=TRUE)
var.elev <- sd.elev^2
curve(dnorm(x, mean=mean.elev, sd=sd.elev), add = TRUE, col = 'red') 
curve(dgamma(x, shape=mean.elev^2/var.elev, rate=mean.elev/var.elev), add = TRUE, col = 'blue') 
hist(data1$MAP.x, prob=TRUE, xlab="MAP", main="Hist of MAP")
map.mean <- mean(data1$MAP.x, na.rm=TRUE)
map.var <- var(data1$MAP.x, na.rm=TRUE)
curve(dexp(x, rate=1/map.mean), add = TRUE, col = 'blue') 
hist(data1$BULK_DENSI, prob=TRUE, breaks=20, xlab="Bulk Density", main="Hist of Bulk Den")  #gamma
dens.mean = mean(data1$BULK_DENSI)
dens.var = var(data1$BULK_DENSI)
curve(dgamma(x, shape=dens.mean^2/dens.var, rate=dens.mean/dens.var), add = TRUE, col = 'blue') 
hist(merged.soil.prod$NITROGEN_C, prob=TRUE, breaks=20, xlab="Soil N", main="Hist Soil N")   #gamma
nitro.mean <- mean(merged.soil.prod$NITROGEN_C, na.rm=TRUE)
nitro.var <- var(merged.soil.prod$NITROGEN_C, na.rm=TRUE)
curve(dgamma(x, shape=nitro.mean^2/nitro.var, rate=nitro.mean/nitro.var), add = TRUE, col = 'blue') 
#histograms part 2
par(mfrow=c(2,3))
hist(merged.soil.prod$PH_OF_SOIL, prob=TRUE, xlab="pH", main="Hist of pH")  #normal?
ph.mean <- mean(data1$PH_OF_SOIL, na.rm=TRUE)
ph.sd <- sd(data1$PH_OF_SOIL, na.rm=TRUE)
curve(dnorm(x, mean=ph.mean, sd=ph.sd), add = TRUE, col = 'red') 
hist(merged.soil.prod$CARBON_CON, prob=TRUE, xlab="Soil C", main="Hist of Soil C")  #gamma
carb.mean = mean(merged.soil.prod$CARBON_CON, na.rm=TRUE)
carb.var = var(merged.soil.prod$CARBON_CON, na.rm=TRUE)
curve(dgamma(x, shape=carb.mean^2/carb.var, rate=carb.mean/carb.var), add = TRUE, col = 'blue') 
hist(leaf$LAI_O, prob=TRUE, xlab="LAI", main="Hist of Canopy LAI")
hist(merged.leaf$SPP_O1_BAS, prob=TRUE, xlab="% BA in Dom Speices", main="Hist of BA in Dom Species")
spp.mean = mean(merged.leaf$SPP_O1_BAS, na.rm=TRUE)
spp.log.mean = mean(log(merged.leaf$SPP_O1_BAS+1),na.rm=TRUE)
spp.var = var(merged.leaf$SPP_O1_BAS, na.rm=TRUE)
spp.log.sd = sd(log(merged.leaf$SPP_O1_BAS+1),na.rm=TRUE)
curve(dgamma(x, shape=spp.mean^2/spp.var, rate=spp.mean/spp.var), add = TRUE, col = 'blue') 
curve(dlnorm(x, meanlog=spp.log.mean, sdlog=spp.log.sd), add=TRUE, col='red')
hist(leaf$LEAF_CN, prob=TRUE, xlab="Leaf C:N", main="Hist of Leaf C:N")  #normal?
cn.mean <- mean(data1$LEAF_CN, na.rm=TRUE)
cn.sd <- sd(data1$LEAF_CN, na.rm=TRUE)
curve(dnorm(x, mean=cn.mean, sd=cn.sd), add = TRUE, col = 'red') 
curve(dgamma(x, shape=cn.mean^2/cn.sd^2, rate=cn.mean/cn.sd^2), add = TRUE, col = 'blue') 
hist(leaf$LEAF_NITRO, prob=TRUE,breaks=15, xlab="Leaf N", main="Hist of Leaf N")  #gamma? lognormal?
nitro.mean <- mean(data1$LEAF_NITRO, na.rm=TRUE)
nitro.var <- var(data1$LEAF_NITRO, na.rm=TRUE)
curve(dgamma(x, shape=nitro.mean^2/nitro.var, rate=nitro.mean/nitro.var), add = TRUE, col = 'blue') 
#Histogram for response variable (is it lognormally distributed?)
par(mfrow=c(1,1))
lin.spp <- log(merged.leaf$SPP_O1_BAS+1)
lin.spp.mean <- mean(lin.spp, na.rm=TRUE)
lin.spp.sd <- sd(lin.spp, na.rm=TRUE)
hist(lin.spp, prob=TRUE, xlab="Log(% BA in Dominant Speices)", main="Hist of %BA in Dominant Species")
curve(dgamma(x, shape=lin.spp.mean^2/lin.spp.sd^2, rate=lin.spp.mean/lin.spp.sd^2), add = TRUE, col = 'blue') 
lin.spp.frac <- lin.spp/100
lin.spp.frac.mean <- mean(lin.spp.frac, na.rm=TRUE)
lin.spp.frac.var <- var(lin.spp.frac, na.rm=TRUE)
shape1.spp <- (lin.spp.frac.mean^2-lin.spp.frac.mean^3-lin.spp.frac.mean*lin.spp.frac.var)/lin.spp.frac.var
shape2.spp <- (lin.spp.frac.mean-2*lin.spp.frac.mean^2+lin.spp.frac.mean^3-lin.spp.frac.var+lin.spp.frac.mean*lin.spp.frac.var)/lin.spp.frac.var
hist(lin.spp.frac, prob=TRUE, xlab="Fraction of BA in Dominant Species", main="Hist of Fraction BA in Dominant Species")
curve(dbeta(x, shape1=shape1.spp, shape2=shape2.spp), add=TRUE, col="green")

#Final histograms to show (lin.spp, lat, mat, soil N, soil pH, bulk den)
par(mfrow=c(2,3))
lin.spp.mean <- mean(lin.spp, na.rm=TRUE)
lin.spp.sd <- sd(lin.spp, na.rm=TRUE)
hist(lin.spp, prob=TRUE, xlab="Log(% BA in Dominant Speices)", 
     main="Log(%BA) in Dominant Species")
curve(dgamma(x, shape=lin.spp.mean^2/lin.spp.sd^2, 
             rate=lin.spp.mean/lin.spp.sd^2), add = TRUE, col = 'blue') 
curve(dnorm(x,mean=lin.spp.mean, sd=lin.spp.sd), add = TRUE, col = 'red')
hist(newdata$LAT, prob=TRUE, xlab="Latitude", main="Latitude", breaks=20)
lat.mean <- mean(newdata$LAT, na.rm=TRUE)
lat.var <- var(newdata$LAT, na.rm=TRUE)
curve(dgamma(x, shape=lat.mean^2/lat.var, rate=lat.mean/lat.var), 
      add = TRUE, col = 'blue') 
hist(data1$MAT_C, prob=TRUE, xlab="Temperature", main="Mean Annual Temperature")   #normal
mean.mat <- mean(data1$MAT_C, na.rm=TRUE)
sd.mat <- sd(data1$MAT_C, na.rm=TRUE)
curve(dnorm(x, mean=mean.mat, sd=sd.mat), add = TRUE, col = 'red') 
hist(merged.soil.prod$NITROGEN_C, prob=TRUE, 
     breaks=20, xlab=expression('Soil N (gN/m'^'2)'), main=expression('Soil N (gN/m'^'2)') )  #gamma
nitro.mean <- mean(merged.soil.prod$NITROGEN_C, na.rm=TRUE)
nitro.var <- var(merged.soil.prod$NITROGEN_C, na.rm=TRUE)
curve(dgamma(x, shape=nitro.mean^2/nitro.var, rate=nitro.mean/nitro.var), add = TRUE, col = 'blue') 
hist(merged.soil.prod$PH_OF_SOIL, prob=TRUE, xlab="pH", main="Soil pH")  #normal?
ph.mean <- mean(newdata$SOILPH, na.rm=TRUE)
ph.sd <- sd(newdata$SOILPH, na.rm=TRUE)
curve(dnorm(x, mean=ph.mean, sd=ph.sd), add = TRUE, col = 'red') 
hist(data1$BULK_DENSI, prob=TRUE, breaks=20, xlab="Bulk Density", 
     main=expression('Bulk Denisty (kg/m'^'2)'))  #gamma
dens.mean = mean(data1$BULK_DENSI)
dens.var = var(data1$BULK_DENSI)
curve(dgamma(x, shape=dens.mean^2/dens.var, rate=dens.mean/dens.var), 
      add = TRUE, col = 'blue') 

#Exploratory plots
data1.small <- data1[,c(5,7,16,20,21,25,28)]
plot(data1.small)

#Exploratory plots for variables as of 11/17
data2.small <- data1[,c(2,5:8,11,12,20,21,23:28,35:42,43,45,46)]
dotplot(SPP_O1_BAS~NITROGEN_C | SPP_O1_ABB, data=data2.small)
barchart(SPP_O1_BAS~NITROGEN_C | SPP_O1_ABB, data=data2.small, xlab="Soil N (gN/m2)", ylab="% BA in Primary Species", main="%BA from primary species by Soil N")

mycols <- c("YEAR","PLOT_ID","LATITUDE","LAI_O","LEAR_CAR_1","LEAF_NITRO","LEAF_CN","ASA","SPP_O1_BAS","AG_BIOMASS","AG_BIOMA_1","AG_PROD_TR")
mycols <- which(names(data1) %in% mycols)
data1.small2 <- data1[,mycols]
plot(data1.small2)

# Exploratory Plots -------------------------------------------------------

#scatter plots
par(mfrow=c(2,3))
plot(data1$NITROGEN_C, data1$LEAF_NITRO, xlab="Soil N",ylab="Leaf N")
plot(data1$NITROGEN_C, data1$LEAF_CN, xlab="Soil N",ylab="Leaf C:N")
plot(data1$NumDom, data1$LEAF_CN, xlab="Number Dominant Species",ylab="Leaf C:N")
plot(data1$LAI_O, data1$LEAF_CN, xlab="LAI",ylab="Leaf C:N")
plot(data1$HEIGHTC_m, data1$LEAF_CN, xlab="Height",ylab="Leaf C:N")
plot(data1$LEAF_NITRO, data1$LEAF_CN, xlab="Leaf N", ylab="Leaf C:N")

par(mfrow=c(2,3))
plot(data1$LEAF_CARBO, data1$LEAF_CN, xlab="Leaf C",ylab="Leaf C:N")
plot(data1$LEAF_CARBO, data1$LEAF_NITRO, xlab="Leaf C",ylab="Leaf N")
plot(data1$SPP_O1_ABB, data1$LEAF_CN, xlab="Dominant Species",ylab="Leaf C:N")   #Dominant name vs. nitrogen limitation
plot(data1$SPP_O1_BAS, data1$LEAF_CN, xlab="%BA in primary species",ylab="Leaf C:N")
plot(data1$SPP_O1_BAS, data1$ASA, xlab="%BA in primary species",ylab="Stand age")
plot(data1$NumDom, data1$ASA, xlab="Number of dominant species",ylab="Stand age")
plot(data1$NumDom, data1$LEAF_NITRO)
par(mfrow=c(1,1))
boxplot(data1$SPP_O1_BAS, data1$NITROGEN_C)
barchart(SPP_O1_BAS~NITROGEN_C | SPP_O1_ABB, data=data1)   # Fewer Dominant species when leaf nitrogen lower within species
barchart(NumDom~LEAF_CN | GENUS, data=data1)   # Different genuses have one Dom species at low CN
barchart(SPP_O1_BAS~LEAF_NITRO | GENUS, data=data1)
bwplot(SPP_O1_BAS~LEAF_CN, data=data1)
bwplot(NumDom~LEAF_NITRO, data=data1)
dotplot(SPP_O1_BAS~LEAF_NITRO | GENUS, data=data1)
barchart(NumDom~LEAF_NITRO | SPP_O1_ABB, data=data1)

table(data1$NumDom)
table(merged.soil.prod$NumDom)
barchart(SPP_O1_BAS~NITROGEN_C | ECOREGION.x, data=merged.soil.prod,xlab="Soil N (gN/m2)",ylab="%BA in primary species",main="%BA of Primary Species vs. Soil N by Ecoregion")
barchart(NumDom~)

#xy plots
xyplot(data1$LATITUDE~data1$LEAF_CN|data1$GENUS, type = "p", auto.key = TRUE)
xyplot(NumDom~LEAF_CN|GENUS, data=data1, type = "p", auto.key = TRUE)
xyplot(SPP_O1_BAS~LEAF_CN|ECOREGION, data=data1, type="p")
xyplot(SPP_O1_BAS~NITROGEN_C|ECOREGION, data=data1, type="p",xlab="Soil N (gN/m2)",ylab="%BA in primary species",main="%BA of Primary Species by Soil N")
xyplot(SPP_O1_BAS~NITROGEN_C|SPP_O1_ABB, data=data1, type="p",xlab="Soil N (gN/m2)",ylab="%BA in primary species",main="%BA of Primary Species by Soil N")

xyplot(SPP_O1_BAS~NITROGEN_C|ECOREGION, data=data1, type="p", xlab="Soil N (gN/m2)",ylab="%BA in primary species")

xyplot(SPP_O1_BAS~NITROGEN_C|factor(FIXER), data=merged.leaf, type="p")
xyplot(SPP_O1_BAS~NITROGEN_C,data=merged.leaf,groups=factor(FIXER,labels=c("0","1")),
       pch=20,auto.key=list(columns=2),type=c("p","g"),xlab="Soil N (gN/m2)", ylab="%BA by Dominant Species")

xyplot(SPP_O1_BAS~NITROGEN_C|factor(ASA_R), data=merged.leaf, type="p", xlab="Soil N (gN/m2)",ylab="%BA in primary species")
xyplot(SPP_O1_BAS~NITROGEN_C|factor(ASA_R1), data=merged.leaf, type="p", xlab="Soil N (gN/m2)",ylab="%BA in primary species")

# Exploratory ggplots -----------------------------------------------------

hp<-qplot(x=merged.leaf$ASA, fill=..count.., geom="histogram") 
hp
# Sequential color scheme
hp+scale_fill_gradient(low="blue", high="red")

# Early plots with ggplot2
qplot(lin.spp, c.SOILN, data=newdata, xlab="Soil N", ylab="%BA", main="%BA in Dominant Species vs. Soil N")

qplot(as.factor(FIXER), lin.spp, facets = . ~ ECOREGION, 
      colour = ECOREGION, geom = "boxplot", data = newdata,
      xlab="Fixer Status", ylab="log(%BA)")

# ggplot with regression
mod1 = lm(SPP_O1_BAS ~ NITROGEN_C, data = merged.soil.prod)
r2 <- summary(mod1)$R.squared
r2.adj <- summary(mod1)$adj.r.squared
int1 <- signif(coef(mod1)[1], digits = 2)
slope1 <- signif(coef(mod1)[2], digits = 2)
textlab <- paste("y = ",slope1,"x + ",int1, sep="")
p1<- ggplot(merged.soil.prod, aes(x=SPP_O1_BAS, y=NITROGEN_C)) 
p1 + geom_point(size=1) + geom_smooth(method=lm, se=TRUE) +  # Hollow circles, Add linear regression line 
  annotate("text", x = 75, y = 800, label = textlab, color="black", size = 5, parse=FALSE) +
  labs(x="%BA from Primary Species", y="Soil Nitrogen")+ggtitle("Soil Nitrogen vs. %BA Dominant")

# ggplot % in Dom Species vs. Lat by Species
p <- ggplot(merged.soil.prod, aes(SPP_O1_BAS, LATITUDE.x))
p+geom_point(aes(color = factor(SPP_O1_ABB)))+labs(x="%BA in Dom Species", y="Latitude")+ggtitle("Latitude vs. %BAs in Dominant Species by Species")

# ggplot with foliar biomass vs. leaf cn by elevation
library(dplyr)
data1 %>%
  group_by(GENUS) %>%
  summarize(mean(data1$LEAF_CN,na.rm=T), sd(data1$LEAF_CN,na.rm=T), mean(data1$AG_BIOMA_1,na.rm=T), sd(data1$AG_BIOMA_1,na.rm=T), cor(data1$LEAF_CN,data1$AG_BIOMA_1))
p = ggplot(data1, aes(data1$LEAF_CN, data1$AG_BIOMA_1, color=data1$ELEVATION)) + geom_point()
p = p + geom_smooth(method = lm, se = FALSE)
p = p + facet_wrap(~data1$GENUS) + labs(x="Leaf C:N", y="Foliar Biomass", title="Foliar Biomass vs. Leaf C:N by Genus")
p = p + labs(color="Elevation")
p

p1 = ggplot(data1, aes(data1$NITROGEN_C, data1$SPP_O1_BAS, color=data1$LATITUDE)) + geom_point()
p1 = p1 + geom_smooth(method = lm, se = FALSE)
p1 = p1 + facet_wrap(~GENUS) + labs(x="Soil N", y="%BA", title="%BA Dom Species vs. Soil N by Genus")
p1 = p1 + labs(color="Latitude")
p1

p2 = ggplot(merged.leaf, aes(merged.leaf$NITROGEN_C, merged.leaf$SPP_O1_BAS, color=merged.leaf$GENUS)) + geom_point()
p2 = p2 + geom_smooth(method = lm, se = FALSE)
p2 = p2 + facet_wrap(~FIXER) + labs(x="Soil N", y="%BA", title="%BA vs. Soil N by Fixer Status")
p2 = p2 + labs(color="Genus")
p2

windows(title="%BA vs. Soil N by Ecoregion")
p3 = ggplot(newdata, aes(newdata$SOILN, newdata$lin.spp, color=newdata$ECOREGION)) + geom_point()
p3 = p3 + geom_smooth(method = lm, se = FALSE)
p3 = p3 + facet_wrap(~GENUS) + labs(x=expression(paste("Soil N (gN/m"^"2)")), y="%BA in Dominant Species", title="%BA vs. Soil N by Ecoregion")
p3 = p3 + labs(color="Ecoregion")
p3
#
GENUS1 <- as.factor(newdata$SPP_O1_ABB)
p3 = ggplot(newdata, aes(newdata$SOILN, newdata$lin.spp, color=newdata$ECOREGION)) + geom_point()
p3 = p3 + geom_smooth(method = lm, se = FALSE)
p3 = p3 + facet_wrap(~as.factor(newdata$SPP_O1_ABB)) + labs(x=expression(paste("Soil N (gN/m"^"2)")), y="%BA in Dominant Species", title="%BA vs. Soil N by Ecoregion")
p3 = p3 + labs(color="Ecoregion")
p3
#
windows(title="%BA vs. Latitude by Ecoregion")
p4 = ggplot(newdata, aes(newdata$LAT, newdata$lin.spp, color=newdata$ECOREGION)) + geom_point()
p4 = p4 + geom_smooth(method = lm, se = FALSE)
p4 = p4 + facet_wrap(~GENUS) + labs(x="Latitude", y="%BA in Dominant Species", title="%BA vs. Soil N by Ecoregion")
p4 = p4 + labs(color="Ecoregion")
p4

p <- ggplot(merged.leaf, aes(SPP_O1_BAS, NITROGEN_C))
p+geom_point(aes(color = factor(FIXER)))+labs(x="% BA in Dom Species", y="Soil N")+ggtitle("Soil N vs. % BA in Dominant Species by Fixer Status")+scale_y_continuous(limits=c(-1,900))+scale_x_continuous(limits=c(-1,100))

p <- ggplot(merged.leaf, aes(SPP_O1_BAS, NITROGEN_C))
p+geom_point(aes(color = (ASA_R)))+labs(x="% BA in Dom Species", y="Soil N (gN/m2)")+ggtitle("Soil N vs. % BA in Dominant Species by Stand Age")+scale_y_continuous(limits=c(-1,900))+scale_x_continuous(limits=c(-1,100))+ labs(color="Stand Age")

########################## 
##
##plotting locations on map
##
###########################

# Mapping Sites -----------------------------------------------------------

##using rworldmap
library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-125,-70), ylim = c(35, 45), asp = 1)
points(merged.soil.prod$lon, merged.soil.prod$lat, col = "red", cex = .6)

##simple ggmap
# visited <- c("SFO", "Chennai", "London", "Melbourne", "Johannesburg, SA")
# ll.visited <- geocode(visited)
# visit.x <- ll.visited$lon
# visit.y <- ll.visited$lat
mp <- NULL
mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
mp <- ggplot() +   mapWorld
#Now Layer the plots on top
mp1 <- mp+ geom_point(aes(x=merged.soil.prod$LONGITUDE.x, y=merged.soil.prod$LATITUDE.x) ,color="blue", size=3) 
mp1 <- mp1 + labs(color="Ecoregion")
mp1

##simple ggmap just usa
# adapted from: http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
usa <- map_data("usa") # we already did this, but we can do it again
mp <- ggplot() + geom_polygon(data = usa, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3)
mp1 <- mp+ geom_point(aes(x=merged.soil.prod$LONGITUDE.x, y=merged.soil.prod$LATITUDE, color=merged.soil.prod$ECOREGION) , size=2) + labs(x=" ", y=" ")
mp1 + labs(color='Ecoregion')
mp1
#deeppink4
## using ggmap
library(ggplot2)
library(ggmap)
# creating a sample data.frame with your lat/lon points
lon <- c(-38.31,-35.5)
lat <- c(40.96, 37.5)
df <- as.data.frame(cbind(lon,lat))
# getting the map
mapgilbert <- get_map(location = c(lon = mean(df$lon), lat = mean(df$lat)), zoom = 4,
                      maptype = "satellite", scale = 2)
# plotting the map with some points on it
ggmap(mapgilbert) +
  geom_point(data = df, aes(x = lon, y = lat, fill = "red", alpha = 0.8), size = 5, shape = 21) +
  guides(fill=FALSE, alpha=FALSE, size=FALSE)

###################################################################################
########################## 
##
##correlation plot
##
###########################

#correlation plot
require(GGally)
ggpairs(merged.soil.prod[, c("SPP_O1_BAS", "NITROGEN_C", "LATITUDE.x", "MAP.x")])
ggpairs(data=merged.soil.prod, # data.frame with variables
        title="productivity data", # title of the plot
        colour = "SPP_O1_ABB") # aesthetics, ggplot2 style

pm = ggpairs(data=merged.soil.prod,
             columns = c(10,14,16,18,52,58),
             upper = list(continuous = "density"),
             lower = list(combo = "facetdensity"),
             title="productivity data")
print(pm)

##############################
########################## 
##
## multiple regression plots
##
###########################

require(datasets)
require(GGally)
require(ggplot2)

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = merged.soil.prod, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

g = ggpairs(merged.soil.prod,columns = c(10,14,16,18,52,58), lower = list(continuous = my_fn))
g

###################################################################
########################## 
##
## glm models
##
###########################
library(lme4)
library(LMERConvenienceFunctions)

#example from class
surv.glm<-glm(sdata$survival~sdata$height.ADJ, family=binomial)
surv.glmer<-glmer(sdata$survival~sdata$height.ADJ+(1|sdata$plot),family=binomial)
# 
## old models --------------
# # Nitrogen only glm model
# glm1 <- glm(merged.soil.prod$SPP_O1_BAS~merged.soil.prod$NITROGEN_C)
# summary(glm1)   #AIC = 1295.8
# confint(glm1) # 95% CI for the coefficients
# exp(coef(glm1)) # exponentiated coefficients
# exp(confint(glm1)) # 95% CI for exponentiated coefficients
# predict(glm1, type="response") # predicted values
# residuals(glm1, type="deviance") # residuals
# windows(title="Nitrogen Only Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm1, ask=FALSE)
# lm1 <- lm(merged.soil.prod$SPP_O1_BAS~merged.soil.prod$NITROGEN_C)
# 
# 
# # Nitrogen + Climate + Plot + Soil
# glm2 <- glm(SPP_O1_BAS~NITROGEN_C+LATITUDE.x+ELEVATION.x+ASA+LAI_O+MAT+MAP.y+BULK_DENSI+CARBON_CON+PH_OF_SOIL, data=merged.soil.prod)
# summary(glm2)   #AIC = 662.88
# confint(glm2)
# exp(coef(glm2))
# exp(confint(glm2))
# windows(title="Nitrogn + Climate + Plot + Soil Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm2, ask=FALSE)
# 
# # Random combo of factors
# glm3 <- glm(SPP_O1_BAS~NITROGEN_C+LATITUDE.x+ASA+LAI_O, data=merged.soil.prod)
# summary(glm3)   # AIC = 1209.2
# confint(glm3)
# exp(coef(glm3))
# windows(title="Random Combo Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm3, ask=FALSE)
# 
# # Nitrogen + climate
# # climate = latitude, mat, map
# # Soil = nitrogen, bulk density, carbon content, pH
# # plot effects = latitude, elevation, ASA
# # LAI
# glm4 <- glm(SPP_O1_BAS~NITROGEN_C+LATITUDE.x+MAT+MAP.y, data=merged.soil.prod)
# summary(glm4) #AIC = 1301.4
# windows(title="Nitrogen + Climate Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm4, ask=FALSE)
# 
# # Nitrogn + Soil
# glm5 <- glm(SPP_O1_BAS~NITROGEN_C+BULK_DENSI+CARBON_CON+PH_OF_SOIL, data=merged.soil.prod)
# summary(glm5)   #AIC=769.08
# windows(title="Nitrogen + Soil Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm5, ask=FALSE)
#   
# # Nitrogen + Plot
# glm6 <- glm(SPP_O1_BAS~NITROGEN_C+LATITUDE.x+ELEVATION.x+ASA, data=merged.soil.prod)
# summary(glm6)   #AIC=1299.8
# windows(title="Nitrogen + Plot Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm6, ask=FALSE)
# 
# # Plot only
# glm7 <- glm(SPP_O1_BAS~LATITUDE.x+ELEVATION.x+ASA, data=merged.soil.prod)
# summary(glm7)   #AIC=1408.8
# windows(title="Plot ONly Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm7, ask=FALSE)
# 
# # Soil only (without nitrogen)
# glm8 <- glm(SPP_O1_BAS~BULK_DENSI+CARBON_CON+PH_OF_SOIL, data=merged.soil.prod)
# summary(glm8)   #AIC=767.12
# windows(title="Soil ONly Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm8, ask=FALSE)
# 
# # Climate only
# glm9 <- glm(SPP_O1_BAS~LATITUDE.x+MAT+MAP.y, data=merged.soil.prod)
# summary(glm9)   #AIC=1845.2
# windows(title="Climate ONly Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm9, ask=FALSE)
# 
# # survival -> SPP_O1_BAS, hieight.ADJ -> NITROGEN_C
# plot(merged.soil.prod$SPP_O1_BAS~merged.soil.prod$NITROGEN_C,ylab="SPP", xlab="nitrogen")
# lines(fitted(glm2)[order(merged.soil.prod$NITROGEN_C)]~merged.soil.prod$NITROGEN_C[order(merged.soil.prod$NITROGEN_C)], col='red')#need to "order" from lowest to largest height in order to plot this as a line.
# legend("right", col=c("black", "red"), c("observed", "predicted"), pch=1)
# 
# plot(merged.soil.prod$SPP_O1_BAS~merged.soil.prod$NITROGEN_C,ylab="survival", xlab="height")
# #plot line for each plot
# for (i in 1:max(sdata$plot)){
#   lines(fitted(surv.glmer)[sdata$plot==i][order(sdata$height.ADJ[sdata$plot==i])]~sdata$height.ADJ[sdata$plot==i][order(sdata$height.ADJ[sdata$plot==i])], col=i)
# }
# 
# ################ playing around with soil N as variable of interest
# glm11 <- glm(NITROGEN_C~LATITUDE.x+MAT+MAP.y+ELEVATION.x+PH_OF_SOIL+CARBON_CON+NumDom+BULK_DENSI+AG_BIOMASS+LAI_O+ASA, data=merged.soil.prod)
# summary(glm11)
# par(mfrow=c(2,2))
# plot(glm11, ask=FALSE)
# 
# glm11 <- glm(NITROGEN_C~LATITUDE.x+MAT+MAP.y+ELEVATION.x+PH_OF_SOIL+CARBON_CON+NumDom+BULK_DENSI+AG_BIOMASS+LAI_O+ASA, data=merged.soil.prod)
# summary(glm11)
# par(mfrow=c(2,2))
# plot(glm11, ask=FALSE)

##################################
windows(title="bodysize vs. survival") # creates a window with title
plot(merged.soil.prod$SPP_O1_BAS,merged.soil.prod$NITROGEN_C, xlab="BAS",ylab="Nitrogen") # plot with body size on x-axis and survival (0 or 1) on y-axis
g=glm(SPP_O1_BAS~NITROGEN_C,data=merged.soil.prod) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
# curve(predict(g,data.frame(NITROGEN_C=x),type="response"),add=TRUE,col='red') # draws a curve based on prediction from logistic regression model
# points(merged.soil.prod$NITROGEN_C,fitted(g),pch=20) 

##################################################################
########################## 
##
## working with N-fixer data
##
###########################

# What environmental characteristics favor N-fixer growth?
# par(mfrow=c(1,1))
# plot(merged.leaf$NITROGEN_C, merged.leaf$AG_BIOMASS)
# lm12 <- lm(AG_BIOMASS~NITROGEN_C, data=merged.leaf)
# r2 <- summary(lm12)$R.squared
# r2.adj <- summary(lm12)$adj.r.squared
# int1 <- signif(coef(lm12)[1], digits = 2)
# slope1 <- signif(coef(lm12)[2], digits = 2)
# textlab <- paste("y = ",slope1,"x + ",int1, sep="")
# p1<- ggplot(merged.leaf, aes(x=NITROGEN_C, y=SPP_O1_BAS)) 
# p1 + geom_point(size=1) + geom_smooth(method=lm, se=TRUE) +  # Hollow circles, Add linear regression line 
#   annotate("text", x = 500, y = 80, label = textlab, color="black", size = 5, parse=FALSE) +
#   labs(x="Soil Nitrogen (gN/m2)", y="%BA in Dominant (%)")+ggtitle("Soil Nitrogen vs. %BA")
# # Nitrogen + Climate + Plot + Soil
# glm12 <- glm(SPP_O1_BAS~NITROGEN_C+LATITUDE.x+ELEVATION.x+ASA+LAI_O+MAT+MAP.y+BULK_DENSI+CARBON_CON+PH_OF_SOIL, data=merged.leaf)
# summary(glm12)   #AIC = 3774.3
# confint(glm12)
# exp(coef(glm12))
# exp(confint(glm12))
# windows(title="Nitrogn + Climate + Plot + Soil Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm12, ask=FALSE)
# # Plot with fixer status
# fixers <- merged.leaf$FIXER
# sum(fixers) #47 fixers total
# lm13 <- lm(SPP_O1_BAS~FIXER, data=merged.leaf)
# r2 <- summary(lm13)$R.squared
# r2.adj <- summary(lm13)$adj.r.squared
# int1 <- signif(coef(lm13)[1], digits = 2)
# slope1 <- signif(coef(lm13)[2], digits = 2)
# textlab <- paste("y = ",slope1,"x + ",int1, sep="")
# ggplot something
 p1<- ggplot(merged.leaf, aes(x=FIXER, y=AG_BIOMASS)) 
 p1 + geom_point(size=1) + geom_smooth(method=lm, se=TRUE) +  # Hollow circles, Add linear regression line 
   annotate("text", x = 0.5, y = 8000, label = textlab, color="black", size = 5, parse=FALSE) +
   labs(x="Fixer Status", y="Biomass (gC/m2)")+ggtitle("Fixer Status vs. Biomass")
# More old Models -------------
 # glm13 <- glm(SPP_O1_BAS~FIXER+NITROGEN_C+LATITUDE.x+ELEVATION.x+ASA+LAI_O+MAT+MAP.y+BULK_DENSI+CARBON_CON+PH_OF_SOIL, data=merged.leaf)
# summary(glm13)   #AIC = 3775.9
# confint(glm13)
# exp(coef(glm13))
# exp(confint(glm13))
# windows(title="Nitrogn +Fixer Status + Climate + Plot + Soil Model") # creates a window with title
# par(mfrow=c(2,2))
# plot(glm13, ask=FALSE)
# #with fixer status and age
# glm14 <- glm(SPP_O1_BAS~FIXER+NITROGEN_C+ASA_R, data=merged.leaf)
# summary(glm14)   #AIC=18805
# glm15 <- glm(SPP_O1_BAS~FIXER+NITROGEN_C+ASA_R, data=merged.leaf)
# summary(glm15)   #7679.8
# glm16 <- glm(SPP_O1_BAS~FIXER+NITROGEN_C+ASA_R+MAT+MAP.y+CARBON_CON, data=merged.leaf)
# summary(glm16)
# anova(glm14,glm16)

##############################
library(lme4)
y <- merged.leaf$SPP_O1_BAS
y1 <- newdata$lin.spp
x <- (merged.leaf$NITROGEN_C-mean(merged.leaf$NITROGEN_C))/sd(merged.leaf$NITROGEN_C)
a <- (merged.leaf$ASA-mean(merged.leaf$ASA))/sd(merged.leaf$ASA)
grouping <- as.factor(merged.leaf$FIXER)
# #not working
# my.lmer <-  lmer(y ~ x + a + x * a + (1+ x | grouping), data = merged.leaf)
# summary(my.lmer)
# #working
# fm2 <- lmer(y ~ 1 + (1|FIXER), merged.leaf)
# summary(fm2)
# AIC(fm2)   #8808
# fm3 <- lmer(y ~ 1 + (1|FIXER) + (1|ASA) + (1|NITROGEN_C), merged.leaf)
# summary(fm3)
# AIC(fm3)   #4152
# fm4 <- lmer(y1 ~ 1 + (1|FIXER) + (1|ASA), merged.leaf)
# summary(fm4)
# AIC(fm4)
# # Make climate variable
# climate <- merged.leaf$LATITUDE*merged.leaf$MAT_C*merged.leaf$MAP.x*merged.leaf$ELEVATION.x
# c.climate <- (climate-mean(climate))/sd(climate)
# soil <- merged.leaf$NITROGEN_C*merged.leaf$CARBON_CON*merged.leaf$PH_OF_SOIL*merged.leaf$BULK_DENSI
# c.soil <- (soil-mean(soil))/sd(soil)
# c.ASA <- (merged.leaf$ASA-mean(merged.leaf$ASA))/sd(merged.leaf$ASA)
# fm4 <- lmer(y ~ 1 + (1|c.ASA) + (1|c.soil))
# 
# # each has its own slope
# lmer(y ~ x + (time | subjects), data=data)
# lmer1 <- lmer(y ~ x + (x|ASA), merged.leaf)

merged.leaf$FIXER = factor(merged.leaf$FIXER, 
                           levels=c(0,1), 
                           labels=c("Nonfixer","Fixer"))
table(merged.leaf$FIXER)

######################################################################

# Make Newdata ------------------------------------------------------------
# simplify dataframe further and do actual glmm models
myvars <- c("PLOT_ID", "ID.1", "LATITUDE", "LONGITUDE","LAI_O","HEIGHTC_m","FIXER","GENUS","SPECIES","LEAF_CN","ECOREGION","ELEVATION.x","MAT_C","MAP.x","ASA","ASA_R","ASA_R1","SPP_O1_ABB","SPP_O1_BAS","BULK_DENSI","CARBON_CON","NITROGEN_C","PH_OF_SOIL","NumDom")
newdata <- merged.leaf[myvars]
colnames(newdata) <- c("ID","ID.1","LAT","LONG","LAI","HEIGHT","FIXER","GENUS","SPECIES","LEAF_CN","ECOREGION","ELEV","MAT","MAP","ASA","ASA_R","ASA_R1","SPP_O1_ABB","SPP_O1_BAS","BULKDEN","SOILC","SOILN","SOILPH","NUMDOM")
head(newdata)
table(newdata$GENUS)

# Making newdata3 from the merged.soil.prod data
myvars <- c("PLOT_ID", "LATITUDE.x", "LONGITUDE.y","LAI_O","HEIGHTC","ECOREGION.y","ELEVATION.x","MAT_C","MAP.x","ASA","ASA_R","ASA_R1","SPP_O1_ABB","SPP_O1_BAS","BULK_DENSI","CARBON_CON","NITROGEN_C","PH_OF_SOIL","NumDom")
newdata3 <- merged.soil.prod[,myvars]
colnames(newdata3) <- c("ID","LAT","LONG","LAI","HEIGHT","ECOREGION","ELEV","MAT","MAP","ASA","ASA_R","ASA_R1","SPP_O1_ABB","SPP_O1_BAS","BULKDEN","SOILC","SOILN","SOILPH","NUMDOM")
head(newdata3)
table(newdata3$SPP_O1_ABB)
newdata3 <- newdata3[complete.cases(newdata3),]

##outliers
mahal = mahalanobis(newdata[ ,c(4,15,18,21)],
                    colMeans(newdata[ ,c(4,15,18,21)],na.rm=TRUE),
                    cov(newdata[ ,c(4,15,18,21)],use="pairwise.complete.obs"))
summary(mahal)
cutoff = qchisq(1-0.001,ncol(newdata[ ,c(4,15,18,21)]))
cutoff
summary(mahal < cutoff)   #16 outliers, 209 NAs, 790 normal
newdata1 = newdata[mahal < cutoff , ]   #Exclude outliers

##multicollinearity
corr <- cor(newdata[, -c(1,2,5,8,9,11,16,17,18,19,24,25)],use="complete")
library(corrplot)
corrplot(corr, method="color", order="AOE")
symnum(corr)
corrplot(corr, method="color", order="AOE",type="lower")

##assumptions
fake = lm(SPP_O1_BAS ~ ., data=newdata)
fitted = scale(fake$fitted.values)
standardized = rstudent(fake)

##linearity
qqnorm(standardized)
abline(0,1)

##normality
hist(standardized)

##intercept only model
library(nlme)
model2 = lme(SPP_O1_BAS ~ 1,
             data = newdata,
             method = "ML",
             na.action = "na.omit",
             random = ~1|ID)
summary(model2)
model4 = lme(SPP_O1_BAS ~ c.SOILN + c.ASA,
             data = newdata,
             method = "ML",
             na.action = "na.omit",
             random = ~1|FIXER)
summary(model4)

# scale predictors----------------------------
##
c.SOILN = scale(newdata$SOILN)
c.ASA = scale(newdata$ASA)
c.SOILC = scale(newdata$SOILC)
c.SOILPH = scale(newdata$SOILPH)
c.BULKDEN = scale(newdata$BULKDEN)
c.ASA = scale(newdata$ASA)
c.LAT = scale(newdata$LAT)
c.LONG = scale(newdata$LONG)
c.MAT = scale(newdata$MAT)
c.MAP = scale(newdata$MAP)
##
c1.SOILN = scale(newdata3$SOILN)
c1.ASA = scale(newdata3$ASA)
c1.SOILC = scale(newdata3$SOILC)
c1.SOILPH = scale(newdata3$SOILPH)
c1.BULKDEN = scale(newdata3$BULKDEN)
c1.ASA = scale(newdata3$ASA)
c1.LAT = scale(newdata3$LAT)
c1.LONG = scale(newdata3$LONG)
c1.MAT = scale(newdata3$MAT)
c1.MAP = scale(newdata3$MAP)

#means for reverse transform (later)
m.SOILN = mean(newdata$SOILN, na.rm=TRUE)
m.SOILPH = mean(newdata$SOILPH, na.rm=TRUE)
m.LAT = mean(newdata$LAT, na.rm=TRUE)
m.MAT = mean(newdata$MAT, na.rm=TRUE)

#sds for reverse transform (later)
s.SOILN = sd(newdata$SOILN, na.rm=TRUE)
s.SOILPH = sd(newdata$SOILPH, na.rm=TRUE)
s.LAT = sd(newdata$LAT, na.rm=TRUE)
s.MAT = sd(newdata$MAT, na.rm=TRUE)

model3 = lmer(SPP_O1_BAS ~ c.SOILN + c.ASA + (1|FIXER), data = newdata)
summary(model3)

tapply(newdata$SOILN, newdata$FIXER, mean, na.rm=T)
#nonfixers have higher soil N (205.3 compared to 146.2)
model4 = lmer(SPP_O1_BAS ~ c.SOILN + c.ASA + (c.SOILC|FIXER), data = newdata)
summary(model4)
anova(model3,model4)

####################################################
#############
## working with newdata
##
#############

# Glm Models --------------------------------------------------------------

lin.spp <- log(newdata$SPP_O1_BAS+1)
lin.spp.mean <- mean(lin.spp, na.rm=TRUE)
lin.spp.sd <- sd(lin.spp, na.rm=TRUE)
hist(lin.spp, prob=TRUE)
curve(dgamma(x, shape=lin.spp.mean^2/lin.spp.sd^2, rate=lin.spp.mean/lin.spp.sd^2), add = TRUE, col = 'blue') 

newdata <- cbind(newdata,lin.spp)
colnames(newdata)

glm1 <- glm(lin.spp ~ c.LAT+c.LONG+c.MAT+c.MAP+c.ASA+c.SOILN+c.SOILC+c.SOILPH+c.BULKDEN, data=newdata, na.action=na.exclude)
summary(glm1) # significant covariates: lat, mat, asa, soilN, soilpH
AIC(glm1) #AIC = 2.283
logLik(glm1)
par(mfrow=c(1,1))
plot(predict(glm1),newdata$lin.spp,col='blue',cex=0.5,xlab="Predicted",ylab="%BA Data",main="Data vs. Predictions for glm1")
coefplot(glm1, vertical = FALSE, intercept = FALSE, var.las = 1, frame.plot = TRUE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
par(mfrow=c(2,2))
plot(glm1, ask=FALSE)
lm1 <- lm(lin.spp ~ c.LAT+c.LONG+c.MAT+c.MAP+c.ASA+c.SOILN+c.SOILC+c.SOILPH+c.BULKDEN, data=newdata, na.action=na.exclude)

#
glm2 <- glm(lin.spp ~ c.SOILN+c.SOILPH, data=newdata, na.action=na.exclude)
summary(glm2)
AIC(glm2) #AIC=137.6
logLik(glm2)
plot(predict(glm2),newdata$lin.spp,col='blue',cex=0.5,xlab="Predicted",ylab="%BA Data",main="Data vs. Predictions for glm2")
coefplot(glm2, vertical = FALSE, intercept = FALSE, var.las = 1, frame.plot = TRUE, newNames = c(c.SOILPH="Soil pH",c.SOILN="Soil N"))
par(mfrow=c(2,2))
plot(glm2, ask=FALSE)
#
glm3 <- glm(lin.spp ~ c.LAT+c.MAT+c.ASA+c.SOILN+c.SOILPH, data=newdata, na.action=na.exclude)
summary(glm3) #all covariates significant, ASA less significant than others
AIC(glm3) #AIC = -2.224
logLik(glm3)
dev.off()
plot(predict(glm3),newdata$lin.spp,col='blue',cex=0.5,xlab="Predicted",ylab="%BA Data",main="Data vs. Predictions for glm3")
coefplot(glm3, vertical = FALSE, intercept = FALSE, var.las = 1, frame.plot = TRUE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LAT="Latitude"))
par(mfrow=c(2,2))
plot(glm3, ask=FALSE)
#
glm4 <- glm(lin.spp ~ c.LAT+c.SOILN+c.ASA+c.SOILPH, data=newdata, na.action=na.exclude)
summary(glm4)
AIC(glm4) #AIC=33.36
logLik(glm4)
display(glm4)
dev.off()
plot(predict(glm4),newdata$lin.spp,col='blue',cex=0.5,xlab="Predicted",ylab="%BA Data",main="Data vs. Predictions for glm4")
coefplot(glm4, vertical = FALSE, intercept = FALSE, var.las = 1, frame.plot = TRUE, newNames = c(c.SOILPH="Soil pH",c.SOILN="Soil N",c.ASA="Age",c.LAT="Latitude"))
par(mfrow=c(2,2))
plot(glm4, ask=FALSE)
#
glm5 <- glm(lin.spp ~ c.LAT+c.SOILN+c.MAT+c.SOILPH, data=newdata, na.action=na.omit)
summary(glm5)
AIC(glm5) #AIC=6.36
logLik(glm5)
display(glm5)
dev.off()
plot(predict(glm5),newdata$lin.spp,col='blue',cex=0.5,xlab="Predicted",ylab="%BA Data",main="Data vs. Predictions for glm4")
coefplot(glm5, vertical = FALSE, intercept = FALSE, var.las = 1, frame.plot = TRUE, newNames = c(c.SOILPH="Soil pH",c.SOILN="Soil N",c.ASA="Age",c.LAT="Latitude"))
par(mfrow=c(2,2))
plot(glm5, ask=FALSE)
# it doesn't matter for the glm's which data set you use (newdata or newdata3)
glm5b <- glm(lin.spp ~ c.LAT+c.SOILN+c.MAT+c.SOILPH, data=newdata3, na.action=na.omit)
summary(glm5b)
#
anova.glm <- anova(glm1,glm2,glm3,glm4) #compare glm models
capture.output(anova.glm,file="test.doc")
multiplot(glm1, glm2, glm3, glm4, glm5, single=FALSE, intercept=FALSE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))

# Lmer Models -------------------------------------------------------------

f.FIXER <- as.factor(newdata$FIXER)
c.ID <- as.character(newdata$ID.1)
f.ID <- as.factor(c.ID)
f.GENUS <- as.factor(newdata$GENUS)
f.ECOREGION <- as.factor(newdata$ECOREGION)

lmer1 <- lmer(lin.spp ~ c.LAT+c.LONG+c.MAT+c.MAP+c.ASA+c.SOILN+c.SOILC+c.SOILPH+c.BULKDEN+(1|f.FIXER), data=newdata, na.action = na.exclude)
#add REML=0 to switch method to maximum likelihood
summary(lmer1) # significant covariates: lat, mat, asa, soilN, soilpH, random = fixer
AIC(lmer1) #AIC = 65.6
logLik(lmer1)
display(lmer1)   #DIC = -81, deviance = -19.7
confint(lmer1)
dev.off()
coefplot(lmer1, vertical = FALSE, var.las = 1, frame.plot = TRUE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
plot(predict(lmer1),newdata$lin.spp,col='blue',cex=0.5,xlab="Predicted",ylab="%BA Data",main="Data vs. Predictions for lmer1")
lm1 <- lm(lin.spp~predict(lmer1), data=newdata)
r2.adj <- summary(lm1)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int1 <- signif(coef(lm1)[1], digits = 2)
slope1 <- signif(coef(lm1)[2], digits = 2)
eqn <- bquote(italic(y) == .(int1) + .(slope1)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(4.4,4,eqn,pos=1)
abline(a=int1,b=slope1,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
par(mfrow=c(3,3))
plotLMER.fnc(lmer1, ylab="%BA")
fitted.values(lmer1)
coef(lmer1) #intercept for each level in SPECIES
ranef(lmer1)
randoms <- REsim(lmer1, n.sims = 500) # looks pretty good
plotREsim(randoms)
dev.off()
hist(randoms$mean, breaks=10)
sjp.lmer(lmer1, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)

# fixer is not significant for intercept
lmer2 <- lmer(lin.spp ~ c.LAT+c.MAT+c.MAP+c.ASA+c.SOILN+c.SOILC+c.SOILPH+c.BULKDEN+(1|f.ECOREGION), data=newdata, na.action = na.omit)
summary(lmer2) # long removed because it is highly correlated, random = ecoregion
AIC(lmer2) #AIC = 52.67
logLik(lmer2)
display(lmer2)  #DIC = -66.9, deviance = -18.1
confint(lmer2)
plot(predict(lmer2),newdata$lin.spp,col='blue',cex=0.5,xlab="Predicted",ylab="%BA Data",main="Data vs. Predictions for Lmer2")
lm <- lm(lin.spp~predict(lmer2), data=newdata)
r2.adj <- summary(lm)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int <- signif(coef(lm)[1], digits = 2)
slope <- signif(coef(lm)[2], digits = 2)
eqn <- bquote(italic(y) == .(int) + .(slope)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(4.4,4,eqn,pos=1)
abline(a=int,b=slope,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
par(mfrow=c(3,3))
plotLMER.fnc(lmer2, ylab="%BA")
coefplot(lmer2, vertical = FALSE, var.las = 1, frame.plot = TRUE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
fitted.values(lmer2)
coef(lmer2) #intercept for each level in SPECIES
ranef(lmer2)
dev.off()
randoms2 <- REsim(lmer2, n.sims = 500) # looks pretty good
plotREsim(randoms2, labs=TRUE)
hist(randoms2$mean, breaks=10)
anova(lmer1,lmer2)
sjp.lmer(lmer2, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)
#
lmer3 <- lmer(lin.spp ~ c.LAT+c.MAT+c.MAP+c.ASA+c.SOILN+c.SOILC+c.SOILPH+c.BULKDEN+(1|f.GENUS), data=newdata, na.action = na.exclude)
summary(lmer3) # long removed because it is highly correlated, random = ecoregion
AIC(lmer3) #AIC = 52.67
logLik(lmer3)
display(lmer3)  #DIC = -88, deviance = -33.9
confint(lmer3)
par(mfrow=c(3,3))
plotLMER.fnc(lmer3, ylab="%BA")
par(mfrow=c(1,1))
coefplot(lmer3, vertical = FALSE, var.las = 1, frame.plot = TRUE, intercept = FALSE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
dev.off()
plot(predict(lmer3),newdata$lin.spp,col='blue',cex=0.5,xlab="Predicted",ylab="%BA Data",main="Data vs. Predictions for lmer3")
lm1 <- lm(lin.spp~predict(lmer3), data=newdata)
r2.adj <- summary(lm1)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int1 <- signif(coef(lm1)[1], digits = 2)
slope1 <- signif(coef(lm1)[2], digits = 2)
eqn <- bquote(italic(y) == .(int1) + .(slope1)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(4.5,4,eqn,pos=2)
abline(a=int1,b=slope1,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
fitted.values(lmer3)
coef(lmer3) #intercept for each level in SPECIES
ranef(lmer3)
randoms3 <- REsim(lmer3, n.sims = 500) # looks pretty good
plotREsim(randoms3, labs=TRUE)
hist(randoms3$mean, breaks=10)
anova(lmer1,lmer2)
sjp.lmer(lmer3, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)
#
lmer4 <- lmer(lin.spp ~ c.LAT+c.MAT+c.MAP+c.ASA+c.SOILN+c.SOILC+c.SOILPH+c.BULKDEN+(1|f.ID), data=newdata, na.action = na.exclude)
summary(lmer4) # long removed because it is highly correlated, random = plot ID new
AIC(lmer4) #AIC = 52.67
logLik(lmer4)
display(lmer4)
confint(lmer4)
#
lmer5 <- lmer(lin.spp ~ c.LAT+c.MAT+c.ASA+c.SOILN+c.SOILPH+(1|f.ECOREGION), data=newdata,na.action = na.omit)
summary(lmer5)
AIC(lmer5)  #AIC = 39.74
display(lmer5)   #DIC=-56.2, deviance=-16.2
confint(lmer5)
logLik(lmer5)
anova(lmer5)
plot(predict(lmer5),newdata$lin.spp,col='blue',cex=0.5, xlab="Predicted", ylab="Data",main="Predicted vs. Data for Lmer5")
lm1 <- lm(lin.spp~predict(lmer5), data=newdata)
r2.adj <- summary(lm1)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int1 <- signif(coef(lm1)[1], digits = 2)
slope1 <- signif(coef(lm1)[2], digits = 2)
eqn <- bquote(italic(y) == .(int1) + .(slope1)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(4.5,4,eqn,pos=2)
abline(a=int1,b=slope1,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
randoms5 <- REsim(lmer5, n.sims = 500) # looks pretty good
plotREsim(randoms5, labs=TRUE)
#
lmer6 <- lmer(lin.spp ~ c.LAT+c.MAT+c.SOILN+c.SOILPH+(1|c.ASA), data=newdata,na.action = na.exclude)
display(lmer6)
confint(lmer6)
#
lmer7 <- lmer(lin.spp ~ c.LAT+c.MAT+c.ASA+c.SOILN+c.SOILPH+(1|f.ECOREGION)+(1|f.GENUS), data=newdata,na.action = na.exclude)
summary(lmer7)
display(lmer7)
logLik(lmer7)
plot(predict(lmer7),newdata$lin.spp,col='blue',cex=0.5, xlab="Predicted", ylab="Data",main="Predicted vs. Data for Lmer7")
lm1 <- lm(lin.spp~predict(lmer7), data=newdata)
r2.adj <- summary(lm1)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int1 <- signif(coef(lm1)[1], digits = 2)
slope1 <- signif(coef(lm1)[2], digits = 2)
eqn <- bquote(italic(y) == .(int1) + .(slope1)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(4.5,4,eqn,pos=2)
abline(a=int1,b=slope1,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
#
lmer8 <- lmer(lin.spp ~ c.LAT+c.MAT+c.ASA+c.SOILN+c.SOILPH + (0+c.SOILN|f.GENUS), data=newdata, na.action = na.exclude)
display(lmer8)
logLik(lmer8)
#
lmer9 <- lmer(lin.spp ~ c.LAT + (c.LAT|f.ECOREGION), data=newdata, na.action = na.exclude)
display(lmer9)
confint(lmer9)

lmer10 <- lmer(lin.spp ~ c.SOILN+c.SOILPH + (c.SOILN|f.GENUS) + (c.SOILN|f.ECOREGION), data=newdata, na.action = na.exclude)
display(lmer10)
confint(lmer10)

lmer11 <- lmer(lin.spp ~ c.SOILN+c.SOILPH + c.MAT + c.SOILPH + (c.SOILN|f.GENUS), data=newdata, na.action = na.exclude)
display(lmer11)

lmer12 <- lmer(lin.spp ~ c.LAT+c.MAT+c.SOILN+c.SOILPH + (0+c.SOILN|f.GENUS), 
               data=newdata, na.action = na.omit)
display(lmer12)
logLik(lmer12)
sjp.lmer(lmer12, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)
par(mfrow=c(2,2))
plotLMER.fnc(lmer13, ylab="%BA")
par(mfrow=c(1,1))
coefplot(lmer12, vertical = FALSE, var.las = 1, frame.plot = TRUE, 
         intercept = FALSE, newNames = c(c.BULKDEN="Bulk Density",
                                         c.SOILPH="Soil pH",c.SOILC="Soil C",
                                         c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",
                                         c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
dev.off()
sjp.lmer(lmer12, y.offset = .4)
sjp.lmer(lmer12, type = "ri.slope")
sjp.lmer(lmer12, type="re.qq")
sjp.lmer(lmer12, type="fe")
plot(predict(lmer12),newdata$lin.spp,col='blue',cex=0.5, xlab="Predicted", ylab="Data",main="Predicted vs. Data for Lmer12")
lm1 <- lm(lin.spp~predict(lmer12), data=newdata)
r2.adj <- summary(lm1)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int1 <- signif(coef(lm1)[1], digits = 2)
slope1 <- signif(coef(lm1)[2], digits = 4)
eqn <- bquote(italic(y) == .(int1) + .(slope1)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(4.5,4,eqn,pos=2)
abline(a=int1,b=slope1,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
legend("bottomright", legend = c("1:1 line", "Regression Line"), col = c("grey", "red"), lty = c("dashed","solid"))
randoms12 <- REsim(lmer12, n.sims = 500) # looks pretty good
plotREsim(randoms12, labs=TRUE)
randoms12 <- REsim(lmer12, n.sims = 500) # looks pretty good
plotREsim(randoms12, labs=TRUE)
pred12 <- predict(lmer12)
plot(exp(pred12)~newdata$SOILN,col='blue',cex=0.5, xlab="Soil N (gN/m2)", ylab="%BA from model",main="Modeled %BA vs. Soil N for Lmer 12")
#
f1.GENUS <- as.factor(newdata3$SPP_O1_ABB)
lmer13 <- lmer(lin.spp ~ c.LAT+c.MAT+c.SOILPH+c.SOILN+(1|f.GENUS), 
               data=newdata,na.action = na.omit)
display(lmer13)
logLik(lmer13)
par(mfrow=c(2,2))
plotLMER.fnc(lmer13, ylab="%BA")
par(mfrow=c(1,1))
coefplot(lmer13, vertical = FALSE, var.las = 1, frame.plot = TRUE, 
         intercept = FALSE, newNames = c(c.BULKDEN="Bulk Density",
                                         c.SOILPH="Soil pH",c.SOILC="Soil C",
                                         c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",
                                         c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
dev.off()
sjp.lmer(lmer13, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)
sjp.lmer(lmer13, type="re.qq")
sjp.lmer(lmer13, type="ri.slope", geom.size=1)
plot(predict(lmer13, interval="confidence"),newdata$lin.spp,col='blue',cex=0.5, xlab="Predicted", ylab="Data",main="Predicted vs. Data for Lmer13")
lm1 <- lm(lin.spp~predict(lmer13), data=newdata)
r2.adj <- summary(lm1)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int1 <- signif(coef(lm1)[1], digits = 2)
slope1 <- signif(coef(lm1)[2], digits = 4)
eqn <- bquote(italic(y) == .(int1) + .(slope1)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(4.5,4,eqn,pos=2)
abline(a=int1,b=slope1,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
legend("bottomright", legend = c("1:1 line", "Regression Line"), col = c("grey", "red"), lty = c("dashed","solid"))
randoms13 <- REsim(lmer13, n.sims = 500) # looks pretty good
plotREsim(randoms13, labs=TRUE)
pred13 <- predict(lmer13)
plot(exp(pred13)~newdata$SOILN,col='blue',cex=0.5, xlab="Soil N (gN/m2)", ylab="%BA from model",main="Modeled %BA vs. Soil N for Lmer 13")
#
f1.GENUS <- as.factor(newdata3$SPP_O1_ABB)
lmer13b <- lmer(log(SPP_O1_BAS+1) ~ c1.LAT+c1.MAT+c1.SOILPH+c1.SOILN+(1|f1.GENUS), 
               data=newdata3,na.action = na.omit)
summary(lmer13b)
display(lmer13b)
dev.off()
plot(log(newdata3$SPP_O1_BAS),predict(lmer13b))
#
lmer14 <- lmer(lin.spp ~ c.LAT+c.MAT+c.SOILN+c.SOILPH+(1+c.SOILN|f.GENUS), 
               data=newdata,na.action = na.omit)
display(lmer14)
logLik(lmer14)
sjp.lmer(lmer14, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)
par(mfrow=c(2,2))
plotLMER.fnc(lmer14, ylab="%BA")
par(mfrow=c(1,1))
coefplot(lmer14, vertical = FALSE, var.las = 1, frame.plot = TRUE, 
         intercept = FALSE, newNames = c(c.BULKDEN="Bulk Density",
                                         c.SOILPH="Soil pH",c.SOILC="Soil C",
                                         c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",
                                         c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
dev.off()
sjp.lmer(lmer14, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)
sjp.lmer(lmer14, type="re.qq")
sjp.lmer(lmer14, type = "rs.ri", show.legend = TRUE)
p14 <- predict(lmer14)
plot(p14,newdata$lin.spp,col='blue',cex=0.5, xlab="Predicted", ylab="Data",main="Predicted vs. Data for Lmer14")
lm1 <- lm(lin.spp~predict(lmer14), data=newdata)
r2.adj <- summary(lm1)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int1 <- signif(coef(lm1)[1], digits = 2)
slope1 <- signif(coef(lm1)[2], digits = 4)
eqn <- bquote(italic(y) == .(int1) + .(slope1)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(4.5,4,eqn,pos=2)
abline(a=int1,b=slope1,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
legend("bottomright", legend = c("1:1 line", "Regression Line"), col = c("grey", "red"), lty = c("dashed","solid"))
randoms14 <- REsim(lmer14, n.sims = 500) # looks pretty good
plotREsim(randoms14, labs=TRUE)
pred14 <- predict(lmer14)
e.pred14 <- exp(pred14)
dev.off()
plot(exp(pred14)~newdata$SOILN,col='blue',cex=0.5, xlab="Soil N (gN/m2)", ylab="%BA from model",main="Modeled %BA vs. Soil N for Lmer 14")
summary(lm(log(e.pred14) ~ time, data = x, subset = Factor)) # I need the summary statistics to compare models
ggplot(x, aes(x = time, y = y, color = Factor)) + geom_point() + geom_smooth(method = "glm", family = gaussian(lin="log"), start=c(5,0))
#
lmer15 <- lmer(exp(lin.spp) ~ c.LAT+c.MAT+c.SOILN+c.SOILPH+(1|f.GENUS), data=newdata,na.action = na.exclude)
sjp.lmer(lmer15, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)
randoms15 <- REsim(lmer15, n.sims = 500) # looks pretty good
plotREsim(randoms15, labs=TRUE)
par(mfrow=c(2,2))
plotLMER.fnc(lmer15, ylab="%BA")
par(mfrow=c(1,1))
coefplot(lmer15, vertical = FALSE, var.las = 1, frame.plot = TRUE, intercept = FALSE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
dev.off()
sjp.lmer(lmer15, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)
sjp.lmer(lmer15, type="re.qq")
sjp.lmer(lmer15, type = "rs.ri", show.legend = TRUE)
plot(predict(lmer15),newdata$SPP_O1_BAS,col='blue',cex=0.5, xlab="Predicted", ylab="Data",main="Predicted vs. Data for Lmer14")
lm1 <- lm(exp(lin.spp)~predict(lmer15), data=newdata)
r2.adj <- summary(lm1)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int1 <- signif(coef(lm1)[1], digits = 2)
slope1 <- signif(coef(lm1)[2], digits = 4)
eqn <- bquote(italic(y) == .(int1) + .(slope1)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(90,45,eqn,pos=2)
abline(a=int1,b=slope1,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
randoms15 <- REsim(lmer15, n.sims = 500) # looks pretty good
plotREsim(randoms15, labs=TRUE)

### get plot ID working as a grouping factor
fm1<-lmer(lin.spp~c.LAT + c.MAT + c.ASA + c.SOILN + c.SOILPH  + (1 | f.ID))
## 2. check singularity
diag.vals <- getME(fm1,"theta")[getME(fm1,"lower") == 0]
any(diag.vals < 1e-6) # FALSE
## 3. recompute gradient and Hessian with Richardson extrapolation
devfun <- update(fm1, devFunOnly=TRUE)
if (isLMM(fm1)) {
  pars <- getME(fm1,"theta")
} else {
  ## GLMM: requires both random and fixed parameters
  pars <- getME(fm1, c("theta","fixef"))
}

if (require("numDeriv")) {
  cat("hess:\n"); print(hess <- hessian(devfun, unlist(pars)))
  cat("grad:\n"); print(grad <- grad(devfun, unlist(pars)))
  cat("scaled gradient:\n")
  print(scgrad <- solve(chol(hess), grad))
}
## compare with internal calculations:
fm1@optinfo$derivs
## 4. restart the fit from the original value (or
## a slightly perturbed value):
fm1.restart <- update(fm1, start=pars)
## 5. try all available optimizers
source(system.file("utils", "allFit.R", package="lme4"))
fm1.all <- allFit(fm1)
ss <- summary(fm1.all)
ss$ fixef ## extract fixed effects
ss$ llik ## log-likelihoods
ss$ sdcor ## SDs and correlations
ss$ theta ## Cholesky factors
ss$ which.OK ## which fits worked

#
lmer16 <- lmer(lin.spp ~ c.LAT+c.MAT+c.SOILPH+c.SOILN+(1|f.GENUS)+(1|f.ECOREGION), 
               data=newdata,na.action = na.omit)
summary(lmer16)
AIC(lmer16)
display(lmer16)
logLik(lmer16)
par(mfrow=c(2,2))
plotLMER.fnc(lmer16, ylab="%BA")
par(mfrow=c(1,1))
coefplot(lmer16, vertical = FALSE, var.las = 1, frame.plot = TRUE, 
         intercept = FALSE, newNames = c(c.BULKDEN="Bulk Density",
                                         c.SOILPH="Soil pH",c.SOILC="Soil C",
                                         c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",
                                         c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
dev.off()
sjp.lmer(lmer16, y.offset = .4,sort.est = "sort.all",facet.grid = FALSE)
sjp.lmer(lmer16, type="re.qq")
sjp.lmer(lmer16, type="ri.slope")
sjp.lmer(lmer16, type="rs.ri")
y <- predict(lmer16)
plot(y,newdata$lin.spp,col='blue',cex=0.5, xlab="Predicted", ylab="Data",main="Predicted vs. Data for Lmer16", na.rm=TRUE)
lm1 <- lm(lin.spp~predict(lmer16), data=newdata)
r2.adj <- summary(lm1)$adj.r.squared
r2.adj.round <- round(r2.adj, digits=3)
int1 <- signif(coef(lm1)[1], digits = 2)
slope1 <- signif(coef(lm1)[2], digits = 4)
eqn <- bquote(italic(y) == .(int1) + .(slope1)*italic(x) * "," ~~ r^2 == .(r2.adj.round))
text(4.5,4,eqn,pos=2)
abline(a=int1,b=slope1,col="red")
abline(a=0,b=1,col="grey",lty="dashed")
legend("bottomright", legend = c("1:1 line", "Regression Line"), col = c("grey", "red"), lty = c("dashed","solid"))
randoms13 <- REsim(lmer16, n.sims = 500) # looks pretty good
plotREsim(randoms16, labs=TRUE)
pred13 <- predict(lmer16)
plot(exp(pred16)~newdata$SOILN,col='blue',cex=0.5, xlab="Soil N (gN/m2)", ylab="%BA from model",main="Modeled %BA vs. Soil N for Lmer 13")
#
f1.ECOREGION <- as.factor(newdata3$ECOREGION)
lmer16b <- lmer(log(SPP_O1_BAS) ~ c1.LAT+c1.MAT+c1.SOILPH+c1.SOILN+(1|f1.GENUS)+(1|f1.ECOREGION), 
               data=newdata3,na.action = na.omit)
summary(lmer16b)
display(lmer16b)
#
f.ABB <- as.factor(newdata$SPP_O1_ABB)
lmer17 <- lmer(lin.spp ~ (c.LAT+c.MAT+c.SOILPH+c.SOILN)+(1|f.ABB), 
               data=newdata, na.action=na.omit)
summary(lmer17)
AIC(lmer17)
display(lmer17)
sjp.lmer(lmer17, type="ri.slope")
plot(predict(lmer17),newdata$lin.spp)
lmer18 <- lmer(lin.spp ~ (c.LAT+c.MAT+c.SOILPH+c.SOILN)+(1+c.SOILN|f.ABB),
               data=newdata, na.action=na.omit)
summary(lmer18)
display(lmer18)

#anova
anova(lmer13,lmer14,lmer16)

drop1(glm5)
model.aic.backward <- step(lmer14, direction = "backward", trace = 1)

# Plotting confidence intervals ----------------
pr01 <- profile(lmer13)
confint(pr01)
plot(pr01)
splom(pr01)

plot(ranef(lmer13))

#code for axis label: x=expression(paste("Soil N (gN/m"^"2)")) 

#show output of each variable effect nicely --------------------
library(effects)
plot(allEffects(lmer12),ylab="%BA")
plot(allEffects(lmer13),ylab="%BA")
plot(allEffects(lmer14),ylab="%BA")
plot(allEffects(lmer15), ylab="%BA")
plot(allEffects(lmer16), ylab="%BA")

## unscale coefficients function (Back transform)
get_real <- function(coef, scaled_covariate){
  # collect mean and standard deviation from scaled covariate
  mean_sd <- unlist(attributes(scaled_covariate)[-1])
  # reverse the z-transformation
  answer <- (coef * mean_sd[2]) + mean_sd[1]
  # this value will have a name, remove it
  names(answer) <- NULL
  # return unscaled coef
  return(answer)
}

coefs12 <- fixef(lmer12)[2]
lmer12_coefs <- get_real(coefs12, using_scale)

exp(posterior.mode(sb$Sol))
inv.logit(HPDinterval(sb$Sol, 0.95))

# Model Evaluation --------------------------------------------------------
## Examine Results
anova.lmer <- anova(lmer1,lmer2,lmer3,lmer4,lmer5,lmer6)
capture.output(anova.lmer,file="test.doc")

## Nice output table of model comparison
sjt.lmer(lmer12,lmer13,lmer14, show.header=TRUE,string.est="Estimate",string.ci="Conf. Int.",
         string.p="p-value", string.dv="Response",string.pred="Predictors",
         depvar.labels=c("Random Slope","Random Intercept","Random Slope & Intercept"),
         pred.labels=c("Latitude","Temperature","Soil N","Soil pH"),
         p.numeric=TRUE, group.pred=TRUE, show.r2=TRUE, show.aic=TRUE)

sjt.lmer(lmer13,lmer14,lmer16, show.header=TRUE,string.est="Estimate",string.ci="Conf. Int.",
         string.p="p-value", string.dv="Response",string.pred="Predictors",
         depvar.labels=c("Random Intercept (Genus)","Random Intercept & Slope (Genus)","Random Intercept (Genus & Ecoregion"),
         pred.labels=c("Latitude","Temperature","Soil N","Soil pH"),
         p.numeric=TRUE, group.pred=TRUE, show.r2=FALSE, show.aic=TRUE)

#plotting covariate plots
multiplot(lmer1, lmer2, lmer3, lmer4, lmer5, single=FALSE, intercept=FALSE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
multiplot(glm1, glm2, glm3, glm4, glm5, single=FALSE, intercept=FALSE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
multiplot(glm1, glm2, glm3, glm4, glm5, single=FALSE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))
multiplot(lmer12, lmer13, lmer14, single=FALSE, intercept=FALSE, newNames = c(c.BULKDEN="Bulk Density",c.SOILPH="Soil pH",c.SOILC="Soil C",c.SOILN="Soil N",c.ASA="Age",c.MAP="Precipitation",c.MAT="Temperature",c.LONG="Longitude",c.LAT="Latitude"))

#compare lmer and glm, not working
exactRLRT(m=lmer1, mA = NULL, m0=lm1)

plot(newdata$lin.spp ~ newdata$LAT, col = f.ECOREGION, las =1)
a=fixef(lmer9)
b=ranef(lmer9, condVar=TRUE)
#and plot each random slope
mean_=(b[[1]][2]+a[2])
abline(a = b[[1]][1,1]+a[1], b= mean_$c.LAT[1], col = "black")
abline(a = b[[1]][2,1]+a[1], b= mean_$c.LAT[2], col = "red")
abline(a = b[[1]][3,1]+a[1], b= mean_$c.LAT[3], col = "green")
abline(a = b[[1]][4,1]+a[1], b= mean_$c.LAT[4], col = "green")
abline(a = b[[1]][5,1]+a[1], b= mean_$c.LAT[5], col = "green")
abline(a = b[[1]][6,1]+a[1], b= mean_$c.LAT[6], col = "green")
abline(a=a[[1]], lty = 2)  #general response

#plot data by color plot
p1 <- predict(lmer12)
plot(p1 ~ newdata$SOILN, col=as.factor(newdata$ID.1), pch=16)

#get marginal and conciditional R2s
# require(piecewiseSEM)
sem.model.fits(lmer13)
sem.model.fits(lmer14)
sem.model.fits(lmer16)
sem.model.fits(lmer13b)
sem.model.fits(lmer16b)

#checking homogenous variance
dev.off()
plot(resid(lmer13) ~ fitted(lmer13),main="residual plot")
abline(h=0)
plot(resid(lmer14) ~ fitted(lmer14),main="residual plot")
abline(h=0)
plot(resid(lmer16) ~ fitted(lmer16),main="residual plot")
abline(h=0)
par(mfrow=c(2,2))
plot(glm5)

#look for outliers
library(influence.ME)
infl.lm0 <- influence(glm5)
dim(infl.lm0$infmat)
cook5 <- cooks.distance(glm5, na.action=na.omit)
plot(cook5)
cook5_outlier <- cook5[cook5 > 0.01]
# outliers = 35-41, 442-446 (12 in total)
cook5_outlier <- as.data.frame(cook5_outlier)
#35-41 = ID 10
newdata[35,]
#442-446 = ID 80
newdata[442,]
plot(newdata$lin.spp, newdata$SOILN)
points(4.615121,310,col='red')
points(3.637586,880,col='blue')
# blue definitely anamalously high nitrogen content
newdata2 <- newdata[-c(35,36,37,38,39,40,41), ]
c1.LAT <- scale(newdata2$LAT)
c1.MAT <- scale(newdata2$MAT)
c1.SOILN <- scale(newdata2$SOILN)
c1.SOILPH <- scale(newdata2$SOILPH)
f1.GENUS <- as.factor(newdata2$GENUS)
lmer13b <- lmer(lin.spp ~ c1.LAT+c1.MAT+c1.SOILPH+c1.SOILN+(1|f1.GENUS), 
               data=newdata2,na.action = na.omit)
display(lmer13b)

# Reverse Transform Coefficients -----------
#Lmer12
coef12<-fixef(lmer12)
int12 <- exp(coef12[1])   #63.45788
lat12 <- exp(coef12[2])*s.LAT+m.LAT   #46.67675
mat12 <- exp(coef12[3])*s.MAT+m.MAT   #10.67462
soiln12 <- exp(coef12[4])*s.SOILN+m.SOILN   #356.6241
soilph12 <- exp(coef12[5])*s.SOILPH+m.SOILPH   #6.724744
rand12 <- colMeans(ranef(lmer12)$f.GENUS)
slope12 <- rand12*s.SOILN+m.SOILN   #not working
plot(c.SOILN,newdata$lin.spp)
points(c.LAT,newdata$SPP_O1_BAS,col='red')
points(c.SOILPH,newdata$SPP_O1_BAS,col='blue')
points((c.LAT+c.SOILPH*c.SOILN),newdata$SPP_O1_BAS,col='green')
x3 <- newdata$SOILN[complete.cases(x,y)==TRUE]
x <- newdata$SOILN
y <- newdata$lin.spp
y3 <- newdata$lin.spp[complete.cases(x,y)==TRUE]
a_start<-50 #param a is the asymptote level
b_start<-300 #b is the (decay rate)a/b is initial slope
m<-nls(y~a*x/(b+x),start=list(a=a_start,b=b_start))
plot(x,y)
lines(x3,predict(m),lty=2,col="red",lwd=3)
re12 <- c(0.0130561267,-0.0024445471,-0.0191511210,0.0130086741,0.0023867456,0.0041741818,
          0.0009873565,0.0012794312,0.0012794312,0.0012794312, -0.0308142132,0.0155475901,
          0.0297414256,0.0025520824,-0.0031595472,0.0041741818,0.0004041615,-0.0343013912)
varCorr(lmer12)
se12 <- se.ranef(lmer12)
se12g <- se.ranef(lmer12)$f.GENUS
se12g <- se12g[1:18]   #all se as a list
se12g.exp <- exp(se12g)
se12.var <- var(se12g.exp)
sum(se12g^2)/18
var(ranef(lmer12)$f.GENUS)
xpred = predict(lmer12)
xpred1 = predict(lmer13)
xpred2 = predict(lmer14)
summary(lm(y3 ~ x$pred, subset = f.GENUS)) # I need the summary statistics to compare models
ggplot(x, aes(x = time, y = y, color = Factor)) + geom_point() + geom_smooth(method = "glm", family = gaussian(lin="log"), start=c(5,0))
plot(y3, predict(lmer12))
xyplot(exp(xpred) ~ exp(y3), type = c("p","r"), col.line = "red", xlab="Observed", ylab="Predicted",main="Predicted vs. Observed for lmer12")
xyplot(exp(xpred1) ~ exp(y3), type = c("p","r"), col.line = "red", xlab="Observed", ylab="Predicted",main="Predicted vs. Observed for lmer13")
xyplot(exp(xpred2) ~ exp(y3), type = c("p","r"), col.line = "red", xlab="Observed", ylab="Predicted",main="Predicted vs. Observed for lmer14")
xyplot(exp(xpred2) ~ x3, type = c("p","r"), col.line = "red", xlab="Soil N (gN/m2)", ylab="Modeled %BA",main="Modeled %BA vs. Soil N")

#Lmer13
coef13<-fixef(lmer13)
int13 <- exp(coef13[1])   #62.48565
lat13 <- exp(coef13[2])*s.LAT+m.LAT   #46.69519
mat13 <- exp(coef13[3])*s.MAT+m.MAT   #10.60877
soiln13 <- exp(coef13[4])*s.SOILN+m.SOILN   #357.6493
soilph13 <- exp(coef13[5])*s.SOILPH+m.SOILPH   #6.72721
rand13 <- colMeans(ranef(lmer13)$f.GENUS)
int13 <- rand13*s.SOILN+m.SOILN   #not working
#Lmer14
coef14<-fixef(lmer14)
int14 <- exp(coef14[1])   #62.485
lat14 <- exp(coef14[2])*s.LAT+m.LAT   #46.69593
mat14 <- exp(coef14[3])*s.MAT+m.MAT   #10.61285
soiln14 <- exp(coef14[4])*s.SOILN+m.SOILN   #358.1204
soilph14 <- exp(coef14[5])*s.SOILPH+m.SOILPH   #6.726795
rand14 <- colMeans(ranef(lmer14)$f.GENUS)
slope14 <- rand14*s.SOILN+m.SOILN   #202.1749
#Lmer16
coef16<-fixef(lmer16)
int16 <- exp(coef16[1])   #62.48565
lat16 <- exp(coef16[2])*s.LAT+m.LAT   #46.69519
mat16 <- exp(coef16[3])*s.MAT+m.MAT   #10.60877
soiln16 <- exp(coef16[4])*s.SOILN+m.SOILN   #396.35
soilph16 <- exp(coef16[5])*s.SOILPH+m.SOILPH   #6.5405
rand16 <- colMeans(ranef(lmer16)$f.GENUS)
rand16b <- colMeans(ranef(lmer16)$f.ECOREGION)

## MCMC plots -------------------------
newdata2 <- cbind(newdata, c.SOILN, c.LAT, c.SOILPH, c.SOILC, c.BULKDEN, c.LONG, c.ASA, c.MAT, c.MAP, f.GENUS, f.ECOREGION, f.FIXER, f.ID)
newdata.nona <- na.omit(newdata2)
mtest.bayes <- MCMCglmm(lin.spp ~ c.SOILN + c.LAT + c.SOILPH, random = ~f.GENUS, data = newdata.nona, 
                           verbose = FALSE, pl = TRUE, nitt = 50000, thin = 1, burnin = 1000)
m3.bayes <- MCMCglmm(lin.spp ~ c.MAP + c.SOILC+ c.BULKDEN +c.SOILN + c.LAT + c.SOILPH + c.ASA + c.MAT, random = ~f.ECOREGION, data = newdata.nona, 
                     verbose = FALSE, pl = TRUE, nitt = 50000, thin = 1, burnin = 1000)
summary(m3.bayes)
plot(m3.bayes)
m5.bayes <- MCMCglmm(lin.spp ~ c.SOILN + c.LAT + c.SOILPH + c.ASA + c.MAT, random = ~f.ECOREGION, data = newdata.nona, 
                     verbose = FALSE, pl = TRUE, nitt = 50000, thin = 1, burnin = 1000)
summary(m5.bayes)
plot(m5.bayes)
m7.bayes <- MCMCglmm(lin.spp ~ c.SOILN + c.LAT + c.SOILPH, random = ~f.GENUS + f.ECOREGION, data = newdata.nona, 
                     verbose = FALSE, pl = TRUE, nitt = 50000, thin = 1, burnin = 1000)
summary(m7.bayes)
windows(title="Posterior Distributions")
plot(m7.bayes,ask=TRUE)
#
m12.bayes <- MCMCglmm(lin.spp ~ c.LAT+c.MAT+c.SOILN+c.SOILPH, random=~idh(c.SOILN):f.GENUS, data = newdata.nona, 
                    verbose = FALSE, pl = TRUE, nitt = 50000, thin = 1, burnin = 1000)
summary(m12.bayes)
plot(m12.bayes)
m12.bayes$DIC
#
lmer13 <- lmer(lin.spp ~ c.LAT+c.MAT+c.SOILN+c.SOILPH+(1|f.GENUS), data=newdata,na.action = na.exclude)
m13.bayes <- MCMCglmm(lin.spp ~ c.LAT+c.MAT+c.SOILN+c.SOILPH, random = ~f.GENUS, data = newdata.nona, 
                      verbose = FALSE, pl = TRUE, nitt = 50000, thin = 1, burnin = 1000)
summary(m13.bayes)
plot(m13.bayes)
#
lmer14 <- lmer(lin.spp ~ c.LAT+c.MAT+c.SOILN+c.SOILPH+(1+c.SOILN|f.GENUS), data=newdata,na.action = na.exclude)
m14.bayes <- MCMCglmm(lin.spp ~ c.LAT+c.MAT+c.SOILN+c.SOILPH, random = ~f.GENUS+idh(c.SOILN):f.GENUS, data = newdata.nona, 
                      verbose = FALSE, pl = TRUE, nitt = 50000, thin = 1, burnin = 1000)
summary(m14.bayes)
plot(m14.bayes)
#
lmer16 <- lmer(lin.spp ~ c.LAT+c.MAT+c.SOILPH+c.SOILN+(1|f.GENUS)+(1|f.ECOREGION), 
               data=newdata,na.action = na.omit)
m16.bayes <- MCMCglmm(lin.spp ~ c.LAT+c.MAT+c.SOILPH+c.SOILN, random = ~f.GENUS+f.ECOREGION, data = newdata.nona, 
                      verbose = FALSE, pl = TRUE, nitt = 50000, thin = 1, burnin = 1000)
summary(m16.bayes)
plot(m16.bayes)

# Percent variance attributed to fixed effects ------------
Fixed12 <- fixef(lmer5)[2] * model.matrix(lmer5)[, 2] # Proportion of variance attributed to fixed effects.
VarF5 <- var(Fixed5) # Calculation of the variance in fitted values
FixedVar <- VarF5/(VarF5 + VarCorr(lmer5)$f.ECOREGION[1] + pi^2/3)
RandVar<-(VarCorr(lmer5)$f.ECOREGION[1] )/(VarF5 + VarCorr(lmer5)$f.ECOREGION[1] + pi^2/3)  #Proportion of variance attributed to random effects.
TotalVar<-(VarF5+VarCorr(lmer5)$f.ECOREGION[1] )/(VarF5 + VarCorr(lmer5)$f.ECOREGION[1] + pi^2/3) # Overall variance explained by fix and random
#12
fixed12 <- fixef(lmer12)[2]*model.matrix(lmer12)[,2]
varF12 <- var(fixed12)
varF12/(varF12+VarCorr(lmer12)$f.GENUS[1]+pi^2/3) #prportion of variance from fixed
(VarCorr(lmer12)$f.GENUS)[1]/(varF12 + VarCorr(lmer12)$f.GENUS[1]+pi^2/3) #prop var from random
(varF12+VarCorr(lmer12)$f.GENUS[1])/(varF12 + VarCorr(lmer12)$f.GENUS[1] + pi^2/3)  #prop var from all factors
#13
fixed13 <- fixef(lmer13)[2]*model.matrix(lmer13)[,2]+fixef(lmer13)[3]*model.matrix(lmer13)[,3]+
  fixef(lmer13)[4]*model.matrix(lmer13)[,4]+fixef(lmer13)[5]*model.matrix(lmer13)[,5]
varF13 <- var(fixed13)
FixedVar<- varF13/(varF13+VarCorr(lmer13)$f.GENUS[1]+attr(VarCorr(lmer13), "sc")^2) #prportion of variance from fixed
RandVar<- (VarCorr(lmer13)$f.GENUS)[1]/(varF13 + VarCorr(lmer13)$f.GENUS[1]+attr(VarCorr(lmer13), "sc")^2) #prop var from random
TotalVar<-(varF13+VarCorr(lmer13)$f.GENUS[1])/(varF13 + VarCorr(lmer13)$f.GENUS[1] + attr(VarCorr(lmer13), "sc")^2)  
FixedVar/TotalVar
RandVar/TotalVar
fixef1 <- var(fixef(lmer13)[2]*model.matrix(lmer13)[,2])
fixef2 <- var(fixef(lmer13)[3]*model.matrix(lmer13)[,3])
fixef3 <- var(fixef(lmer13)[4]*model.matrix(lmer13)[,4])
fixef4 <- var(fixef(lmer13)[5]*model.matrix(lmer13)[,5])
FixedVar1 <- fixef1+fixef2+fixef3+fixef4
FixedVar <- FixedVar1/(FixedVar1+VarCorr(lmer13)$f.GENUS[1]+attr(VarCorr(lmer13), "sc")^2)
f1<-fixef1/FixedVar
f2<-fixef2/FixedVar
f3<-fixef3/FixedVar
f4<-fixef4/FixedVar
f1+f2+f3+f4
#14
fixed14 <- fixef(lmer14)[2]*model.matrix(lmer14)[,2]+fixef(lmer14)[3]*model.matrix(lmer14)[,3]+
  fixef(lmer14)[4]*model.matrix(lmer14)[,4]+fixef(lmer14)[5]*model.matrix(lmer14)[,5]
varF14 <- var(fixed14)
FixedVar<- varF14/(varF14+VarCorr(lmer14)$f.GENUS[1]+pi^2/3) #prportion of variance from fixed
RandVar<- (VarCorr(lmer14)$f.GENUS)[1]/(varF14 + VarCorr(lmer14)$f.GENUS[1]+pi^2/3) #prop var from random
TotalVar<-(varF14+VarCorr(lmer14)$f.GENUS[1])/(varF14 + VarCorr(lmer14)$f.GENUS[1] + pi^2/3)  
FixedVar/TotalVar
RandVar/TotalVar
#16
fixed16 <- fixef(lmer16)[2]*model.matrix(lmer16)[,2]+fixef(lmer16)[3]*model.matrix(lmer16)[,3]+
  fixef(lmer16)[4]*model.matrix(lmer16)[,4]+fixef(lmer16)[5]*model.matrix(lmer16)[,5]
varF16 <- var(fixed16)
FixedVar<- varF16/(varF16+VarCorr(lmer16)$f.GENUS[1]+pi^2/3) #prportion of variance from fixed
rand <- VarCorr(lmer16)$f.GENUS[1]+VarCorr(lmer16)$f.ECOREGION[1]
RandVar<- (rand)/(varF16 + rand+pi^2/3) #prop var from random
TotalVar<-(varF16+rand)/(varF16 + rand + pi^2/3)  
FixedVar/TotalVar
RandVar/TotalVar
#14 again (not finished)
varNitro = VarCorr(lmer14)$Fami[1]
varRep = VarCorr(slme4)$Block[1]
h2 = 4*varFami/(varFami + varRep + 3.29)
# 12 another way
Var_Random_effect <- as.numeric(VarCorr(lmer12)) #we extract the variance of the random effect (which capture individual variation not captured by fixed effects)
Var_Residual <- attr(VarCorr(lmer12), "sc")^2 #we extract the residual variance (which corresponds to within individual variation in this case where no fixed effects differ within indiv)
Var_Fix_effect <- var(predict(lm(data+sex~sex))) #we extract the variance explained by the fixed effect
Var_Random_effect+Var_Residual+Var_Fix_effect

# Ecoregions table -----------
#requires library plyr
count(newdata, 'ECOREGION')
count(merged.soil.prod, 'ECOREGION.x')

count(newdata, 'FIXER')

# Misc ---------------
x$pred = exp(predict(lmer12))
ggplot(newdata, aes(x = lin.spp)) + geom_density() + facet_wrap(FIXER ~ ECOREGION) + labs(x="%BA in Dom Species")
summary(lm(log(newdata$SPP_O1_BAS) ~ x$pred, subset = f.GENUS)) # I need the summary statistics to compare models
ggplot(x, aes(x = time, y = y, color = Factor)) + geom_point() + geom_smooth(method = "glm", family = gaussian(lin="log"), start=c(5,0))

#not working
glm.lis <- lmList(lin.spp~c.SOILN|f.ECOREGION,data=newdata)
print(plot.lmList(glm.lis,scale=list(x=list(relation="free"))))

#plot for each genus the soilN vs lin.spp by ecoregion
print(stripplot(lin.spp ~ c.SOILN|f.GENUS,
                data=newdata,
                groups=f.ECOREGION,
                auto.key=list(lines=TRUE, columns=2),
                strip=strip.custom(strip.names=c(TRUE,TRUE)),
                type=c('p','a'), layout=c(5,2,1),
                scales=list(x=list(rot=90)), main="Genus"))
