################################################################################
##################       Data Analysis: Fiskelinie Cod Parasites      ##########
################################################################################

# Load relevant libraries
library("readxl")
library(bbmle)
library(dplyr)
library(ggplot2)
library(car)
library(ggcorrplot)
library(RColorBrewer)

# Notes:
#  - Maturity scale C
#  - Data used only from BITS Q1 - Q4 was removed since maturity staging is not 
#    done in a regular basis
#  - Baltica Q2 also removed  

# Load data extracted from database
datc <- read.csv("TOR_parasite_omr25_MaturityCSex.csv")
str(datc)

# Transform liver category to standard format
datc$liverCat <- datc$parasiteCode -3
summary(datc$liverCat)

# Modify variables
datc$liverCatF <- as.factor(datc$liverCat)
datc$maturityIndex <- as.factor(datc$maturityIndex)
datc$maturityIndexMethod <- as.factor(datc$maturityIndexMethod)
datc$dfuArea <- as.factor(datc$dfuArea)
datc$yearF <- as.factor(datc$year)
datc$cruise <- as.factor(datc$cruise)
datc$tripType <- as.factor(datc$tripType)
datc$speciesCode <- as.factor(datc$speciesCode)
datc$statisticalRectangle <- as.factor(datc$statisticalRectangle)
datc$treatment <- as.factor(datc$treatment)
datc$liverCatF <- as.factor(datc$liverCat)
datc$sexCode <- as.factor(datc$sexCode)
dim(datc)
summary(datc)
str(datc)

# Subset so that we keep only fish of >30 cm with gonad weight, liver weight, parasite code
# and BITS cruise from Q1
datc <- subset(datc, (length >= 300) & (!is.na(liverCat)) &
                  (!is.na(maturityIndex)) & (!is.na(weightLiver)) &
                  (!is.na(weightGonads)) & (!is.na(weightGutted)) & (cruise=="BITS-1"))

# Convert gutted weight to kg
datc$weightGutted <- datc$weightGutted/1000

# Compute GSI (calculated with gutted weight) and HSI (calculated with gutted weight), (see Parmeswaran et al. 1974 and DTU Aqua manual)
datc$K <- 100*(datc$weight*1000)/(datc$length/10)^3
datc$gsi <- 100*datc$weightGonads/(datc$weightGutted*1000)
datc$hsi <- 100*datc$weightLiver/(datc$weightGutted*1000)
datc$hsio <- 100*datc$weightLiver/(datc$weight*1000-datc$weightGonads)

# Select cod of maturity 4
datc1 <- subset(datc, maturityIndex==4)
nrow(datc1)
#nrow(subset(datc1, liverCat==3)) #549 fish with 0 parasites

#Claculate overall prevalence
(prevalence <- nrow(subset(datc1, liverCat>0))/nrow(datc1)) #63.76%
(prevalence <- nrow(subset(datc1, (liverCat>0 & sexCode=="F")))/nrow(subset(datc1, sexCode=="F")))
#61.72% for females
(prevalence <- nrow(subset(datc1, (liverCat>0 & sexCode=="M")))/nrow(subset(datc1, sexCode=="M")))
#65.22% males

# Percentage of cod of maturity stage 4
nrow(datc1)/nrow(datc) # 45.9% (with baltica cruise was 35,9%)

# Outliers that should be removed - Gonad weight was mistakenly reported
plot(datc1$gsi, ylab="GSI")
datc2 <- subset(datc1, gsi<50)
nrow(datc2) #lost 12 fish, all around line 400
plot(datc2$gsi)
outliers <- subset(datc1, gsi>50)
out <- outliers[,-c(1,4,5,8:17,19,20,22,24:34,36:40)] # For table in report

# Is gsi normally distributed?
par(mfrow=c(1,2))
hist(datc2$gsi[datc1$sexCode=="F"], main="Females", xlab="GSI")
hist(datc2$gsi[datc1$sexCode=="M"], main="Males", xlab="GSI")
# 640 vs 300
shapiro.test(datc2$gsi[datc2$sexCode=="F"])
shapiro.test(datc2$gsi[datc2$sexCode=="M"])
par(mfrow=c(1,2))
qqnorm(datc2$gsi[datc2$sexCode=="F"], main='Females')
qqline(datc2$gsi[datc2$sexCode=="F"])
qqnorm(datc2$gsi[datc2$sexCode=="M"], main='Males')
qqline(datc2$gsi[datc2$sexCode=="M"])
##
par(mfrow=c(1,2))
qqnorm(datc2$gsi, main='All')
qqline(datc2$gsi)
hist(datc2$gsi, main="All", xlab="GSI")
shapiro.test(datc2$gsi)

# Is there a difference in GSI between the 2 sexes?
mean(datc2$gsi[datc2$sexCode=="F"])
mean(datc2$gsi[datc2$sexCode=="M"])
sd(datc2$gsi[datc2$sexCode=="F"])
sd(datc2$gsi[datc2$sexCode=="M"])
t.test(datc2$gsi[datc2$sexCode=="M"],datc2$gsi[datc2$sexCode=="F"])
plot(datc2$gsi[datc2$sexCode=="F"], col="red")
points(datc2$gsi[datc2$sexCode=="M"], col="blue")
boxplot(datc2$gsi~datc2$sexCode)

# how is weight distributed?
par(mfrow=c(1,2))
hist(datc2$weight[datc2$sexCode=="F"], main="Females", xlab="Weight (kg)")
hist(datc2$weight[datc2$sexCode=="M"], main="Males", xlab="Weight (kg)")
# 640 vs 300

# how is gutted weight distributed?
par(mfrow=c(1,2))
hist(datc2$weightGutted[datc2$sexCode=="F"], main="Females", xlab="Gutted weight (kg)")
hist(datc2$weightGutted[datc2$sexCode=="M"], main="Males", xlab="Gutted weight (kg)")
# 640 vs 300

# how is parasite code distributed?
par(mfrow=c(1,2))
hist(datc2$liverCat[datc2$sexCode=="F"], main="Females", xlab="Parasite Code")
hist(datc2$liverCat[datc2$sexCode=="M"], main="Males", xlab="Parasite Code")
# 640 vs 300

# Checking correlation between variables
#Simplifying dataset by only keeping used columns
dcorr <- datc2[,-c(1,3:14,16:17,19,20,22:23,25:33,37,39:41)]
dcorr <- dcorr[,-c(3)]
corgg <- ggcorrplot(cor(dcorr))
#ggsave(filename = "correlation plot of dcorr.png", plot = corgg, bg="white", width = 5.3, height = 4, dpi = 700)
pairs(dcorr)


################################################################################
#######         Investigate effects of C. osculatum infection on GSI    ########
################################################################################

#Start with a complete model
zalls <- step(lm(gsi ~ (hsi+weight+liverCat+sexCode+year+K)^3, datc2), trace = FALSE)
summary(zalls)

# Reduce - backwards selection
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -weight:liverCat:year)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -weight:year:K)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -weight:year)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -sexCode:year)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -hsi:weight:K)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -hsi:liverCat:K)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -hsi:liverCat)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -hsi:K)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -liverCat:year:K)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -liverCat:year)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -liverCat:K)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -weight:liverCat:sexCode)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -liverCat:sexCode)
drop1(zalls, test = "F")
zalls <- update(zalls, ~ . -weight:sexCode:K)
drop1(zalls, test = "F")
# can't drop anything else

##### Check transformations of x variables
par(mfrow=c(1,1))
boxCox(zalls, lambda =seq(0.1, 0.6, by = 0.01))


# Try 3rd root, and also try log - do this for all variables that could have non-linear relationships
z3ralls <- lm(gsi ~ I(hsi^1/3) + weight + liverCat + sexCode + year + K +
                I(hsi^1/3):weight + I(hsi^1/3):sexCode + I(hsi^1/3):year + weight:liverCat + 
                weight:sexCode + weight:K + year:K + I(hsi^1/3):weight:sexCode, data = datc2)

z3ralls2 <- lm(gsi ~ hsi + I(weight^1/3) + liverCat + sexCode + year + K +
                hsi:I(weight^1/3) + hsi:sexCode + hsi:year + I(weight^1/3):liverCat + 
                I(weight^1/3):sexCode + I(weight^1/3):K + year:K + hsi:I(weight^1/3):sexCode, data = datc2)

z3ralls3 <- lm(gsi ~ hsi + weight + I(liverCat^1/3) + sexCode + year + K +
                 hsi:weight + hsi:sexCode + hsi:year + weight:I(liverCat^1/3) + 
                 weight:sexCode + weight:K + year:K + hsi:weight:sexCode, data = datc2)

z3ralls4 <- lm(gsi ~ hsi + weight + liverCat + sexCode + year + I(K^1/3) +
                 hsi:weight + hsi:sexCode + hsi:year + weight:liverCat + 
                 weight:sexCode + weight:I(K^1/3) + year:I(K^1/3) + hsi:weight:sexCode, data = datc2)

zlogalls <- lm(gsi ~ log(hsi) + weight + liverCat + sexCode + year + K +
                 log(hsi):weight + log(hsi):sexCode + log(hsi):year + weight:liverCat + 
                 weight:sexCode + weight:K + year:K + log(hsi):weight:sexCode, data = datc2)
                 
zlogalls2 <- lm(gsi ~ hsi + log(weight) + liverCat + sexCode + year + K +
                  hsi:log(weight) + hsi:sexCode + hsi:year + log(weight):liverCat + 
                  log(weight):sexCode + log(weight):K + year:K + hsi:log(weight):sexCode, data = datc2)
  
zlogalls3 <- lm(gsi ~ log(hsi) + log(weight) + liverCat + sexCode + year + K +
                  log(hsi):log(weight) + log(hsi):sexCode + log(hsi):year + log(weight):liverCat + 
                  log(weight):sexCode + log(weight):K + year:K + log(hsi):log(weight):sexCode, data = datc2)

zlogalls4 <- lm(gsi ~ log(hsi) + log(weight) + log(liverCat) + sexCode + year + K +
                  log(hsi):log(weight) + log(hsi):sexCode + log(hsi):year + log(weight):log(liverCat) + 
                  log(weight):sexCode + log(weight):K + year:K + log(hsi):log(weight):sexCode, data = datc2)

zlogalls5 <- lm(gsi ~ log(hsi) + weight + log(liverCat) + sexCode + year + K +
                  log(hsi):weight + log(hsi):sexCode + log(hsi):year + weight:log(liverCat) + 
                  weight:sexCode + weight:K + year:K + log(hsi):weight:sexCode, data = datc2)

zlogalls6 <- lm(gsi ~ hsi + log(weight) + log(liverCat) + sexCode + year + K +
                  hsi:log(weight) + hsi:sexCode + hsi:year + log(weight):log(liverCat) + 
                  log(weight):sexCode + log(weight):K + year:K + hsi:log(weight):sexCode, data = datc2)

zlogalls7 <- lm(gsi ~ hsi + weight + log(liverCat) + sexCode + year + K +
                  hsi:weight + hsi:sexCode + hsi:year + weight:log(liverCat) + 
                  weight:sexCode + weight:K + year:K + hsi:weight:sexCode, data = datc2)

zlogalls8 <- lm(gsi ~ hsi + weight + liverCat + sexCode + year + log(K) +
                  hsi:weight + hsi:sexCode + hsi:year + weight:liverCat + 
                  weight:sexCode + weight:log(K) + year:log(K) + hsi:weight:sexCode, data = datc2)

zlogalls9 <- lm(gsi ~ log(hsi) + weight + liverCat + sexCode + year + log(K) +
                  log(hsi):weight + log(hsi):sexCode + log(hsi):year + weight:liverCat + 
                  weight:sexCode + weight:log(K) + year:log(K) + log(hsi):weight:sexCode, data = datc2)


# Perform an AIC test to select most robust model
AICtab(z3ralls, z3ralls2, z3ralls3, z3ralls4, zlogalls, zlogalls2, zlogalls3, zlogalls4, zlogalls5, zlogalls6, zlogalls7, zlogalls8, zlogalls9)
summary(zlogalls) # Log is best transformation

# Simplify model again
drop1(zlogalls, test = "F")
zlogallst <- update(zlogalls, ~ . -weight:K)
drop1(zlogallst, test = "F")
zlogallst <- update(zlogallst, ~ . -log(hsi):year)
drop1(zlogallst, test = "F")


### Restart model selection with log transformation to check if the same model is obtained
final <- step(lm(gsi ~ (log(hsi)+weight+K+liverCat+sexCode+year)^3, datc2), trace = FALSE)
final
summary(final)
drop1(final, test = "F")
final <- update(final, ~ . -K:liverCat:year)
drop1(final, test = "F")
final <- update(final, ~ . -liverCat:year)
drop1(final, test = "F")
final <- update(final, ~ . -weight:liverCat:sexCode)
drop1(final, test = "F")
final <- update(final, ~ . -K:liverCat:sexCode)
drop1(final, test = "F")
final <- update(final, ~ . -liverCat:sexCode)
drop1(final, test = "F")
final <- update(final, ~ . -K:liverCat)
drop1(final, test = "F")
final <- update(final, ~ . -weight:K:sexCode)
drop1(final, test = "F")
summary(final)
#This is the same model as the zlogalls, minus the log(hsi):year
# in fact it became the same model as when i was using different GSI (from total not gutted weight)

# Perform AIC
AICtab(final,zlogalls,z3ralls,zalls,zlogallst)
#final is most robust model

# Check final
summary(final)
par(mfrow=c(2,2))
plot(final, col=datc2$sexCode) # Look good - save 500 x 500 res

# Check fistributiuon of residuals
legend("bottomright", legend = unique(datc2$sexCode), fill = unique(datc2$sexCode), cex=0.5)
par(mfrow=c(1,1))
hist(resid(final))#normal
shapiro.test(resid(final))
# 300 x 250


################################################################################
##################                    PLOT PREDICTION                 ##########
################################################################################

# Plot predictions - HSI for sexes and liver category high and low
new_flowp <- expand.grid("sexCode"= c("F"), "year"= c(2021), "liverCat" = c(0), "weight" = mean(datc2$weight[datc2$sexCode=="F"]), "K" = mean(datc2$K[datc2$sexCode=="F"]), "hsi" = seq(min(datc2$hsi[datc2$sexCode=="F"]),max(datc2$hsi[datc2$sexCode=="F"]), length = 1000))
new_fhighp <- expand.grid("sexCode"= c("F"), "year"= c(2021), "liverCat" = c(4), "weight" = mean(datc2$weight[datc2$sexCode=="F"]), "K" = mean(datc2$K[datc2$sexCode=="F"]), "hsi" = seq(min(datc2$hsi[datc2$sexCode=="F"]),max(datc2$hsi[datc2$sexCode=="F"]), length = 1000))
new_mlowp <- expand.grid("sexCode"= c("M"), "year"= c(2021), "liverCat" = c(0), "weight" = mean(datc2$weight[datc2$sexCode=="M"]), "K" = mean(datc2$K[datc2$sexCode=="M"]), "hsi" = seq(min(datc2$hsi[datc2$sexCode=="M"]),max(datc2$hsi[datc2$sexCode=="M"]), length = 1000))
new_mhighp <- expand.grid("sexCode"= c("M"), "year"= c(2021), "liverCat" = c(4), "weight" = mean(datc2$weight[datc2$sexCode=="M"]), "K" = mean(datc2$K[datc2$sexCode=="M"]), "hsi" = seq(min(datc2$hsi[datc2$sexCode=="M"]),max(datc2$hsi[datc2$sexCode=="M"]), length = 1000))

new_flowp$sexCode <- as.factor(new_flowp$sexCode)
new_fhighp$sexCode <- as.factor(new_fhighp$sexCode)
new_mlowp$sexCode <- as.factor(new_mlowp$sexCode)
new_mhighp$sexCode <- as.factor(new_mhighp$sexCode)


c_flowp <- predict(final, new_flowp, interval = "confidence")
c_fhighp <- predict(final, new_fhighp, interval = "confidence")
c_mlowp <- predict(final, new_mlowp, interval = "confidence")
c_mhighp <- predict(final, new_mhighp, interval = "confidence")


par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$hsi[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "HSI", ylab = "GSI")
matlines(new_flowp$hsi, c_flowp[,1:3],
         lty = c(1,2,2), lw = 2, col = "blue")
matlines(new_fhighp$hsi, c_fhighp[,1:3],
         lty = c(1,2,2), lw = 2, col = "red")
legend("topright", legend = c("Liver Category = 0","Liver Category = 4"), col = c("blue","red"), lty = 1, lwd=2, cex=0.6)

plot(datc2$hsi[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = "HSI", ylab = "")
matlines(new_mlowp$hsi, c_mlowp[,1:3],
         lty = c(1,2,2), lw = 2, col = "blue")
matlines(new_mhighp$hsi, c_mhighp[,1:3],
         lty = c(1,2,2), lw = 2, col = "red")
# 800 x 500

# Plot predictions - Liver categories for sexes and HSI conditions
new_flowp2 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = 0:4, "weight" = mean(datc2$weight[datc2$sexCode=="F"]), "hsi" = 5.569) # summary(datc2$hsi[datc2$sexCode=="F"])
new_fhighp2 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = 0:4, "weight" = mean(datc2$weight[datc2$sexCode=="F"]), "hsi" = 8.354)
new_mlowp2 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = 0:4, "weight" = mean(datc2$weight[datc2$sexCode=="M"]), "hsi" = 4.168)
new_mhighp2 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = 0:4, "weight" = mean(datc2$weight[datc2$sexCode=="M"]), "hsi" = 6.508)

new_flowp2$sexCode <- as.factor(new_flowp2$sexCode)
new_fhighp2$sexCode <- as.factor(new_fhighp2$sexCode)
new_mlowp2$sexCode <- as.factor(new_mlowp2$sexCode)
new_mhighp2$sexCode <- as.factor(new_mhighp2$sexCode)

c_flowp2 <- predict(final, new_flowp2, interval = "confidence")
c_fhighp2 <- predict(final, new_fhighp2, interval = "confidence")
c_mlowp2 <- predict(final, new_mlowp2, interval = "confidence")
c_mhighp2 <- predict(final, new_mhighp2, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$liverCat[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "Liver Category", ylab = "GSI")
matlines(new_flowp2$liverCat, c_flowp2[,1:3],
         lty = c(1,3,3), lw = 2, col = "gold2")
matlines(new_fhighp2$liverCat, c_fhighp2[,1:3],
         lty = c(1,3,3), lw = 2, col = "darkolivegreen4")
legend(x = 4.2, y = 33.5, legend = c("HSI = 1st Qu.","HSI = 3rd Qu."), col = c("gold2","darkolivegreen4"), lty = 1, lwd=2, cex=0.6, bty="n")

plot(datc2$liverCat[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = "Liver Category", ylab = "")
matlines(new_mlowp2$liverCat, c_mlowp2[,1:3],
         lty = c(1,3,3), lw = 2, col = "gold2")
matlines(new_mhighp2$liverCat, c_mhighp2[,1:3],
         lty = c(1,3,3), lw = 2, col = "darkolivegreen4")
# 800 x 500


# Plot predictions - year (time) for sexes and liver category high and low
new_flowp3 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= 2017:2023, "liverCat" = c(0), "weight" = mean(datc2$weight[datc2$sexCode=="F"]), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
new_fhighp3 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= 2017:2023, "liverCat" = c(4), "weight" = mean(datc2$weight[datc2$sexCode=="F"]), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
new_mlowp3 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= 2017:2023, "liverCat" = c(0), "weight" = mean(datc2$weight[datc2$sexCode=="M"]), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
new_mhighp3 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= 2017:2023, "liverCat" = c(4), "weight" = mean(datc2$weight[datc2$sexCode=="M"]), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))

new_flowp3$sexCode <- as.factor(new_flowp3$sexCode)
new_fhighp3$sexCode <- as.factor(new_fhighp3$sexCode)
new_mlowp3$sexCode <- as.factor(new_mlowp3$sexCode)
new_mhighp3$sexCode <- as.factor(new_mhighp3$sexCode)

c_flowp3 <- predict(final, new_flowp3, interval = "confidence")
c_fhighp3 <- predict(final, new_fhighp3, interval = "confidence")
c_mlowp3 <- predict(final, new_mlowp3, interval = "confidence")
c_mhighp3 <- predict(final, new_mhighp3, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$year[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "Year", ylab = "GSI")
matlines(new_flowp3$year, c_flowp3[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightblue2")
matlines(new_fhighp3$year, c_fhighp3[,1:3],
         lty = c(1,3,3), lw = 2, col = "purple3")

plot(datc2$year[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = "Year", ylab = "")
matlines(new_mlowp3$year, c_mlowp3[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightblue2")
matlines(new_mhighp3$year, c_mhighp3[,1:3],
         lty = c(1,3,3), lw = 2, col = "purple3")
legend(x = 2020.5, y = 37, legend = c("Liver category = 0","Liver category = 4"), col = c("lightblue2","purple3"), lty = 1, lwd=2, cex=0.6, bty="n")

# 800 x 500


# Plot predictions - weight for sexes and liver category high and low
new_flowp4 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = c(0), "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
new_fhighp4 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = c(4), "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
new_mlowp4 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = c(0), "weight" = seq(min(datc2$weight[datc2$sexCode=="M"]),max(datc2$weight[datc2$sexCode=="M"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
new_mhighp4 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = c(4), "weight" = seq(min(datc2$weight[datc2$sexCode=="M"]),max(datc2$weight[datc2$sexCode=="M"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))

new_flowp4$sexCode <- as.factor(new_flowp4$sexCode)
new_fhighp4$sexCode <- as.factor(new_fhighp4$sexCode)
new_mlowp4$sexCode <- as.factor(new_mlowp4$sexCode)
new_mhighp4$sexCode <- as.factor(new_mhighp4$sexCode)

c_flowp4 <- predict(final, new_flowp4, interval = "confidence")
c_fhighp4 <- predict(final, new_fhighp4, interval = "confidence")
c_mlowp4 <- predict(final, new_mlowp4, interval = "confidence")
c_mhighp4 <- predict(final, new_mhighp4, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$weight[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "Weight (kg)", ylab = "GSI")
matlines(new_flowp4$weight, c_flowp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "orange2")
matlines(new_fhighp4$weight, c_fhighp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "grey")

plot(datc2$weight[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = "Weight (kg)", ylab = "")
matlines(new_mlowp4$weight, c_mlowp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "orange2")
matlines(new_mhighp4$weight, c_mhighp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "grey")
legend(x = 1, y = 35, legend = c("Liver category = 0","Liver category = 4"), col = c("orange2","grey"), lty = 1, lwd=2, cex=0.6, bty="n")

# 800 x 500

#### GGPLOT VERSION
mpp <- data.frame("weight"=new_flowp4$weight, "mpp"=c_flowp4[,1], "mppmin"=c_flowp4[,2], "mppmax"=c_flowp4[,3])
mpp2 <- data.frame("weight"=new_fhighp4$weight, "mpp"=c_fhighp4[,1], "mppmin"=c_fhighp4[,2], "mppmax"=c_fhighp4[,3])
at4 <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode=="F",], aes(x = weight, y = gsi),
             size = 2, color="black", shape = 19, alpha=0.5) +
  geom_ribbon(data = mpp, aes(x =weight, ymin = mppmin, ymax = mppmax), alpha=0.5, fill="#EC9706") + ##FA8128
  geom_line(data = mpp, aes(x = weight, y= mpp), color="orange2", lwd=1.2) +
  geom_ribbon(data = mpp2, aes(x =weight, ymin = mppmin, ymax = mppmax), alpha=0.6, fill="lightgrey") +
  geom_line(data = mpp2, aes(x = weight, y= mpp), color="grey70", lwd=1.2) +
  labs(y = "GSI", x = "Weight (kg)", title="Females") +
  theme_classic()
#ggsave(filename = "GSI vs weight for parasites ggplot females.png", plot = at4, bg="white", width = 4, height = 4, dpi = 500)
#with legend for same proportions
at4fl <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode=="F",], aes(x = weight, y = gsi, color = "Observations"),
             size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = mpp, aes(x =weight, ymin = mppmin, ymax = mppmax, fill = "Not infected"), alpha=0.5) + ##FA8128
  geom_line(data = mpp, aes(x = weight, y= mpp, color = "Not infected"), lwd=1.2) +
  geom_ribbon(data = mpp2, aes(x =weight, ymin = mppmin, ymax = mppmax, fill = "Liver category = 4"), alpha=0.6) +
  geom_line(data = mpp2, aes(x = weight, y= mpp, color = "Liver category = 4"), lwd=1.2) +
  labs(y = "GSI", x = "Weight (kg)", title="Females") +
  scale_color_manual(name = "Legend", values = c("Observations" = "black", "Not infected" = "orange2", "Liver category = 4" = "grey70")) +
  scale_fill_manual(name = "Legend", values = c("Not infected" = "#EC9706", "Liver category = 4" = "lightgrey"), guide = "none") +
  theme(legend.position = "bottom") +
  theme_classic()
#ggsave(filename = "GSI vs weight for parasites ggplot females legend.png", plot = at4fl, bg="white", width = 5.3, height = 4, dpi = 500)


mppm <- data.frame("weight"=new_mlowp4$weight, "mpp"=c_mlowp4[,1], "mppmin"=c_mlowp4[,2], "mppmax"=c_mlowp4[,3])
mppm2 <- data.frame("weight"=new_mhighp4$weight, "mpp"=c_mhighp4[,1], "mppmin"=c_mhighp4[,2], "mppmax"=c_mhighp4[,3])
at4m <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode == "M", ], aes(x = weight, y = gsi, color = "Observations"), #datc2[(datc2$sexCode == "M" & datc2$liverCat %in% c(3, 7)), ] to only pic parasite codes of 3 and 7
             size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = mppm, aes(x =weight, ymin = mppmin, ymax = mppmax, fill = "Not infected"), alpha=0.5) +
  geom_line(data = mppm, aes(x = weight, y= mpp, color = "Not infected"), lwd=1.2) +
  geom_ribbon(data = mppm2, aes(x =weight, ymin = mppmin, ymax = mppmax, fill = "Liver category = 4"), alpha=0.6) +
  geom_line(data = mppm2, aes(x = weight, y= mpp, color = "Liver category = 4"), lwd=1.2) +
  labs(y = "GSI", x = "Weight (kg)", title="Males") +
  scale_color_manual(name = "Legend", values = c("Observations" = "black", "Not infected" = "orange2", "Liver category = 4" = "grey70")) +
  scale_fill_manual(name = "Legend", values = c("Not infected" = "#EC9706", "Liver category = 4" = "lightgrey"), guide = "none") +
  theme(legend.position = "bottom") +
  theme_classic()
#ggsave(filename = "GSI vs weight for parasites ggplot males.png", plot = at4m, bg="white", width = 5.3, height = 4, dpi = 500)


# # Plot predictions - HSI for sexes and liver category high and low - althernative
par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$weight[datc2$sexCode=="F" & datc2$liverCat!=0 & datc2$liverCat!=4], datc2$gsi[datc2$sexCode=="F"& datc2$liverCat!=0 & datc2$liverCat!=4], col="white",pch=1, main="Females", xlab = "Weight (kg)", ylab = "GSI")
points (datc2$weight[datc2$sexCode=="F" & datc2$liverCat==0], datc2$gsi[datc2$sexCode=="F"& datc2$liverCat==0], col = "orange2",pch=1)
points (datc2$weight[datc2$sexCode=="F" & datc2$liverCat==4], datc2$gsi[datc2$sexCode=="F"& datc2$liverCat==4], col = "grey", pch=1)
matlines(new_flowp4$weight, c_flowp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "orange2")
matlines(new_fhighp4$weight, c_fhighp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "grey")

plot(datc2$weight[datc2$sexCode=="M"& datc2$liverCat!=0 & datc2$liverCat!=4], datc2$gsi[datc2$sexCode=="M"& datc2$liverCat!=0 & datc2$liverCat!=4], col="white",pch=1, main="Males", xlab = "Weight (kg)", ylab = "")
points (datc2$weight[datc2$sexCode=="M" & datc2$liverCat==0], datc2$gsi[datc2$sexCode=="M"& datc2$liverCat==0], col = "orange2")
points (datc2$weight[datc2$sexCode=="M" & datc2$liverCat==4], datc2$gsi[datc2$sexCode=="M"& datc2$liverCat==4], col = "grey")
matlines(new_mlowp4$weight, c_mlowp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "orange2")
matlines(new_mhighp4$weight, c_mhighp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "grey")
legend(x = 1, y = 35, legend = c("Parasite Code = 3","Parasite Code = 7"), col = c("orange2","grey"), lty = 1, lwd=2, cex=0.6, bty="n")

# 800 x 500

# # Plot predictions - HSI for sexes and liver category high and low - liver cat of 3 instead of 4
new_flowp4 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = c(0), "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
new_fhighp4 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = c(3), "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
new_mlowp4 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = c(0), "weight" = seq(min(datc2$weight[datc2$sexCode=="M"]),max(datc2$weight[datc2$sexCode=="M"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
new_mhighp4 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = c(3), "weight" = seq(min(datc2$weight[datc2$sexCode=="M"]),max(datc2$weight[datc2$sexCode=="M"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))

new_flowp4$sexCode <- as.factor(new_flowp4$sexCode)
new_fhighp4$sexCode <- as.factor(new_fhighp4$sexCode)
new_mlowp4$sexCode <- as.factor(new_mlowp4$sexCode)
new_mhighp4$sexCode <- as.factor(new_mhighp4$sexCode)

c_flowp4 <- predict(final, new_flowp4, interval = "confidence")
c_fhighp4 <- predict(final, new_fhighp4, interval = "confidence")
c_mlowp4 <- predict(final, new_mlowp4, interval = "confidence")
c_mhighp4 <- predict(final, new_mhighp4, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$weight[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "Weight (kg)", ylab = "GSI")
matlines(new_flowp4$weight, c_flowp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "orange2")
matlines(new_fhighp4$weight, c_fhighp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "grey")

plot(datc2$weight[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = "Weight (kg)", ylab = "")
matlines(new_mlowp4$weight, c_mlowp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "orange2")
matlines(new_mhighp4$weight, c_mhighp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "grey")
legend(x = 1, y = 35, legend = c("Liver category = 3","Liver category = 6"), col = c("orange2","grey"), lty = 1, lwd=2, cex=0.6, bty="n")

# 800 x 500

# # Plot predictions - liver cat for sexes and weights
new_flowp5 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = 0:4, "weight" = 0.3, "hsi" = mean(datc2$hsi[datc2$sexCode=="F"])) # summary(datc2$hsi[datc2$sexCode=="F"])
new_fhighp5 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = 0:4, "weight" = 1.25, "hsi" = mean(datc2$hsi[datc2$sexCode=="F"])) # also tried close to max
new_mlowp5 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = 0:4, "weight" = 0.3, "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
new_mhighp5 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = 0:4, "weight" = 1.25, "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))

new_flowp5$sexCode <- as.factor(new_flowp5$sexCode)
new_fhighp5$sexCode <- as.factor(new_fhighp5$sexCode)
new_mlowp5$sexCode <- as.factor(new_mlowp5$sexCode)
new_mhighp5$sexCode <- as.factor(new_mhighp5$sexCode)

c_flowp5 <- predict(final, new_flowp5, interval = "confidence")
c_fhighp5 <- predict(final, new_fhighp5, interval = "confidence")
c_mlowp5 <- predict(final, new_mlowp5, interval = "confidence")
c_mhighp5 <- predict(final, new_mhighp5, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$liverCat[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "Liver category", ylab = "GSI")
matlines(new_flowp5$liverCat, c_flowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightgreen")
matlines(new_fhighp5$liverCat, c_fhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "darkred")
legend(x = 4.2, y = 33.5, legend = c("Weight = 0.3 kg.","Weight = 1.25 kg"), col = c("lightgreen","darkred"), lty = 1, lwd=2, cex=0.6, bty="n")

plot(datc2$liverCat[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = "Liver category", ylab = "")
matlines(new_mlowp5$liverCat, c_mlowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightgreen")
matlines(new_mhighp5$liverCat, c_mhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "darkred")
# 800 x 500
mpp1 <- data.frame("liverCat"=new_flowp5$liverCat, "mpp"=c_flowp5[,1], "mppmin"=c_flowp5[,2], "mppmax"=c_flowp5[,3])
mpp12 <- data.frame("liverCat"=new_fhighp5$liverCat, "mpp"=c_fhighp5[,1], "mppmin"=c_fhighp5[,2], "mppmax"=c_fhighp5[,3])

at1fl <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode=="F",], aes(x = liverCat, y = gsi, color = "Observations"),
             size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = mpp1, aes(x = liverCat, ymin = mppmin, ymax = mppmax, fill = "Weight = 0.3 kg"), alpha=0.5) + 
  geom_line(data = mpp1, aes(x = liverCat, y= mpp, color = "Weight = 0.3 kg"), lwd=1.2) +
  geom_ribbon(data = mpp12, aes(x = liverCat, ymin = mppmin, ymax = mppmax, fill = "Weight = 1.25 kg"), alpha=0.6) +
  geom_line(data = mpp12, aes(x = liverCat, y= mpp, color = "Weight = 1.25 kg"), lwd=1.2) +
  labs(y = "GSI", x = "Liver category", title="Females") +
  scale_color_manual(name = "Legend", values = c("Observations" = "black", "Weight = 0.3 kg" = "lightgreen", "Weight = 1.25 kg" = "darkred")) +
  scale_fill_manual(name = "Legend", values = c("Weight = 0.3 kg" = "lightgreen", "Weight = 1.25 kg" = "darkred"), guide = "none") +
  theme(legend.position = "bottom") +
  theme_classic()
at1fl
#ggsave(filename = "GSI vs liverCat for weights ggplot females legend.png", plot = at1fl, bg="white", width = 5.3, height = 4, dpi = 500)


mppm1 <- data.frame("liverCat"=new_mlowp5$liverCat, "mpp"=c_mlowp5[,1], "mppmin"=c_mlowp5[,2], "mppmax"=c_mlowp5[,3])
mppm12 <- data.frame("liverCat"=new_mhighp5$liverCat, "mpp"=c_mhighp5[,1], "mppmin"=c_mhighp5[,2], "mppmax"=c_mhighp5[,3])
at1m <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode=="M",], aes(x = liverCat, y = gsi, color = "Observations"),
             size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = mppm1, aes(x =liverCat, ymin = mppmin, ymax = mppmax, fill = "Weight = 0.3 kg"), alpha=0.5) +
  geom_line(data = mppm1, aes(x = liverCat, y= mpp, color = "Weight = 0.3 kg"), lwd=1.2) +
  geom_ribbon(data = mppm12, aes(x =liverCat, ymin = mppmin, ymax = mppmax, fill = "Weight = 1.25 kg"), alpha=0.6) +
  geom_line(data = mppm12, aes(x = liverCat, y= mpp, color = "Weight = 1.25 kg"), lwd=1.2) +
  labs(y = "GSI", x = "Liver category", title="Males") +
  scale_color_manual(name = "Legend", values = c("Observations" = "black", "Weight = 0.3 kg" = "lightgreen", "Weight = 1.25 kg" = "darkred")) +
  scale_fill_manual(name = "Legend", values = c("Weight = 0.3 kg" = "lightgreen", "Weight = 1.25 kg" = "darkred"), guide = "none") +
  theme(legend.position = "bottom") +
  theme_classic()
at1m
#ggsave(filename = "GSI vs liverCat for weights ggplot males.png", plot = at1m, bg="white", width = 5.3, height = 4, dpi = 500)



# # Plot predictions - liver category for sexes and weights
new_flowp5 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = 0:4, "weight" = 0.4738, "hsi" = mean(datc2$hsi[datc2$sexCode=="F"])) # summary(datc2$hsi[datc2$sexCode=="F"])
new_fhighp5 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = 0:4, "weight" = 0.8790, "hsi" = mean(datc2$hsi[datc2$sexCode=="F"])) # also tried close to max
new_mlowp5 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = 0:4, "weight" = 0.3420, "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
new_mhighp5 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = 0:4, "weight" = 0.6025, "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))

new_flowp5$sexCode <- as.factor(new_flowp5$sexCode)
new_fhighp5$sexCode <- as.factor(new_fhighp5$sexCode)
new_mlowp5$sexCode <- as.factor(new_mlowp5$sexCode)
new_mhighp5$sexCode <- as.factor(new_mhighp5$sexCode)

c_flowp5 <- predict(final, new_flowp5, interval = "confidence")
c_fhighp5 <- predict(final, new_fhighp5, interval = "confidence")
c_mlowp5 <- predict(final, new_mlowp5, interval = "confidence")
c_mhighp5 <- predict(final, new_mhighp5, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$liverCat[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "Liver category", ylab = "GSI")
matlines(new_flowp5$liverCat, c_flowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightgreen")
matlines(new_fhighp5$liverCat, c_fhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "darkred")
legend(x = 4.2, y = 33.5, legend = c("Weight = 1st qt.","Weight = 3rd qt."), col = c("lightgreen","darkred"), lty = 1, lwd=2, cex=0.6, bty="n")

plot(datc2$liverCat[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = "Liver category", ylab = "")
matlines(new_mlowp5$liverCat, c_mlowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightgreen")
matlines(new_mhighp5$liverCat, c_mhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "darkred")
# 800 x 500

# # Plot predictions - K for sexes and liver category high and low
new_flowp5 <- expand.grid("sexCode"= c("F"), "K" = seq(min(datc2$K[datc2$sexCode=="F"]), max(datc2$K[datc2$sexCode=="F"]), length=1000), "year"= c(2021), "liverCat" = 0, "weight" = mean(datc2$weight[datc2$sexCode=="F"]), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
new_fhighp5 <- expand.grid("sexCode"= c("F"), "K" = seq(min(datc2$K[datc2$sexCode=="F"]), max(datc2$K[datc2$sexCode=="F"]), length=1000), "year"= c(2021), "liverCat" = 4, "weight" = mean(datc2$weight[datc2$sexCode=="F"]), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"])) 
new_mlowp5 <- expand.grid("sexCode"= c("M"), "K" = seq(min(datc2$K[datc2$sexCode=="M"]), max(datc2$K[datc2$sexCode=="M"]), length=1000), "year"= c(2021), "liverCat" = 0, "weight" = mean(datc2$weight[datc2$sexCode=="M"]), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
new_mhighp5 <- expand.grid("sexCode"= c("M"), "K" = seq(min(datc2$K[datc2$sexCode=="M"]), max(datc2$K[datc2$sexCode=="M"]), length=1000), "year"= c(2021), "liverCat" = 4, "weight" = mean(datc2$weight[datc2$sexCode=="M"]), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))

new_flowp5$sexCode <- as.factor(new_flowp5$sexCode)
new_fhighp5$sexCode <- as.factor(new_fhighp5$sexCode)
new_mlowp5$sexCode <- as.factor(new_mlowp5$sexCode)
new_mhighp5$sexCode <- as.factor(new_mhighp5$sexCode)

c_flowp5 <- predict(final, new_flowp5, interval = "confidence")
c_fhighp5 <- predict(final, new_fhighp5, interval = "confidence")
c_mlowp5 <- predict(final, new_mlowp5, interval = "confidence")
c_mhighp5 <- predict(final, new_mhighp5, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$K[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "Fulton's K", ylab = "GSI")
matlines(new_flowp5$K, c_flowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightblue4")
matlines(new_fhighp5$K, c_fhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "yellow3")

plot(datc2$K[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = "Fulton's K", ylab = "")
matlines(new_mlowp5$K, c_mlowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightblue4")
matlines(new_mhighp5$K, c_mhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "yellow3")
legend(x = 0.7, y = 35, legend = c("Liver category = 0","Liver category = 4"), col = c("lightblue4","yellow3"), lty = 1, lwd=2, cex=0.6, bty="n")

# 800 x 500

# GGPLOT VERSION
# Plot predictions - K for sexes and liver category high and low
mpp5 <- data.frame("K"=new_flowp5$K, "mpp"=c_flowp5[,1], "mppmin"=c_flowp5[,2], "mppmax"=c_flowp5[,3])
mpp52 <- data.frame("K"=new_fhighp5$K, "mpp"=c_fhighp5[,1], "mppmin"=c_fhighp5[,2], "mppmax"=c_fhighp5[,3])

at5fl <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode=="F",], aes(x = K, y = gsi, color = "Observations"),
             size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = mpp5, aes(x = K, ymin = mppmin, ymax = mppmax, fill = "Not infected"), alpha=0.5) + 
  geom_line(data = mpp5, aes(x = K, y= mpp, color = "Not infected"), lwd=1.2) +
  geom_ribbon(data = mpp52, aes(x = K, ymin = mppmin, ymax = mppmax, fill = "Liver category = 4"), alpha=0.6) +
  geom_line(data = mpp52, aes(x = K, y= mpp, color = "Liver category = 4"), lwd=1.2) +
  labs(y = "GSI", x = "Fulton's K", title="Females") +
  scale_color_manual(name = "Legend", values = c("Observations" = "black", "Not infected" = "lightblue4", "Liver category = 4" = "yellow3")) +
  scale_fill_manual(name = "Legend", values = c("Not infected" = "lightblue4", "Liver category = 4" = "yellow3"), guide = "none") +
  theme(legend.position = "bottom") +
  theme_classic()
at5fl
#ggsave(filename = "GSI vs k for parasites ggplot females legend.png", plot = at5fl, bg="white", width = 5.3, height = 4, dpi = 500)


mppm5 <- data.frame("K"=new_mlowp5$K, "mpp"=c_mlowp5[,1], "mppmin"=c_mlowp5[,2], "mppmax"=c_mlowp5[,3])
mppm52 <- data.frame("K"=new_mhighp5$K, "mpp"=c_mhighp5[,1], "mppmin"=c_mhighp5[,2], "mppmax"=c_mhighp5[,3])
at5m <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode=="M",], aes(x = K, y = gsi, color = "Observations"),
             size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = mppm5, aes(x =K, ymin = mppmin, ymax = mppmax, fill = "Not infected"), alpha=0.5) +
  geom_line(data = mppm5, aes(x = K, y= mpp, color = "Not infected"), lwd=1.2) +
  geom_ribbon(data = mppm52, aes(x =K, ymin = mppmin, ymax = mppmax, fill = "Liver category = 4"), alpha=0.6) +
  geom_line(data = mppm52, aes(x = K, y= mpp, color = "Liver category = 4"), lwd=1.2) +
  labs(y = "GSI", x = "Fulton's K", title="Males") +
  scale_color_manual(name = "Legend", values = c("Observations" = "black", "Not infected" = "lightblue4", "Liver category = 4" = "yellow3")) +
  scale_fill_manual(name = "Legend", values = c("Not infected" = "lightblue4", "Liver category = 4" = "yellow3"), guide = "none") +
  theme(legend.position = "bottom") +
  theme_classic()
at5m
#ggsave(filename = "GSI vs k for parasites ggplot males.png", plot = at5m, bg="white", width = 5.3, height = 4, dpi = 500)

# Plot predictions - K for sexes and liver category high and low -  but for large fish (1 kg)
new_flowp5 <- expand.grid("sexCode"= c("F"), "K" = seq(min(datc2$K[datc2$sexCode=="F"]), max(datc2$K[datc2$sexCode=="F"]), length=1000), "year"= c(2021), "liverCat" = 0, "weight" = 1.4, "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
new_fhighp5 <- expand.grid("sexCode"= c("F"), "K" = seq(min(datc2$K[datc2$sexCode=="F"]), max(datc2$K[datc2$sexCode=="F"]), length=1000), "year"= c(2021), "liverCat" = 4, "weight" = 1.4, "hsi" = mean(datc2$hsi[datc2$sexCode=="F"])) 
new_mlowp5 <- expand.grid("sexCode"= c("M"), "K" = seq(min(datc2$K[datc2$sexCode=="M"]), max(datc2$K[datc2$sexCode=="M"]), length=1000), "year"= c(2021), "liverCat" = 0, "weight" = 1.1, "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
new_mhighp5 <- expand.grid("sexCode"= c("M"), "K" = seq(min(datc2$K[datc2$sexCode=="M"]), max(datc2$K[datc2$sexCode=="M"]), length=1000), "year"= c(2021), "liverCat" = 4, "weight" = 1.1, "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))

new_flowp5$sexCode <- as.factor(new_flowp5$sexCode)
new_fhighp5$sexCode <- as.factor(new_fhighp5$sexCode)
new_mlowp5$sexCode <- as.factor(new_mlowp5$sexCode)
new_mhighp5$sexCode <- as.factor(new_mhighp5$sexCode)

c_flowp5 <- predict(final, new_flowp5, interval = "confidence")
c_fhighp5 <- predict(final, new_fhighp5, interval = "confidence")
c_mlowp5 <- predict(final, new_mlowp5, interval = "confidence")
c_mhighp5 <- predict(final, new_mhighp5, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$K[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females (1.4 kg)", xlab = "Fulton's K", ylab = "GSI")
matlines(new_flowp5$K, c_flowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightblue4")
matlines(new_fhighp5$K, c_fhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "yellow3")

plot(datc2$K[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males (1.1 kg)", xlab = "Fulton's K", ylab = "")
matlines(new_mlowp5$K, c_mlowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightblue4")
matlines(new_mhighp5$K, c_mhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "yellow3")
legend(x = 0.7, y = 35, legend = c("Liver catgory = 0","Liver catgor = 4"), col = c("lightblue4","yellow3"), lty = 1, lwd=2, cex=0.6, bty="n")

# 800 x 500

################################################################################
############         PLOT ALLPARASITE CATEGORIES : female fish
################################################################################

names_ndpcf <- vector("list", 5)
cnames_ndpcf <- vector("list", 5)
n_mppf <- vector("list", 5)

#names_ndpcm <- paste0("ndpc_m_", 0:4)
pcats <- 0:4
#cnames_ndpcm <- paste0("c_m_ndpc_", 0:4)
#n_mppm <- paste0("mppm_", 0:4)


for (i in 1:5){
  names_ndpcf[[i]] <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCat" = pcats[i], "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
  names_ndpcf[[i]]$sexCode <- as.factor(names_ndpcf[[i]]$sexCode)
  c_f_ndpc2 <- predict(final, names_ndpcf[[i]], interval = "confidence")
  cnames_ndpcf[[i]] <- c_f_ndpc2
  n_mpp_x <- data.frame("weight"=names_ndpcf[[i]]$weight, "mpp"=c_f_ndpc2[,1], "mppmin"=c_f_ndpc2[,2], "mppmax"=c_f_ndpc2[,3])
  n_mppf[[i]] <- n_mpp_x
}

#install.packages("RColorBrewer")
library(RColorBrewer)
display.brewer.pal(n = 5, name = 'YlOrRd')

lab_pc <- c("Not infected","Liver category = 1","Liver category = 2","Liver category = 3","Liver category = 4")
col_pal1 <- c(brewer.pal(n = 4, name = "YlOrRd"), "grey", "black")
#col_pal1 <- c("darkred", "green", rep("darkred",4))
col_pal2 <- c(brewer.pal(n = 4, name = "YlOrRd"),"grey")

at4allf <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode == "F", ], aes(x = weight, y = gsi, color = "Observations"), size = 2, shape = 19, alpha=0.5) +
  
  geom_ribbon(data = mppf_0, aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[1]), alpha = 0.2) +
  geom_line(data = mppf_0, aes(x = weight, y = mpp, color = lab_pc[1]), lwd = 1.2, alpha = 0.4)+
  geom_ribbon(data = mppf_1, aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[2]), alpha = 0.2) +
  geom_line(data = mppf_1, aes(x = weight, y = mpp, color = lab_pc[2]), lwd = 1.2, alpha = 0.5)+
  geom_ribbon(data = mppf_2, aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[3]), alpha = 0.3) +
  geom_line(data = mppf_2, aes(x = weight, y = mpp, color = lab_pc[3]), lwd = 1.2, alpha = 0.6)+
  geom_ribbon(data = mppf_3, aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[4]), alpha = 0.3) +
  geom_line(data = mppf_3, aes(x = weight, y = mpp, color = lab_pc[4]), lwd = 1.2, alpha = 0.7)+
  geom_ribbon(data = mppf_4, aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[5]), alpha = 0.5) +
  geom_line(data = mppf_4, aes(x = weight, y = mpp, color = lab_pc[5]), lwd = 1.2, alpha = 0.8)+
  
  scale_color_manual(name = "Legend", values = col_pal1) +
  scale_fill_manual(name = "Legend", values = col_pal2, guide = "none") +
  labs(y = "GSI", x = "Weight (kg)", title="Females") +
  theme(legend.position = "bottom") +
  theme_classic()
at4allf
#ggsave(filename = "GSI vs weight for all liver cat ggplot females.png", plot = at4allf, bg="white", width = 5.3, height = 4, dpi = 500)

#version 2
v2at4allf <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode == "F", ], aes(x = weight, y = gsi, color = "Observations"), size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = n_mppf[[1]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[1]), alpha = 0.6) +
  geom_line(data = n_mppf[[1]], aes(x = weight, y = mpp, color = lab_pc[1]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = n_mppf[[2]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[2]), alpha = 0.6) +
  geom_line(data = n_mppf[[2]], aes(x = weight, y = mpp, color = lab_pc[2]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = n_mppf[[3]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[3]), alpha = 0.6) +
  geom_line(data = n_mppf[[3]], aes(x = weight, y = mpp, color = lab_pc[3]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = n_mppf[[4]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[4]), alpha = 0.5) +
  geom_line(data = n_mppf[[4]], aes(x = weight, y = mpp, color = lab_pc[4]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = n_mppf[[5]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[5]), alpha = 0.5) +
  geom_line(data = n_mppf[[5]], aes(x = weight, y = mpp, color = lab_pc[5]), lwd = 1.2, alpha = 0.9)+
  scale_color_manual(name = "Legend", values = col_pal1) +
  scale_fill_manual(name = "Legend", values = col_pal2, guide = "none") +
  labs(y = "GSI", x = "Weight (kg)", title="Females") +
  theme(legend.position = "bottom") +
  theme_classic()
v2at4allf
#ggsave(filename = "v2 GSI vs weight for all liver cat ggplot females.png", plot = v2at4allf, bg="white", width = 5.3, height = 4, dpi = 500)


# Extract the slopes

coef_parp <- data.frame("sexCode"=rep("F",5), "coef_gsiw"=rep(0,5), "liverCat"=0:4)

coef_parp$coef_gsiw[1] <- coef(lm(mpp ~ weight, data = n_mppf[[1]]))[2]
coef_parp$coef_gsiw[2] <- coef(lm(mpp ~ weight, data = n_mppf[[2]]))[2]
coef_parp$coef_gsiw[3] <- coef(lm(mpp ~ weight, data = n_mppf[[3]]))[2]
coef_parp$coef_gsiw[4] <- coef(lm(mpp ~ weight, data = n_mppf[[4]]))[2]
coef_parp$coef_gsiw[5] <- coef(lm(mpp ~ weight, data = n_mppf[[5]]))[2]

coef_parp

col_pal3 <- c("grey", brewer.pal(n = 4, name = "YlOrRd"))

coefs_lc <- ggplot(coef_parp, aes(x=liverCat, y=coef_gsiw, fill = factor(liverCat)))+
  geom_bar(stat="identity", color = "black", linewidth = 0.4)+
  scale_fill_manual(values = col_pal3) +
  labs(y = "Slope (GSI ~ weight (kg))", x = "Liver category", title="Females") +
theme_classic()+
guides(fill = "none")
#ggsave(filename = "all lc gsi vs w coefs females.png", plot = coefs_lc, bg="white", width = 4, height = 4, dpi = 500)


#   Repeat but with liver cat as a factor
summary(final)
finalF <- lm(formula = gsi ~ log(hsi) + weight + K + liverCatF + sexCode + 
               year + log(hsi):weight + log(hsi):sexCode + weight:K + weight:liverCatF + 
               weight:sexCode + K:sexCode + K:year + log(hsi):weight:sexCode,data = datc2)

Fnames_ndpcf <- vector("list", 5)
Fcnames_ndpcf <- vector("list", 5)
Fn_mppf <- vector("list", 5)
Fpcats <- 0:4

for (i in 1:5){
  Fnames_ndpcf[[i]] <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "liverCatF" = Fpcats[i], "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
  Fnames_ndpcf[[i]]$sexCode <- as.factor(Fnames_ndpcf[[i]]$sexCode)
  Fnames_ndpcf[[i]]$liverCatF <- as.factor(Fnames_ndpcf[[i]]$liverCatF)
  Fc_f_ndpc2 <- predict(finalF, Fnames_ndpcf[[i]], interval = "confidence")
  Fcnames_ndpcf[[i]] <- Fc_f_ndpc2
  Fn_mpp_x <- data.frame("weight"=Fnames_ndpcf[[i]]$weight, "mpp"=Fc_f_ndpc2[,1], "mppmin"=Fc_f_ndpc2[,2], "mppmax"=Fc_f_ndpc2[,3])
  Fn_mppf[[i]] <- Fn_mpp_x
}
Fcoef_parp <- data.frame("sexCode"=rep("F",5), "coef_gsiw"=rep(0,5), "liverCatF"=0:4)

Fcoef_parp$coef_gsiw[1] <- coef(lm(mpp ~ weight, data = Fn_mppf[[1]]))[2]
Fcoef_parp$coef_gsiw[2] <- coef(lm(mpp ~ weight, data = Fn_mppf[[2]]))[2]
Fcoef_parp$coef_gsiw[3] <- coef(lm(mpp ~ weight, data = Fn_mppf[[3]]))[2]
Fcoef_parp$coef_gsiw[4] <- coef(lm(mpp ~ weight, data = Fn_mppf[[4]]))[2]
Fcoef_parp$coef_gsiw[5] <- coef(lm(mpp ~ weight, data = Fn_mppf[[5]]))[2]

Fcoef_parp

Fcol_pal3 <- c("grey", brewer.pal(n = 4, name = "YlOrRd"))

Fcoefs_lc <- ggplot(Fcoef_parp, aes(x=liverCatF, y=coef_gsiw, fill = factor(liverCatF)))+
  geom_bar(stat="identity", color = "black", linewidth = 0.4)+
  scale_fill_manual(values = Fcol_pal3) +
  labs(y = "Slope (GSI ~ weight (kg))", x = "Liver category (factor)", title="Females") +
  theme_classic()+
  guides(fill = "none")
Fcoefs_lc
#ggsave(filename = "all lc gsi vs w coefs females for livercat factor.png", plot = Fcoefs_lc, bg="white", width = 4, height = 4, dpi = 500)

################################################################################
############         PLOT ALLPARASITE CATEGORIES : male fish
################################################################################


pcats <- 0:4
names_ndpcm <- vector("list", 5)
cnames_ndpcm <- vector("list", 5)
n_mppm <- vector("list", 5)


for (i in 1:5){
  names_ndpcm[[i]] <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "liverCat" = pcats[i], "weight" = seq(min(datc2$weight[datc2$sexCode=="M"]),max(datc2$weight[datc2$sexCode=="M"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
  names_ndpcm[[i]]$sexCode <- as.factor(names_ndpcm[[i]]$sexCode)
  c_m_ndpc2 <- predict(final, names_ndpcm[[i]], interval = "confidence")
  cnames_ndpcm[[i]] <- c_m_ndpc2
  n_mpp_x <- data.frame("weight"=names_ndpcm[[i]]$weight, "mpp"=c_m_ndpc2[,1], "mppmin"=c_m_ndpc2[,2], "mppmax"=c_m_ndpc2[,3])
  n_mppm[[i]] <- n_mpp_x
}

lab_pc <- c("Not infected","Liver category = 1","Liver category = 2","Liver category = 3","Liver category = 4")
col_pal1 <- c(brewer.pal(n = 4, name = "YlOrRd"), "grey", "black")
#col_pal1 <- c("darkred", "green", rep("darkred",4))
col_pal2 <- c(brewer.pal(n = 4, name = "YlOrRd"),"grey")

v2at4allm <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode == "M", ], aes(x = weight, y = gsi, color = "Observations"), size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = n_mppm[[1]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[1]), alpha = 0.6) +
  geom_line(data = n_mppm[[1]], aes(x = weight, y = mpp, color = lab_pc[1]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = n_mppm[[2]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[2]), alpha = 0.6) +
  geom_line(data = n_mppm[[2]], aes(x = weight, y = mpp, color = lab_pc[2]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = n_mppm[[3]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[3]), alpha = 0.6) +
  geom_line(data = n_mppm[[3]], aes(x = weight, y = mpp, color = lab_pc[3]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = n_mppm[[4]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[4]), alpha = 0.5) +
  geom_line(data = n_mppm[[4]], aes(x = weight, y = mpp, color = lab_pc[4]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = n_mppm[[5]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = lab_pc[5]), alpha = 0.5) +
  geom_line(data = n_mppm[[5]], aes(x = weight, y = mpp, color = lab_pc[5]), lwd = 1.2, alpha = 0.9)+
  scale_color_manual(name = "Legend", values = col_pal1) +
  scale_fill_manual(name = "Legend", values = col_pal2, guide = "none") +
  labs(y = "GSI", x = "Weight (kg)", title="Males") +
  theme(legend.position = "bottom") +
  theme_classic()
v2at4allm
#ggsave(filename = "v2 GSI vs weight for all liver cat ggplot males.png", plot = v2at4allm, bg="white", width = 5.3, height = 4, dpi = 500)


################################################################################
############             ESTIMATE THE NUMBER OF WORMS              #############
################################################################################

#### Thist step is based on the model presented by Ryberg et al. 2020

# Need HSI and length in cm
# length is in mm
datc2$TL <- datc2$length/10

# Liver category
datc2$livercat <- datc2$liverCat

datc2$log_nworms <- NA

# Use parameters from Ryberg et al 2020
for (i in 1:nrow(datc2)){
  if (datc2$livercat[i]==0) {
    datc2$log_nworms[i] <- -0.043*datc2$hsi[i] + 0.049*datc2$TL[i] - 0.731
  } 
  if (datc2$livercat[i]==1) {
    datc2$log_nworms[i] <- -0.043*datc2$hsi[i] + 0.049*datc2$TL[i] + 0.787
  }
  if (datc2$livercat[i]==2) {
      datc2$log_nworms[i] <- -0.043*datc2$hsi[i] + 0.049*datc2$TL[i] + 1.625
  }
  if (datc2$livercat[i]==3) {
    datc2$log_nworms[i] <- -0.043*datc2$hsi[i] + 0.049*datc2$TL[i] + 1.971
  }
  if (datc2$livercat[i]==4) {
    datc2$log_nworms[i] <- -0.043*datc2$hsi[i] + 0.049*datc2$TL[i] + 2.395
  }
}

datc2$nworms <- exp(datc2$log_nworms)
summary(datc2$nworms)

# Claculate worms per gram liver - infection density
datc2$wormPgl <- datc2$nworms/datc2$weightLiver

# Create boxplots of number of worms
boxplot(datc2$nworms~datc2$livercat)

my_colors <- c("darkgreen", "lightgreen", "#f1a340", "red", "red3")
lw1 <- ggplot(datc2, aes(x = livercat, y = nworms, group=livercat, fill = as.factor(livercat))) +
  geom_boxplot() +
  scale_fill_manual(values = my_colors) +
  labs(x = "Liver category", y = "Predicted number of worms") +
  theme_classic()+
  theme(legend.position = "none",axis.text = element_text(size = 12),axis.title = element_text(size = 14))
#ggsave(filename = "livworm livercat from ryberg model.png", plot = lw1, bg="white", width = 5.3, height = 4, dpi = 500)

#Looking for obvious relationship between worm, and wormpgl on gsi
plot(datc2$gsi, datc2$nworms)
hist(datc2$nworms) # count data follows Poisson
#try log
hist(log(datc2$nworms))
plot(datc2$gsi, log(datc2$nworms))


# Infection density
hist(datc2$wormPgl)
hist(log(datc2$wormPgl))
#more normally distributed
plot(datc2$gsi, log(datc2$wormPgl))


#### Invetsigate relationship between variables

ggplot(datc2, aes(x = wormPgl, y = K, color = sexCode)) +
  geom_point() +
  geom_smooth(aes(group = sexCode), method = lm, formula = y ~ log(x))+
  labs(x = expression("Infection density (worm g liver"^-1~")"), y = "Fulton's K") +
  scale_color_manual(values = c("M" = "darkgreen", "F" = "gold")) +  # Set colors
  labs(color = "Sex")+
  theme_minimal()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 14))
#ggsave(filename = "livworm vs k bits1.png", plot = p5, bg="white", width = 5.3, height = 4, dpi = 500)

lw2 <- ggplot(datc2, aes(x = log(wormPgl), y = K)) +
  geom_point() +
  geom_smooth(method = lm)+
  labs(x = expression("ln(Infection density (worm g liver"^-1~"))"), y = "Fulton's K") +
  theme_classic()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 14))
lw2
#ggsave(filename = "nworm vs k classic.png", plot = lw2, bg="white", width = 5.3, height = 4, dpi = 500)
kmlw2 <- lm(weight~log(wormPgl), datc2)
summary(kmlw2) 


lw3 <- ggplot(datc2, aes(x = wormPgl, y = K)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ log(x))+
  labs(x = expression("Infection density (worm g liver"^-1~")"), y = "Fulton's K") +
  theme_minimal()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 14))
#ggsave(filename = "nworm vs k log curve.png", plot = lw3, bg="white", width = 5.3, height = 4, dpi = 500)

lw4 <- ggplot(datc2, aes(x = wormPgl, y = gsi, color=sexCode)) +
  geom_point() +
  geom_smooth(method = lm)+
  labs(x = expression("Infection density (worm g liver"^-1~")"), y = "GSI") +
  scale_color_manual(values = c("M" = "lightblue3", "F" = "purple"))+
  labs(color = "Sex")+
  theme_minimal()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 14),legend.text = element_text(size = 12))
#ggsave(filename = "gsi vs infden.png", plot = lw4, bg="white", width = 5.3, height = 4, dpi = 500)

lw4test <- ggplot(datc2, aes(x = wormPgl, y = gsi)) +
  geom_point() +
  geom_smooth(method = lm)+
  labs(x = expression("Infection density (worm g liver"^-1~")"), y = "GSI") +
  theme_minimal()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 14))

gsilw1 <- lm(gsi~log(wormPgl), datc2)
summary(gsilw1) #best fits data
gsilw2 <- lm(gsi~wormPgl, datc2)
summary(gsilw2)
gsilw3 <- lm(gsi~wormPgl+sexCode, datc2)
summary(gsilw3)
gsilw4 <- lm(gsi~wormPgl*sexCode, datc2)
summary(gsilw4)
gsilw5 <- lm(gsi~log(wormPgl)*sexCode, datc2)
summary(gsilw5)
gsilw6 <- lm(gsi~log(wormPgl)+sexCode, datc2)
summary(gsilw6)
AICtab(gsilw1,gsilw2,gsilw3,gsilw4,gsilw5,gsilw6)
BICtab(gsilw1,gsilw2,gsilw3,gsilw4,gsilw5,gsilw6)
par(mfrow=c(2,2))
plot(gsilw4)

# Against weight
wlw1 <- lm(weight~log(wormPgl), datc2)
summary(wlw1)
wlw2 <- lm(weight~wormPgl, datc2)
summary(wlw2)
wlw3 <- lm(weight~wormPgl+sexCode, datc2)
summary(wlw3)
wlw4 <- lm(weight~wormPgl*sexCode, datc2)
summary(wlw4)
wlw5 <- lm(weight~log(wormPgl)*sexCode, datc2)
summary(wlw5)
wlw6 <- lm(weight~log(wormPgl)+sexCode, datc2)
summary(wlw6) ##best fits data
AICtab(wlw1,wlw2,wlw3,wlw4,wlw5,wlw6)
#BICtab(wlw1,wlw2,wlw3,wlw4,wlw5,wlw6) #lets rather use AIC
par(mfrow=c(2,2))
plot(wlw5)


lw5 <- ggplot(datc2, aes(x = wormPgl, y = weight, color=sexCode)) +
  geom_point() +
  geom_smooth(aes(group = sexCode), method = lm, formula = y ~ log(x))+
  labs(x = expression("Infection density (worm g liver"^-1~")"), y = "Weight (kg)") +
  scale_color_manual(values = c("M" = "darkgreen", "F" = "gold"))+
  labs(color = "Sex")+
  theme_minimal()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 14),legend.text = element_text(size = 12))
#ggsave(filename = "weight vs infden.png", plot = lw5, bg="white", width = 5.3, height = 4, dpi = 500)

lw6 <- ggplot(datc2, aes(x = log(wormPgl), y = weight, color=sexCode)) +
  geom_point() +
  geom_smooth(aes(group = sexCode), method = lm, formula = y ~ x)+
  labs(x = expression("ln(Infection density (worm g liver"^-1~"))"), y = "Weight (kg)") +
  scale_color_manual(values = c("M" = "darkgreen", "F" = "gold"))+
  labs(color = "Sex")+
  theme_classic()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key = element_blank())
lw6
#ggsave(filename = "weight vs infden new both sexes.png", plot = lw6, bg="white", width = 5.3, height = 4, dpi = 500)

lw7 <- ggplot(datc2, aes(x = liverCat, y = weight, color=sexCode)) +
  geom_point() +
  geom_smooth(aes(group = sexCode), method = lm, formula = y ~ x)+
  labs(x = "Liver Category", y = "Weight (kg)") +
  scale_color_manual(values = c("M" = "darkgreen", "F" = "gold"))+
  labs(color = "Sex")+
  theme_classic()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key = element_blank())
lw7

# Against weight
hsiifd <- lm(hsi~wormPgl, datc2)
summary(hsiifd)
hsiifd1 <- lm(hsi~wormPgl*sexCode, datc2)
summary(hsiifd1)
hsiifd2 <- lm(hsi~wormPgl+sexCode, datc2)
summary(hsiifd2)
hsiifd3 <- lm(hsi~log(wormPgl)+sexCode, datc2)
summary(hsiifd3)
hsiifd4 <- lm(hsi~log(wormPgl)*sexCode, datc2)
summary(hsiifd4)
hsiifd5 <- lm(hsi~log(wormPgl), datc2)
summary(hsiifd5)
AICtab(hsiifd, hsiifd1, hsiifd2, hsiifd3, hsiifd4, hsiifd5)
#BICtab(wlw1,wlw2,wlw3,wlw4,wlw5,wlw6) #lets rather use AIC
par(mfrow=c(2,2))
plot(hsiifd4)

ifd1 <- ggplot(datc2, aes(x = wormPgl, y = hsi, color=sexCode)) +
  geom_point() +
  geom_smooth(aes(group = sexCode), method = lm, formula = y ~ x)+
  labs(x = expression("Infection density (worm g liver"^-1~")"), y = "HSI") +
  scale_color_manual(values = c("M" = "blue", "F" = "red"))+
  labs(color = "Sex")+
  theme_minimal()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 14),legend.text = element_text(size = 12))
#ggsave(filename = "hsi vs ifd.png", plot = ifd1, bg="white", width = 5.3, height = 4, dpi = 500)

hsivsk <- ggplot(datc2, aes(x = hsi, y = K)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ x)+
  labs(x = "HSI", y = "Fulton's K") +
  theme_minimal()
#ggsave(filename = "hsi vs k.png", plot = hsivsk, bg="white", width = 5.3, height = 4, dpi = 500)

klw1 <- lm(K~log(wormPgl), datc2)
summary(klw1) # highly significant but only 0.15 R-squared
klw2 <- lm(K~wormPgl, datc2)
summary(klw2) # only 0.11 R-squared
AICtab(klw1,klw2) #log miuch better


################################################################################
#######          Investigate effects of infection density on GSI        ########
################################################################################

# Start from complete model
mwpgl <- step(lm(gsi ~ (hsi+weight+log(wormPgl)+sexCode+year+K)^3, datc2), trace = FALSE)
drop1(mwpgl, test = "F")
mwpgl <- update(mwpgl, ~ . -log(wormPgl):sexCode)
drop1(mwpgl, test = "F")
mwpgl <- update(mwpgl, ~ . -hsi:sexCode:year)
drop1(mwpgl, test = "F")
mwpgl <- update(mwpgl, ~ . -hsi:weight:year)
drop1(mwpgl, test = "F")
mwpgl <- update(mwpgl, ~ . -hsi:year)
drop1(mwpgl, test = "F")
mwpgl <- update(mwpgl, ~ . -weight:sexCode:year)
drop1(mwpgl, test = "F")
mwpgl <- update(mwpgl, ~ . -weight:year)
drop1(mwpgl, test = "F")
mwpgl <- update(mwpgl, ~ . -sexCode:year)
drop1(mwpgl, test = "F")
mwpgl <- update(mwpgl, ~ . -weight:sexCode:K)
drop1(mwpgl, test = "F")
mwpgl <- update(mwpgl, ~ . -weight:K)
drop1(mwpgl, test = "F")

mwpglv4 <- step(lm(gsi ~ (log(hsi)+weight+wormPgl+sexCode+year+K)^3, datc2), trace = FALSE)
drop1(mwpglv4, test = "F")
mwpglv4 <- update(mwpglv4, ~ . -year:K)
drop1(mwpglv4, test = "F")
mwpglv4 <- update(mwpglv4, ~ . -wormPgl:K)
drop1(mwpglv4, test = "F")
mwpglv4 <- update(mwpglv4, ~ . -weight:sexCode:K)
drop1(mwpglv4, test = "F")
mwpglv4 <- update(mwpglv4, ~ . -weight:K)
drop1(mwpglv4, test = "F")

mwpglv2 <- lm(gsi ~ hsi + weight + wormPgl + sexCode + year + 
                K + hsi:weight + hsi:wormPgl + hsi:sexCode + weight:wormPgl + 
                weight:sexCode + sexCode:K + year:K + hsi:weight:sexCode, 
              data = datc2)
mwpglv3 <- lm(gsi ~ log(hsi) + weight + wormPgl + sexCode + year + 
                K + log(hsi):weight + log(hsi):wormPgl + log(hsi):sexCode + weight:wormPgl + 
                weight:sexCode + sexCode:K + year:K + log(hsi):weight:sexCode, 
              data = datc2)
AICtab(mwpgl,mwpglv2,mwpglv3,mwpglv4) # linear better

drop1(mwpglv2, test = "F")
mwpglv2 <- update(mwpglv2, ~ . -year:K)
drop1(mwpglv2, test = "F")

summary(mwpglv4)
par(mfrow=c(2,2))
plot(mwpglv4, col=datc2$sexCode) #500 x 500
legend("bottomright", legend = unique(datc2$sexCode), fill = unique(datc2$sexCode), cex=0.5)
par(mfrow=c(1,1))
hist(resid(mwpglv4)) #300 x 250

## Without taking the log of infection density
mwpgl3 <- step(lm(gsi ~ (wormPgl+weight+sexCode+year+K)^3, datc2), trace = FALSE)
drop1(mwpgl3, test = "F")
mwpgl3 <- update(mwpgl3, ~ . -wormPgl:year:K)
drop1(mwpgl3, test = "F")
mwpgl3 <- update(mwpgl3, ~ . -weight:sexCode:K)
drop1(mwpgl3, test = "F")
mwpgl3 <- update(mwpgl3, ~ . -weight:K)
drop1(mwpgl3, test = "F")
mwpgl3 <- update(mwpgl3, ~ . -wormPgl:weight:year)
drop1(mwpgl3, test = "F")
mwpgl3 <- update(mwpgl3, ~ . -wormPgl:year)
drop1(mwpgl3, test = "F")
mwpgl3 <- update(mwpgl3, ~ . -year:K)
drop1(mwpgl3, test = "F")
mwpgl3 <- update(mwpgl3, ~ . -weight:sexCode:year)
drop1(mwpgl3, test = "F")
mwpgl3 <- update(mwpgl3, ~ . -weight:year)
drop1(mwpgl3, test = "F")
mwpgl3 <- update(mwpgl3, ~ . -sexCode:year)
drop1(mwpgl3, test = "F")
summary(mwpgl3)

AICtab(mwpgl,mwpgl2,mwpgl3)
BICtab(mwpgl,mwpgl2,mwpgl3) #most complex model is best


################################################################################
############                 PREDICT FROM  mwpglv4                 #############
################################################################################

# How is infection density distributed
par(mfrow=c(1,2))
hist(datc2$wormPgl[datc2$sexCode=="F"], main="Females", xlab="Infection density")
hist(datc2$wormPgl[datc2$sexCode=="M"], main="Males", xlab="Infection density")
# 640 vs 300


# Plot - gsi against weight for infection density and sexes
summary(datc2$wormPgl)

new_flowp4 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "wormPgl" = c(0.07254), "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
new_fhighp4 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "wormPgl" = c(0.50268), "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"])) #0.50268 3rd q
new_mlowp4 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "wormPgl" = c(0.13805), "weight" = seq(min(datc2$weight[datc2$sexCode=="M"]),max(datc2$weight[datc2$sexCode=="M"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
new_mhighp4 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "wormPgl" = c(0.90152), "weight" = seq(min(datc2$weight[datc2$sexCode=="M"]),max(datc2$weight[datc2$sexCode=="M"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"])) # 0.90152

new_flowp4$sexCode <- as.factor(new_flowp4$sexCode)
new_fhighp4$sexCode <- as.factor(new_fhighp4$sexCode)
new_mlowp4$sexCode <- as.factor(new_mlowp4$sexCode)
new_mhighp4$sexCode <- as.factor(new_mhighp4$sexCode)

c_flowp4 <- predict(mwpglv4, new_flowp4, interval = "confidence")
c_fhighp4 <- predict(mwpglv4, new_fhighp4, interval = "confidence")
c_mlowp4 <- predict(mwpglv4, new_mlowp4, interval = "confidence")
c_mhighp4 <- predict(mwpglv4, new_mhighp4, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$weight[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "Weight (kg)", ylab = "GSI")
matlines(new_flowp4$weight, c_flowp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "orange2")
matlines(new_fhighp4$weight, c_fhighp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "grey")

plot(datc2$weight[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = "Weight (kg)", ylab = "")
matlines(new_mlowp4$weight, c_mlowp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "orange2")
matlines(new_mhighp4$weight, c_mhighp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "grey")
legend(x = 0.9, y = 30, legend = c("Infection density = 1st Q.","Infection density = 3rd Q."), col = c("orange2","grey"), lty = 1, lwd=2, cex=0.6, bty="n")

# 800 x 500
mppif <- data.frame("weight"=new_flowp4$weight, "mpp"=c_flowp4[,1], "mppmin"=c_flowp4[,2], "mppmax"=c_flowp4[,3])
mppif2 <- data.frame("weight"=new_fhighp4$weight, "mpp"=c_fhighp4[,1], "mppmin"=c_fhighp4[,2], "mppmax"=c_fhighp4[,3])

if1 <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode=="F",], aes(x = weight, y = gsi, color = "Observations"),
             size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = mppif, aes(x = weight, ymin = mppmin, ymax = mppmax, fill = "Infection density = 1st Q."), alpha=0.5) + 
  geom_line(data = mppif, aes(x = weight, y= mpp, color = "Infection density = 1st Q."), lwd=1.2) +
  geom_ribbon(data = mppif2, aes(x = weight, ymin = mppmin, ymax = mppmax, fill = "Infection density = 3rd Q."), alpha=0.6) +
  geom_line(data = mppif2, aes(x = weight, y= mpp, color = "Infection density = 3rd Q."), lwd=1.2) +
  labs(y = "GSI", x = "Weight (kg)", title="Females") +
  scale_color_manual(name = "Legend", values = c("Observations" = "black", "Infection density = 1st Q." = "orange", "Infection density = 3rd Q." = "grey60")) +
  scale_fill_manual(name = "Legend", values = c("Infection density = 1st Q." = "orange", "Infection density = 3rd Q." = "grey60"), guide = "none") +
  theme(legend.position = "bottom") +
  theme_classic()
#ggsave(filename = "GSI vs liverCat for weights ggplot females legend.png", plot = at1fl, bg="white", width = 5.3, height = 4, dpi = 500)


mppm1 <- data.frame("liverCat"=new_mlowp5$liverCat, "mpp"=c_mlowp5[,1], "mppmin"=c_mlowp5[,2], "mppmax"=c_mlowp5[,3])
mppm12 <- data.frame("liverCat"=new_mhighp5$liverCat, "mpp"=c_mhighp5[,1], "mppmin"=c_mhighp5[,2], "mppmax"=c_mhighp5[,3])
at1m <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode=="M",], aes(x = liverCat, y = gsi, color = "Observations"),
             size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = mppm1, aes(x =liverCat, ymin = mppmin, ymax = mppmax, fill = "Weight = 0.3 kg"), alpha=0.5) +
  geom_line(data = mppm1, aes(x = liverCat, y= mpp, color = "Weight = 0.3 kg"), lwd=1.2) +
  geom_ribbon(data = mppm12, aes(x =liverCat, ymin = mppmin, ymax = mppmax, fill = "Weight = 1.25 kg"), alpha=0.6) +
  geom_line(data = mppm12, aes(x = liverCat, y= mpp, color = "Weight = 1.25 kg"), lwd=1.2) +
  labs(y = "GSI", x = "Parasite Code", title="Males") +
  scale_color_manual(name = "Legend", values = c("Observations" = "black", "Weight = 0.3 kg" = "lightgreen", "Weight = 1.25 kg" = "darkred")) +
  scale_fill_manual(name = "Legend", values = c("Weight = 0.3 kg" = "lightgreen", "Weight = 1.25 kg" = "darkred"), guide = "none") +
  theme(legend.position = "bottom") +
  theme_classic()
#ggsave(filename = "GSI vs liverCat for weights ggplot males.png", plot = at1m, bg="white", width = 5.3, height = 4, dpi = 500)



# Plot - GSI against infection density for sexes and for infection densities
new_flowp5 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "wormPgl" = seq(min(datc2$wormPgl[datc2$sexCode=="F"]),max(datc2$wormPgl[datc2$sexCode=="F"]), length=1000), "weight" = 0.3, "hsi" = mean(datc2$hsi[datc2$sexCode=="F"])) # summary(datc2$hsi[datc2$sexCode=="F"])
new_fhighp5 <- expand.grid("sexCode"= c("F"), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "wormPgl" = seq(min(datc2$wormPgl[datc2$sexCode=="F"]),max(datc2$wormPgl[datc2$sexCode=="F"]), length=1000), "weight" = 1, "hsi" = mean(datc2$hsi[datc2$sexCode=="F"])) # also tried close to max
new_mlowp5 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "wormPgl" = seq(min(datc2$wormPgl[datc2$sexCode=="M"]),max(datc2$wormPgl[datc2$sexCode=="M"]), length=1000), "weight" = 0.3, "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
new_mhighp5 <- expand.grid("sexCode"= c("M"), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "wormPgl" = seq(min(datc2$wormPgl[datc2$sexCode=="M"]),max(datc2$wormPgl[datc2$sexCode=="M"]), length=1000), "weight" = 1, "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))

new_flowp5$sexCode <- as.factor(new_flowp5$sexCode)
new_fhighp5$sexCode <- as.factor(new_fhighp5$sexCode)
new_mlowp5$sexCode <- as.factor(new_mlowp5$sexCode)
new_mhighp5$sexCode <- as.factor(new_mhighp5$sexCode)

c_flowp5 <- predict(mwpglv4, new_flowp5, interval = "confidence")
c_fhighp5 <- predict(mwpglv4, new_fhighp5, interval = "confidence")
c_mlowp5 <- predict(mwpglv4, new_mlowp5, interval = "confidence")
c_mhighp5 <- predict(mwpglv4, new_mhighp5, interval = "confidence")

par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$wormPgl[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = expression("Infection density (worm g liver"^-1~")"), ylab = "GSI")
matlines(new_flowp5$wormPgl, c_flowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightgreen")
matlines(new_fhighp5$wormPgl, c_fhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "darkred")
legend(x = 0.5, y = 33.5, legend = c("Weight = 0.3 kg.","Weight = 1 kg"), col = c("lightgreen","darkred"), lty = 1, lwd=2, cex=0.6, bty="n")

plot(datc2$wormPgl[datc2$sexCode=="M"], datc2$gsi[datc2$sexCode=="M"], pch=1, main="Males", xlab = expression("Infection density (worm g liver"^-1~")"), ylab = "")
matlines(new_mlowp5$wormPgl, c_mlowp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "lightgreen")
matlines(new_mhighp5$wormPgl, c_mhighp5[,1:3],
         lty = c(1,3,3), lw = 2, col = "darkred")
# 800 x 500



###### GGPLOT
par(mfrow=c(1,2), mar=c(5, 4, 2, 2))
plot(datc2$weight[datc2$sexCode=="F"], datc2$gsi[datc2$sexCode=="F"], pch=1, main="Females", xlab = "Weight (kg)", ylab = "GSI")
matlines(new_flowp4$weight, c_flowp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "orange2")
matlines(new_fhighp4$weight, c_fhighp4[,1:3],
         lty = c(1,3,3), lw = 2, col = "grey")

mpp <- data.frame("weight"=new_flowp4$weight, "mpp"=c_flowp4[,1], "mppmin"=c_flowp4[,2], "mppmax"=c_flowp4[,3])
mpp2 <- data.frame("weight"=new_fhighp4$weight, "mpp"=c_fhighp4[,1], "mppmin"=c_fhighp4[,2], "mppmax"=c_fhighp4[,3])

at <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode=="F",], aes(x = weight, y = gsi),
             size = 2, color="black", shape = 19, alpha=0.5) +
  geom_ribbon(data = mpp, aes(x =weight, ymin = mppmin, ymax = mppmax), alpha=0.5, fill="orange4") +
  geom_line(data = mpp, aes(x = weight, y= mpp), color="orange2") +
  geom_ribbon(data = mpp2, aes(x =weight, ymin = mppmin, ymax = mppmax), alpha=0.5, fill="lightgrey") +
  geom_line(data = mpp2, aes(x = weight, y= mpp), color="grey50") +
  labs(y = "GSI", x = "Weight (kg)") +
  theme_classic()


###     FINAL PLOTS: FEMALES

# Plot for all seven parasites 
infc_lab_exp <- expression("Infection density (worm g liver"^-1~")")
# split infection density
par(mfrow=c(1,2))
hist(datc2$wormPgl[datc2$sexCode=="F"], main="Females", xlab=infc_lab_exp)
hist(datc2$wormPgl[datc2$sexCode=="M"], main="Males", xlab=infc_lab_exp)
# 640 vs 400
sum_ifcd_f <- summary(datc2$wormPgl[datc2$sexCode=="F"])
sum_ifcd_f[1]

ifcd_pcats_f <- c(0,sum_ifcd_f[[2]], sum_ifcd_f[[3]], sum_ifcd_f[[5]], sum_ifcd_f[[6]])
ifcd_names_ndpcf <- vector("list", 5)
ifcd_cnames_ndpcf <- vector("list", 5)
ifcd_mppf <- vector("list", 5)


for (i in 1:5){
  ifcd_names_ndpcf[[i]] <- expand.grid("sexCode"= as.factor(c("F")), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "wormPgl" = ifcd_pcats_f[i], "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
  c_m_ndpc2 <- predict(mwpglv4, ifcd_names_ndpcf[[i]], interval = "confidence")
  ifcd_cnames_ndpcf[[i]] <- c_m_ndpc2
  n_mpp_x <- data.frame("weight"=ifcd_names_ndpcf[[i]]$weight, "mpp"=c_m_ndpc2[,1], "mppmin"=c_m_ndpc2[,2], "mppmax"=c_m_ndpc2[,3])
  ifcd_mppf[[i]] <- n_mpp_x
}


#install.packages("RColorBrewer")
library(RColorBrewer)
#display.brewer.pal(n = 5, name = 'YlOrRd')
fake_labs <- c("Not infected","1st Q","median","3rd Q","max")
#infd_lab <- c("Not infected",expression("Infection density (worm g liver"^-1~") 1"^st~"Q"), expression("Mean infection density (worm g liver"^-1~")"),expression("Infection density (worm g liver"^-1~") 1"^st~"Q"),expression("Maximum nfection density (worm g liver"^-1~")"))
#col_pal1 <- c(brewer.pal(n = 4, name = "YlOrRd"), "grey", "black")
#col_pal1 <- c("darkred", "green", rep("darkred",4))
#col_pal2 <- c(brewer.pal(n = 4, name = "YlOrRd"),"grey")

col_pal2 <- c("#FFFFB2","#FECC5C","#E31A1C","#FD8D3C","grey")
col_pal1 <- c("#FFFFB2","#FECC5C","#E31A1C","#FD8D3C","grey","black")

infd_lab <- c(expression("Infection density (worm g liver"^-1~") 1"^st~"Q") ,expression("Infection density (worm g liver"^-1~") 3"^rd~"Q"),expression("Maximum nfection density (worm g liver"^-1~")"), expression("Mean infection density (worm g liver"^-1~")"),"Not infected", "Observations")

infdf <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode == "F", ], aes(x = weight, y = gsi, color = "Observations"), size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = ifcd_mppf[[1]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[1]), alpha = 0.6) +
  geom_line(data = ifcd_mppf[[1]], aes(x = weight, y = mpp, color = fake_labs[1]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppf[[2]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[2]), alpha = 0.6) +
  geom_line(data = ifcd_mppf[[2]], aes(x = weight, y = mpp, color = fake_labs[2]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppf[[3]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[3]), alpha = 0.6) +
  geom_line(data = ifcd_mppf[[3]], aes(x = weight, y = mpp, color = fake_labs[3]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppf[[4]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[4]), alpha = 0.5) +
  geom_line(data = ifcd_mppf[[4]], aes(x = weight, y = mpp, color = fake_labs[4]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppf[[5]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[5]), alpha = 0.5) +
  geom_line(data = ifcd_mppf[[5]], aes(x = weight, y = mpp, color = fake_labs[5]), lwd = 1.2, alpha = 0.9)+
  scale_color_manual(name = "Legend", values=col_pal1, labels=infd_lab) +
  scale_fill_manual(name = "Legend", values=col_pal2, guide = "none") +
  labs(y = "GSI", x = "Weight (kg)", title="Females") +
  theme(legend.position = "bottom") +
  theme_classic()+
  theme(legend.text = element_text(hjust = 0))
infdf
#ggsave(filename = "v2 GSI vs weight for all liver cat ggplot males.png", plot = v2at4allm, bg="white", width = 5.3, height = 4, dpi = 500)




###   MALES 
sum_ifcd_m <- summary(datc2$wormPgl[datc2$sexCode=="M"])
#take mean wormPgl of all liver cats
pcats <- 0:4
mean_ifd_lv_m <- data.frame("livercat"=c(0:4), "sex"=rep("M", 5), "mean_ifd"=rep(0,5))
for (i in 1:5){
  mean_ifd_lv_m$mean_ifd[i] <- mean(datc2$wormPgl[datc2$sexCode=="M" & datc2$liverCat==pcats[i]])
}
mean_ifd_lv_m

#mean_ifd_lv_m2 <- aggregate(wormPgl~liverCat, subset(datc2, sexCode=="M"), mean)

ifcd_pcats_m <- c(mean_ifd_lv_m$mean_ifd[1],mean_ifd_lv_m$mean_ifd[2],mean_ifd_lv_m$mean_ifd[3],mean_ifd_lv_m$mean_ifd[4],mean_ifd_lv_m$mean_ifd[5])
ifcd_names_ndpcm <- vector("list", 5)
ifcd_cnames_ndpcm <- vector("list", 5)
ifcd_mppm <- vector("list", 5)


for (i in 1:5){
  ifcd_names_ndpcm[[i]] <- expand.grid("sexCode"= as.factor(c("M")), "K" = mean(datc2$K[datc2$sexCode=="M"]), "year"= c(2021), "wormPgl" = ifcd_pcats_m[i], "weight" = seq(min(datc2$weight[datc2$sexCode=="M"]),max(datc2$weight[datc2$sexCode=="M"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="M"]))
  c_m_ndpc2m <- predict(mwpglv4, ifcd_names_ndpcm[[i]], interval = "confidence")
  ifcd_cnames_ndpcm[[i]] <- c_m_ndpc2m
  n_mpp_xm <- data.frame("weight"=ifcd_names_ndpcm[[i]]$weight, "mpp"=c_m_ndpc2m[,1], "mppmin"=c_m_ndpc2m[,2], "mppmax"=c_m_ndpc2m[,3])
  ifcd_mppm[[i]] <- n_mpp_xm
}

col_pal2 <- c("#FFFFB2","#FECC5C","#FD8D3C","#E31A1C","grey")
col_pal1 <- c("#FFFFB2","#FECC5C","#FD8D3C","#E31A1C","grey","black")

infd_lab <- c(expression("Infection density (worm g liver"^-1~") 1"^st~"Q") ,expression("Infection density (worm g liver"^-1~") 3"^rd~"Q"),expression("Maximum nfection density (worm g liver"^-1~")"), expression("Mean infection density (worm g liver"^-1~")"),"Not infected", "Observations")

fake_labs <- c("Not infected","Liver category = 1","Liver category = 2","Liver category = 3","Liver category = 4")


infdm <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode == "M", ], aes(x = weight, y = gsi, color = "Observations"), size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = ifcd_mppm[[1]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[1]), alpha = 0.6) +
  geom_line(data = ifcd_mppm[[1]], aes(x = weight, y = mpp, color = fake_labs[1]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppm[[2]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[2]), alpha = 0.6) +
  geom_line(data = ifcd_mppm[[2]], aes(x = weight, y = mpp, color = fake_labs[2]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppm[[3]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[3]), alpha = 0.6) +
  geom_line(data = ifcd_mppm[[3]], aes(x = weight, y = mpp, color = fake_labs[3]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppm[[4]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[4]), alpha = 0.5) +
  geom_line(data = ifcd_mppm[[4]], aes(x = weight, y = mpp, color = fake_labs[4]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppm[[5]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[5]), alpha = 0.5) +
  geom_line(data = ifcd_mppm[[5]], aes(x = weight, y = mpp, color = fake_labs[5]), lwd = 1.2, alpha = 0.9)+
  scale_color_manual(name = expression("Infection density
      (worm g liver"^-1~")"), values=col_pal1) +
  scale_fill_manual(name = "Legend", values=col_pal2, guide = "none") +
  labs(y = "GSI", x = "Weight (kg)", title="Males") +
  theme(legend.position = "bottom") +
  theme_classic()+
  theme(legend.text = element_text(hjust = 0))
infdm
#ggsave(filename = "GSI vs weight for all liver cats but with ID ggplot males.png", plot = infdm, bg="white", width = 5.3, height = 4, dpi = 500)




### FEMALES - AGAIN
pcats <- 0:4
mean_ifd_lv_f <- data.frame("livercat"=c(0:4), "sex"=rep("F", 5), "mean_ifd"=rep(0,5))
for (i in 1:5){
  mean_ifd_lv_f$mean_ifd[i] <- mean(datc2$wormPgl[datc2$sexCode=="F" & datc2$liverCat==pcats[i]])
}
mean_ifd_lv_f

ifcd_pcats_f <- c(mean_ifd_lv_f$mean_ifd[1],mean_ifd_lv_f$mean_ifd[2],mean_ifd_lv_f$mean_ifd[3],mean_ifd_lv_f$mean_ifd[4],mean_ifd_lv_f$mean_ifd[5])
ifcd_names_ndpcf <- vector("list", 5)
ifcd_cnames_ndpcf <- vector("list", 5)
ifcd_mppf <- vector("list", 5)


for (i in 1:5){
  ifcd_names_ndpcf[[i]] <- expand.grid("sexCode"= as.factor(c("F")), "K" = mean(datc2$K[datc2$sexCode=="F"]), "year"= c(2021), "wormPgl" = ifcd_pcats_f[i], "weight" = seq(min(datc2$weight[datc2$sexCode=="F"]),max(datc2$weight[datc2$sexCode=="F"]), length=1000), "hsi" = mean(datc2$hsi[datc2$sexCode=="F"]))
  c_m_ndpc2f <- predict(mwpglv4, ifcd_names_ndpcf[[i]], interval = "confidence")
  ifcd_cnames_ndpcf[[i]] <- c_m_ndpc2f
  n_mpp_xf <- data.frame("weight"=ifcd_names_ndpcf[[i]]$weight, "mpp"=c_m_ndpc2f[,1], "mppmin"=c_m_ndpc2f[,2], "mppmax"=c_m_ndpc2f[,3])
  ifcd_mppf[[i]] <- n_mpp_xf
}

col_pal2 <- c("#FFFFB2","#FECC5C","#FD8D3C","#E31A1C","grey")
col_pal1 <- c("#FFFFB2","#FECC5C","#FD8D3C","#E31A1C","grey","black")

infd_lab <- c(expression("Infection density (worm g liver"^-1~") 1"^st~"Q") ,expression("Infection density (worm g liver"^-1~") 3"^rd~"Q"),expression("Maximum nfection density (worm g liver"^-1~")"), expression("Mean infection density (worm g liver"^-1~")"),"Not infected", "Observations")

fake_labs <- c("Not infected","Liver category = 1","Liver category = 2","Liver category = 3","Liver category = 4")

infdf2 <- ggplot(data = NULL) +
  geom_point(data = datc2[datc2$sexCode == "F", ], aes(x = weight, y = gsi, color = "Observations"), size = 2, shape = 19, alpha=0.5) +
  geom_ribbon(data = ifcd_mppf[[1]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[1]), alpha = 0.6) +
  geom_line(data = ifcd_mppf[[1]], aes(x = weight, y = mpp, color = fake_labs[1]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppf[[2]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[2]), alpha = 0.6) +
  geom_line(data = ifcd_mppf[[2]], aes(x = weight, y = mpp, color = fake_labs[2]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppf[[3]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[3]), alpha = 0.6) +
  geom_line(data = ifcd_mppf[[3]], aes(x = weight, y = mpp, color = fake_labs[3]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppf[[4]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[4]), alpha = 0.5) +
  geom_line(data = ifcd_mppf[[4]], aes(x = weight, y = mpp, color = fake_labs[4]), lwd = 1.2, alpha = 0.9)+
  geom_ribbon(data = ifcd_mppf[[5]], aes(x = weight, ymin = mppmin, ymax = mppmax, fill = fake_labs[5]), alpha = 0.5) +
  geom_line(data = ifcd_mppf[[5]], aes(x = weight, y = mpp, color = fake_labs[5]), lwd = 1.2, alpha = 0.9)+
  scale_color_manual(name = expression("Infection density
      (worm g liver"^-1~")"), values=col_pal1) +
  scale_fill_manual(name = expression("Infection density (worm g liver"^-1~")"), values=col_pal2, guide = "none") +
  labs(y = "GSI", x = "Weight (kg)", title="Females") +
  theme(legend.position = "bottom") +
  theme_classic()+
  theme(legend.text = element_text(hjust = 0))
infdf2
#ggsave(filename = "GSI vs weight for all liver cats but with ID ggplot females.png", plot = infdf2, bg="white", width = 5.3, height = 4, dpi = 500)
