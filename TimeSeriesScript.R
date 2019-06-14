###############################################################################
# Title:                                          Time Series Weather Analysis
# Version:                                                                    1
# Author:                                                     Kieran Billingham
# Date Created:                                                        22/04/19
# Document width:                                                 80 Characters
# Document Length: @£$%^&*&*()
# -----------------------------------------------------------------------------
# Notes:
# Ensure that all weather data is within a 'Data' subdirectory before executing
# script. 
# Access to CRAN is required for package installs. 
###############################################################################
# Clear Enviroment and Console ------------------------------------------------
rm(list=ls())
cat("\014")  

# Install And Load Packages ---------------------------------------------------
Packages <- c('data.table', 'tidyverse', 'readr', 'TSA', 'mgcv', 'forecast')

# Find packages that are not currently installed and download from CRAN
.NotInstalled <- Packages[!(Packages %in% installed.packages()[,'Package'])]
install.packages(.NotInstalled)

# load and attach add-on packages
lapply(Packages, require, character.only = TRUE)
rm(Packages)
# Set Paramaters --------------------------------------------------------------
.home.path <- '~/Documents/MSc/Analytics for Data Science/Assessment Two'
.data.path <- '/Data'
.model.path <- '/Model'
.image.path <- '/Images'

setwd(.home.path) 
dir.create(file.path(.home.path,.model.path), showWarnings = FALSE)
dir.create(file.path(.home.path,.image.path), showWarnings = FALSE)


.I.NA <- c('NA', 'na', '', '.', 'N/A', '#N/A', '..','Nan', '___', '---')
.k.NA <- -999999

# Read Data -------------------------------------------------------------------

#Get list of Files
filelist = list.files(path = paste(.home.path,.data.path, sep=''),
                      pattern = "*data.txt")

# Get Colum Headers
.cnames <- fread(paste(.home.path,.data.path,'/','aberporthdata.txt', sep=''),
                 sep='auto', skip = 'yyyy', na.strings = .I.NA, nrows=0)
ColHeaders <- colnames(.cnames)

#Import Data
data <- sapply(filelist, function(x) 
  fread(paste(.home.path,.data.path,'/',x, sep=''), 
        sep='auto', skip = 'sun', na.strings = .I.NA))

#Name Tables
names(data) <- (gsub('data.txt.*','', filelist))

#Witby is oddly formatted and so is removed to avoid future errors
data$whitby <- NULL

#Assign Col Headers
for (i in seq_along(data)){
  colnames(data[[i]]) <- ColHeaders
}

rm(filelist)
###############################################################################
# Clean and tidy --------------------------------------------------------------
# Create vector of first and last year data was recorded at each station
firstYear <- replicate(length(data), 0)
lastYear <- replicate(length(data), 0)

for (i in seq_along(data)) {
  firstYear[i] <- min(data[i][[1]][[1]])
  lastYear[i] <- max(data[i][[1]][[1]])
}

# Remove any that closed prior to 2018 (as per instructions)
logical <- lastYear != 2018
closedStations <- names(data)[logical]
data[closedStations] <- NULL
rm(logical)
# Find year all open stations were first active 
OpenYear <- max(firstYear)  #1978 ( Add one to ensure Month = Jan)
OpenYear <- OpenYear + 1

# Subset all data to match active years
for (i in seq_along(data)) {
  data[[i]] <- subset(data[[i]], data[[i]][[1]] >= OpenYear)
}

#Month Active - One station is not active until September. For complete years
#               Select 1979 as starter year (then have all starting at mm=1)
FirstMonth <- replicate(length(data), 0)

for (i in seq_along(data)) {
  FirstMonth[i] <- data[i][[1]][[2]][[1]]
}

# Remove markers for estimated and automatic readings 
for (i in seq_along(data)) {
  data[[i]] <- as.data.frame(sapply(data[[i]],function(x) 
    as.numeric(gsub('[#*]','',as.character(x)))))
}

###############################################################################
# Task One --------------------------------------------------------------------
MinTemps <- list()
MaxTemps <- list()
# Create Time Series Objects for Max and Min Temperatures
for (i in seq_along(data)) {
  MinTemps[[i]] <- ts(data[[i]][[4]], start = OpenYear, frequency = 12 )
  MaxTemps[[i]] <- ts(data[[i]][[3]], start = OpenYear, frequency = 12)
}

names(MinTemps) <- names(data)
names(MaxTemps) <- names(data)

# Exploritory Analysis of min and max temps
Range <- list()
Average <- list()

for (i in seq_along(MinTemps)) {
  Range[[i]] <- MaxTemps[[i]] - MinTemps[[i]]
  Average[[i]] <- as.data.frame(cbind(max(Range[[i]]),
                                      min(Range[[i]]),
                                      mean(Range[[i]]),
                                      max(MaxTemps[[i]]),
                                      mean(MaxTemps[[i]]),
                                      min(MinTemps[[i]]),
                                      mean(MinTemps[[i]])))
  
  colnames(Average[[i]]) <- c('MaxDifference',
                              'MinDifference',
                              'MeanDifference',
                              'Max',
                              'MeanMax',
                              'Min',
                              'MeanMin'
                              )
}

names(Range) <- names(data)
names(Average) <- names(data)


# UK minimum and max Temperature
max_T <- list()
min_T <- list()
for (i in seq_along(Average)) {
  max_T[i] <- Average[[i]][[4]]
  min_T[i] <- Average[[i]][[6]]
}
max_T <- max(unlist(max_T))
min_T <- min(unlist(min_T))

###############################################################################
# Task Two : Max Temps Trends -------------------------------------------------
# Max temperature trend estimates for each station
# PNG is created and saved for each in images directory
# Trend fits are saved to a list all.fits under station name
#
# Set subdirectory for any plots
dir.create(file.path(paste0(.home.path,.image.path,'/Trend')),
           showWarnings = FALSE)
all.max.fits <- list()

for (i in seq_along(MaxTemps)) {
  # Time points
  time.pts <- c(1:length(MaxTemps[[i]]))
  time.pts <- c(time.pts - min(time.pts))/max(time.pts)

  # Moving Average 
  mov_av.fit <- ksmooth(time.pts, MaxTemps[[i]], kernel='box')
  temp.mov_av.fit <- ts(mov_av.fit$y, 
                        start = OpenYear, frequency = 12)

  # Linear Regression
  linear.fit <- lm( MaxTemps[[i]] ~ time.pts)
  temp.linear.fit <- ts(fitted(linear.fit), start = OpenYear, frequency = 12)

  # Fit a quadratic and quartic polynomial 
  x2 <- time.pts^2
  x3 <- time.pts^3
  x4 <- time.pts^4
  poly2fit <- lm(MaxTemps[[i]] ~ time.pts + x2)
  poly4fit <- lm(MaxTemps[[i]] ~ time.pts + x2 +x3 +x4)
  temp.poly2.fit <- ts(fitted(poly2fit), start = OpenYear, frequency = 12)
  temp.poly4.fit <- ts(fitted(poly4fit), start = OpenYear, frequency = 12)

  # Local polynomial model
  local.fit <- loess(MaxTemps[[i]] ~ time.pts)
  temp.local.fit <- ts(fitted(local.fit), start = OpenYear, frequency = 12)

  # Splines
  gam.fit <- gam(MaxTemps[[i]] ~ s(time.pts))
  temp.gam.fit <- ts(fitted(gam.fit), start = OpenYear, frequency = 12)

  # Collate all trends
  all.temp.fit <- list(temp.mov_av.fit,
               temp.linear.fit, 
               temp.poly2.fit, 
               temp.poly4.fit, 
               temp.local.fit,
               temp.gam.fit)
  names(all.temp.fit) <- c("MAV","LINEAR","QUAD","QUART","LOCAL", "SPLINES")
  
  all.max.fits[[i]] <- all.temp.fit
  
  #Open file for plot 
  png(filename= paste0(.home.path,.image.path,'/Trend','/',
                       names(MaxTemps[i]),"_MaxTemp_Trends.png"))
  #set Ylim
  ylim <- c(min(unlist(all.temp.fit)), max(unlist(all.temp.fit)))
  # Plot trends
  ts.plot(temp.mov_av.fit, col = "coral",ylim = ylim, ylab = "Temperature")
  abline(temp.mov_av.fit[1], 0, lwd = 2, col='darkorange3')
  lines(temp.linear.fit, lwd = 2, col = "darkgreen")
  lines(temp.poly2.fit, lwd = 2, col = "cyan2")
  lines(temp.poly4.fit, lwd = 2, col = "deepskyblue4")
  lines(temp.local.fit, lwd = 2, col = "darkorchid")
  lines(temp.gam.fit, lwd = 2, col = "darkred")
  title(paste(names(MaxTemps[i]), 'maximum temperature trends'))
  legend('topleft',
         legend=c("MAV",'BASE',"LINEAR","QUAD","QUART","LOCAL", "SPLINES"),
         col=c("coral","darkorange3","darkgreen","cyan2","deepskyblue4",
               "darkorchid",'darkred'), lty=1, cex=0.8, inset=0.02)
  dev.off()
  
  #Clear enviroment for next loop
  rm(temp.mov_av.fit, temp.linear.fit, temp.poly2.fit, temp.poly4.fit, 
     temp.local.fit, temp.gam.fit, mov_av.fit, linear.fit, poly2fit,
     poly4fit, local.fit, gam.fit, time.pts, x2, x3, x4, ylim)
}

# Add names to the trend fits
names(all.max.fits) <- names(MaxTemps)


# Task Two : Min Temps Trends -------------------------------------------------
# Min temperature trend estimates for each station
# PNG is created and saved for each in images directory
# Trend fits are saved to a list all.fits under station name

all.min.fits <- list()

for (i in seq_along(MinTemps)) {
  # Time points
  time.pts <- c(1:length(MinTemps[[i]]))
  time.pts <- c(time.pts - min(time.pts))/max(time.pts)
  
  # Moving Average 
  mov_av.fit <- ksmooth(time.pts, MinTemps[[i]], kernel='box')
  temp.mov_av.fit <- ts(mov_av.fit$y, 
                        start = OpenYear, frequency = 12)
  
  # Linear Regression
  linear.fit <- lm( MinTemps[[i]] ~ time.pts)
  temp.linear.fit <- ts(fitted(linear.fit), start = OpenYear, frequency = 12)
  
  # Fit a quadratic and quartic polynomial 
  x2 <- time.pts^2
  x3 <- time.pts^3
  x4 <- time.pts^4
  poly2fit <- lm(MinTemps[[i]] ~ time.pts + x2)
  poly4fit <- lm(MinTemps[[i]] ~ time.pts + x2 +x3 +x4)
  temp.poly2.fit <- ts(fitted(poly2fit), start = OpenYear, frequency = 12)
  temp.poly4.fit <- ts(fitted(poly4fit), start = OpenYear, frequency = 12)
  
  # Local polynomial model
  local.fit <- loess(MinTemps[[i]] ~ time.pts)
  temp.local.fit <- ts(fitted(local.fit), start = OpenYear, frequency = 12)
  
  # Splines
  gam.fit <- gam(MinTemps[[i]] ~ s(time.pts))
  temp.gam.fit <- ts(fitted(gam.fit), start = OpenYear, frequency = 12)
  
  # Collate all trends
  all.temp.fit <- list(temp.mov_av.fit,
                       temp.linear.fit, 
                       temp.poly2.fit, 
                       temp.poly4.fit, 
                       temp.local.fit,
                       temp.gam.fit)
  names(all.temp.fit) <- c("MAV","LINEAR","QUAD","QUART","LOCAL", "SPLINES")

  all.min.fits[[i]] <- all.temp.fit
  
  #Open file for plot 
  png(filename= paste0(.home.path,.image.path,'/Trend','/',
                       names(MinTemps[i]),"_MinTemp_Trends.png"))
  #set Ylim
  ylim <- c(min(unlist(all.temp.fit)), max(unlist(all.temp.fit)))
  # Plot trends
  ts.plot(temp.mov_av.fit, col = "coral",ylim = ylim, ylab = "Temperature")
  abline(temp.mov_av.fit[1], 0, lwd = 2, col='darkorange3')
  lines(temp.linear.fit, lwd = 2, col = "darkgreen")
  lines(temp.poly2.fit, lwd = 2, col = "cyan2")
  lines(temp.poly4.fit, lwd = 2, col = "deepskyblue4")
  lines(temp.local.fit, lwd = 2, col = "darkorchid")
  lines(temp.gam.fit, lwd = 2, col = "darkred")
  title(paste(names(MinTemps[i]), 'minimum temperature trends'))
  legend('topleft',
         legend=c("MAV",'BASE',"LINEAR","QUAD","QUART","LOCAL", "SPLINES"),
         col=c("coral","darkorange3","darkgreen","cyan2","deepskyblue4",
               "darkorchid",'darkred'), lty=1, cex=0.8, inset=0.02)
  dev.off()
  
  #Clear enviroment for next loop
  rm(temp.mov_av.fit, temp.linear.fit, temp.poly2.fit, temp.poly4.fit, 
     temp.local.fit, temp.gam.fit, mov_av.fit, linear.fit, poly2fit,
     poly4fit, local.fit, gam.fit, time.pts, x2, x3, x4, ylim, all.temp.fit)
}

# Add names to the trend fits
names(all.min.fits) <- names(MinTemps)


# Compare fits between different stations 
plot(all.max.fits[[1]][[6]], col = 'red' , ylim = c(3.5,16),
     ylab = 'Temperature')
for (i in seq_along(all.max.fits)) {
  lines(all.max.fits[[i]][[6]], lwd = 2, col='red')
  lines(all.min.fits[[i]][[6]], lwd = 2, col='blue')
}
legend('topleft', legend=c("Max Trend", "Min Trend"),col=c("red","blue"),
       lty=1, cex=0.8, inset=0.02)
title('Maximum and Minimum Temperature Trends')


# Compare max and min difference 
plot(all.max.fits[[1]][[6]]-all.min.fits[[1]][[6]], col='green', ylim=c(4,9), ylab = 'Temperature Range')
for (i in seq_along(all.max.fits)) {
  r <-  all.max.fits[[i]][[6]]-all.min.fits[[i]][[6]]
  if (last(r) <= first(r)) 
  lines(r, lwd = 2, col = 'green') 
  else  lines(r, lwd = 2, col = 'red')
}
legend('topleft', legend=c("Decreasing", "Increasing"),col=c("Green","red"),
       lty=1, cex=0.8, inset=0.02)
title('trend ranges for each station')

# Task Two : Max Temp Seasonality ---------------------------------------------
# Set subdirectory for any plots
dir.create(file.path(paste0(.home.path,.image.path,'/Seasonality')),
           showWarnings = FALSE)

# Seasonality using ANOVA
 month <- season(MaxTemps[[1]])
 season.ANOVA <- lm(MaxTemps[[1]] ~ month - 1)

# Seasonality using Sin-Cos
 harm <- harmonic(MaxTemps[[1]], 2)
 season.harmonic <- lm(MaxTemps[[1]] ~ harm)

# Compare Methods
 season1 <- coef(season.ANOVA)
 season2 <- fitted(season.harmonic)[1:12]
#
 estimates <- data.frame(ANOVA = season1, HARMONIC = season2)
#
 plot(1:12, season1, lwd = 2, type = "l", col = "green",
      xlab = "Month", ylab = "Season Temperature")
 lines(1:12, season2, lwd = 2, col = "red")
 legend('topleft', legend=c("ANOVA", "harmoic"),col=c("green","red"),
        lty=1, cex=0.8, inset=0.02)
 title('Seasonality fits')

# Add Trend GAM trend found earlier
# Haromic seems to fit slightly better so continue with this.
# Create time points to replace dates
# time.pts <- c(1:length(MaxTemps[[1]])) # creates the time variable
# time.pts<-c((time.pts-min(time.pts))/max(time.pts))
# # Harmonic
# harm.gam.fit <- gam(MaxTemps[[1]] ~ s(time.pts) + harm)
#
# temp.harm.gam.fit <- ts(fitted(harm.gam.fit),
#                         start = OpenYear, frequency = 12)
#ANOVA
# anova.gam.fit <- gam(MaxTemps[[1]] ~ s(time.pts) + month)
#
# temp.anova.gam.fit <- ts(fitted(anova.gam.fit),
#                         start = OpenYear, frequency = 12)
# plot(MaxTemps[[1]])
# lines(temp.harm.gam.fit, col='blue')
# lines(temp.anova.gam.fit, col='red')

MaxTempTrends <- list()
MaxTempEstimates <- list()
SumSquaredResMaxT <- list()

for (i in seq_along(MaxTemps)) {
  
  harm <- harmonic(MaxTemps[[i]], 2)
  season.harmonic <- lm(MaxTemps[[i]] ~ harm)
  
  MaxTempEstimates[[i]] <- fitted(season.harmonic)[1:12]
  
  time.pts <- c(1:length(MaxTemps[[i]])) # creates the time variable
  time.pts<-c((time.pts-min(time.pts))/max(time.pts))
  
  harm.gam.fit <- gam(MaxTemps[[i]] ~ s(time.pts) + harm)
  temp.harm.gam.fit <- ts(fitted(harm.gam.fit),
                          start = OpenYear, frequency = 12)
  
  MaxTempTrends[[i]] <- temp.harm.gam.fit
  
  #Get SSR for each station
  SumSquaredResMaxT[[i]] <- sum((MaxTemps[[i]]-temp.harm.gam.fit)^2)
  
  #Open file for plot 
  png(filename= paste0(.home.path,.image.path,'/Seasonality','/',
                       names(MaxTemps[i]),"_MaxTemp_Seasonality.png"))
  #set Ylim
  ylim <- c(min(unlist(MaxTemps[[i]])), max(unlist(MaxTemps[[i]])))
  # Plot trends
  ts.plot(MaxTemps[[i]], col = "black",ylim = ylim, ylab = "Temperature")
  lines(temp.harm.gam.fit, lwd = 2, col = 'red')
  title(paste(names(MaxTemps[i]),'maximum temperature trend with seasonality'))
  legend('topleft',
         legend=c('Max Temp', 'Estimate'),
         col=c("black", "deepskyblue3"), lty=1, cex=0.8, inset=0.02)
  dev.off()
  
  #Clear enviroment for next loop
  rm(harm, time.pts, harm.gam.fit, temp.harm.gam.fit, ylim)
}

names(MaxTempTrends)<- names(MaxTemps)
names(MaxTempEstimates) <- names(MaxTemps)
names(SumSquaredResMaxT) <- names(MaxTemps)

# Task Two : Min Temp Seasonality ---------------------------------------------
# Repeate above processes for the minimum temperatures 
MinTempTrends <- list()
MinTempEstimates <- list()
SumSquaredResMinT <- list()

for (i in seq_along(MinTemps)) {
  
  harm <- harmonic(MinTemps[[i]], 2)
  season.harmonic <- lm(MinTemps[[i]] ~ harm)
  
  MinTempEstimates[[i]] <- fitted(season.harmonic)[1:12]
  
  time.pts <- c(1:length(MinTemps[[i]])) # creates the time variable
  time.pts<-c((time.pts-min(time.pts))/max(time.pts))
  
  harm.gam.fit <- gam(MinTemps[[i]] ~ s(time.pts) + harm)
  temp.harm.gam.fit <- ts(fitted(harm.gam.fit),
                          start = OpenYear, frequency = 12)
  
  MinTempTrends[[i]] <- temp.harm.gam.fit
  
  #Get SSR for each station
  SumSquaredResMinT[[i]] <- sum((MinTemps[[i]]-temp.harm.gam.fit)^2)
  
  #Open file for plot 
  png(filename= paste0(.home.path,.image.path,'/Seasonality','/',
                       names(MinTemps[i]),"_MinTemp_Seasonality.png"))
  #set Ylim
  ylim <- c(min(unlist(MinTemps[[i]])), max(unlist(MinTemps[[i]])))
  # Plot trends
  ts.plot(MinTemps[[i]], col = "black",ylim = ylim, ylab = "Temperature")
  lines(temp.harm.gam.fit, lwd = 2, col = 'deepskyblue3')
  title(paste(names(MinTemps[i]),'minimum temperature trend with seasonality'))
  legend('topleft',
         legend=c('Min Temp', 'Estimate'),
         col=c("black", "deepskyblue3"), lty=1, cex=0.8, inset=0.02)
  dev.off()
}

names(MinTempTrends)<- names(MinTemps)
names(MinTempEstimates) <- names(MinTemps)
names(SumSquaredResMinT) <- names(MinTemps)


###############################################################################
# Task 3: Model fitting and forecasting ---------------------------------------
# UK Max Temperature ----
UKAvMax <- ts(rowMeans(as.data.frame(MaxTemps)), 
              start = OpenYear, frequency = 12)

#Extract time series using difference function - only 1st level needed to get 
#stationary data. 1st order. 
diff.Max.ts <- diff(UKAvMax)

plot(diff.Max.ts, ylab = 'Differenced Temperature')
abline(a=0, b=0)

# Removing the trend and seasonality effectively differneces the data. 
# The seasonality can be seen in the acf and pacf plots and so will need
# removing
acf(diff.Max.ts)
pacf(diff.Max.ts) 

# Decompose - showing the trend and seasonality of the data as seen before
MaxComponents <- stl(UKAvMax, 'periodic')
plot(MaxComponents)

#Create linear model
Maxlinear <- tslm(UKAvMax ~ trend + season -1, data = UKAvMax)
linear.Max.fit <- ts(fitted(Maxlinear), start = OpenYear, frequency = 12)
# Plot linear predictions
plot(predict(linear.Max.fit, h=12,) ,ylab='Seasonality and Trend Temp', xlim= c(2017, 2020))

#remove linear trend from data 
MaxResidual <- UKAvMax - linear.Max.fit
# Check lag now
acf(MaxResidual)
pacf(MaxResidual) 

#fit model on residuals
Maxmodel <- auto.arima(MaxResidual, d=0, seasonal = FALSE)
tsdisplay(residuals(Maxmodel), lag.max=45, main='(1,1,1) Model Residuals')


#Get Forecasts for linear fit and ARMA Model
MaxpredRes <- forecast(Maxmodel, h=12)
plot(MaxpredRes, xlim = c(2017,2020), ylab ='Residual ARMA Prediction')
MaxNewRes <- MaxpredRes$mean
MaxpredTrend <- forecast(Maxlinear, h=12)
MaxNewTrend <- MaxpredTrend$mean

#Cobine to get predictions
NewMaxTemp <- MaxNewRes + MaxNewTrend

###
# UK Min Temperature ----
UKAvMin <- ts(rowMeans(as.data.frame(MinTemps)),
              start = OpenYear, frequency = 12)

#Extract time series using difference function - only 1st level needed to get
#stationary data. 1st order. 
diff.Min.ts <- diff(UKAvMin)

plot(diff.Min.ts, ylab = 'Differenced Temperature')
abline(a=0, b=0)

# Removing the trend and seasonality effectively differneces the data. 
# The seasonality can be seen in the acf and pacf plots and so will need 
# removing
acf(diff.Min.ts)
pacf(diff.Min.ts) 

# Decompose - showing the trend and seasonality of the data as seen before
MinComponents <- stl(UKAvMin, 'periodic')
plot(MinComponents)
Minlinear <- tslm(UKAvMin ~ trend + season -1, data = UKAvMin)

linear.Min.fit <- ts(fitted(Minlinear), start = OpenYear, frequency = 12)

#remove linear trend from data 
MinResidual <- UKAvMin - linear.Min.fit
#fit model on residuals
Minmodel <- auto.arima(MinResidual,d=0, seasonal = FALSE)
tsdisplay(residuals(Minmodel), lag.max=45, main='(1,1,1) Model Residuals')

#Get Forecasts for linear fit and ARMA Model
MinpredRes <- forecast(Minmodel, h=12)
MinNewRes <- MinpredRes$mean
MinpredTrend <- forecast(Minlinear, h=12)
MinNewTrend <- MinpredTrend$mean

#Combine to get predictions
NewMinTemp <- MinNewRes + MinNewTrend

#Plot Minimum and minimum temperature predicitons 
ylim <- c(min(c(NewMaxTemp,NewMinTemp)),max(c(NewMaxTemp,NewMinTemp)))

png(filename= paste0(.home.path,.image.path,'/',"2019Predictions.png"))

ts.plot(NewMaxTemp, col='red',ylim = ylim, 
        ylab = 'Predicted Temperature', xlab = 'Year.Month')
lines(NewMinTemp, col='blue')
legend('topleft',
       legend=c('Min Temp', 'Max Temp'),
       col=c("blue", "red"), lty=1, cex=0.8, inset=0.02)
grid()
title('Predicted UK Temperatures ')

dev.off()

100*5/24


###############################################################################
# Task 4: Visualising Results -------------------------------------------------

# install.packages('ggmap')
# require(ggmap)
# 
# #Enable google map services 
# register_google(key = 'XXXXXX')
# # Get UK Map 
# UKmap <- get_map('United Kingdom', zoom = 5)
# ggmap(UKmap) 

