# Install the necessary packages (run only once)
install.packages("readxl")
install.packages("quantmod")
install.packages("lmtest")
install.packages("urca")
install.packages("stargazer")
install.packages("dyn")

# Load the required libraries
library(readxl)
library(quantmod)
library(lmtest)
library(urca)
library(stargazer)
library(dyn)

##############################################################################
# 0.a. Preparing the data
##############################################################################
# load the data set from this link: https://www.princeton.edu/~mwatson/Stock-Watson_3u/Students/EE_Datasets/us_macro_quarterly.xlsx
# Set the path of your working directory (where the excel file is located)

setwd("C:\\Users\\Brend\\OneDrive\\Neoma Business School\\Msc Finance & Big Data\\Econometrics and time series\\Projet de Groupes")


# Import the data (do not take into account the warning message)

data <- as.data.frame(read_xlsx("us_macro_quarterly.xlsx",sheet = 1,col_types = c("text", rep("numeric", 9))))

# Get rid off unecessary variables

data <- data[,-which(colnames(data) %in% c("JAPAN_IP","PCECTPI","GS1","UNRATE","EXUSUK","CPIAUCSL"))]

# Change name of the GDP variable

names(data)[names(data) == "GDPC96"] <- "GDP"

# format the date column

names(data)[1]<- "Date"
data$Date <- as.yearqtr(data$Date, format = "%Y:0%q")

# Create GDP quarterly annualized growth (add NA to be the same length than nrow(data))

# Ussing the lag and diff functions (which require ts class)

GDP <- ts(data$GDP,start = c(1957,1),end = c(2013,4),frequency = 4)
data$GDPGrowth <- c(NA,400*(log(GDP)-lag(log(GDP),k=1)))# Equivalently c(NA,400*diff(log(GDP),lag=1))

# Create Term Spread variable (10Y-3m)

data$term_spread <- data$GS10 - data$TB3MS

# Compute the first difference of the Term Spread

data$delta_term_spread <- c(NA,data$term_spread[-1]-data$term_spread[-length(data$term_spread)])

# Create a dummy variable when the Term Spread turns negative

data$dummy_term_spread <- ifelse(data$term_spread < 0,1,0)

# Load NBER recessions data from FRED database

nber <- getSymbols("USRECQM", src="FRED",auto.assign = FALSE) # NBER US Recession
data$nber <- window(nber,start = "1957-01-01", end = "2013-10-01")

# Compute the two-quarter ahead NBER (National Bureau of Economic Recession) recessions

data$nber_lead2 <- c(lag(ts(data$nber,start = c(1957,1),end = c(2013,4),frequency = 4),k=-2)[-c(1:2)],rep(NA,2))

# Convert data to a ts object at quarterly frequency

data_ts <- ts(data, start = c(1957, 1), end = c(2013, 4), frequency = 4)

head(data_ts,5)# Eyeball the first five rows, or View(data_ts) to see the full dataset

# Remove Date column (no longer used as the data has been transformed to time series)

data_ts <- data_ts[,-which(colnames(data) %in% "Date")]
data_ts <- na.omit(data_ts) # get rid of the first row with NA

# Define time window for training and testing

start_training_date <- time(data_ts)[1] # first date
end_training_date <- time(data_ts)[dim(data_ts)[1]]-3 # final date minus three years

# Create the training dataset (all the data minus three years)
data_insample <- window(data_ts, start = start_training_date, end = end_training_date)

# Create the test dataset (2011, 2012 and 2013)
data_outofsample <- window(data_ts, start = end_training_date + 0.25, end = time(data_ts)[dim(data_ts)[1]])

##############################################################################
# 0.b. Plotting the data
##############################################################################
# Plot GDP Growth
plot.ts(data_insample[,"GDPGrowth"], ylab = "",xlab = "", main = "GDP growth",type="l", col = "steelblue")

# Plot the term spread
plot.ts(data_insample[,"term_spread"], ylab = "",xlab = "", main = "Term spread",type="l", col = "steelblue")

# Plot the term spread and recession events
plot(data_insample[,"term_spread"],col = "steelblue",lwd = 1, xlab = "",ylab = "Percent", main = "Term spread and recessions")
abline(h=0,col="red")
par(new = T)
plot(data_insample[,"nber"], axes=FALSE, xlab=NA, ylab=NA, cex=1.2,col="orange")

##############################################################################
# 0.c. Testing the stationarity
##############################################################################
# Augmented-Dickey-Fuller Unit Root Test
# H0: a unit root is present
summary(ur.df(data_insample[,"GDPGrowth"],type = "none", lags = 0))
summary(ur.df(data_insample[,"term_spread"],type = "drift", lags = 0))
# Results: |test statistic| > |critical values| at 5% => we reject the null that the unit root is present => the time series are stationary

##############################################################################
# 1. Univariate one-step ahead forecast (the benchmark): predicting GDP growth with its own lags
##############################################################################
# Determining the optimal lag number with the BIC
# BIC (Bayesian Information Criterion) is a useful statistic to assess the fit of a model by taking into account the number of parameters (remember that the fit of a model increases mechanically with the number of parameters)

# Create function that computes BIC for AR model
bic_ar <- function(fit){
  RSS <- sum(fit$residuals^2) # Residual Sum of Squares
  T <- length(fit$residuals)  # Number of observations
  K <- length(fit$coef)       # Number of parameters including the intercept
  return(c("Nb. of lags" = round(K-1,0),
           "BIC" = T*log(RSS/(T-K)) + K * log(T)))
}

# Apply BIC function over models of different  orders (from 1 to 8)
bic_ar_output <- sapply(1:8, function(x) "AR" = bic_ar(dynlm(GDPGrowth ~ L(GDPGrowth,1:x),data = data_insample)))

# select the AR model with the smallest BIC
p_best <- bic_ar_output[, which.min(bic_ar_output["BIC",])][1]

# Estimate the model with 4 lags
fit_ar<-dyn$lm(GDPGrowth ~ lag(GDPGrowth,-1) + lag(GDPGrowth,-2) + lag(GDPGrowth,-3) + lag(GDPGrowth,-4) ,data = data_insample)
summary(fit_ar)

predicted_values_ar_insample <- predict(fit_ar,data_insample)
predicted_values_ar_outofsample <-predict(fit_ar,data_outofsample)

# Plot predicted (in-sample one step-ahead prediction) vs actual values
plot(data_insample[,"GDPGrowth"], ylab = "", main = "Predicted vs actual GDP growth",type="l", col = "steelblue",cex.main =0.8)
lines(ts(predicted_values_ar_insample, start = start_training_date, end = end_training_date, frequency = 4),col="red")

# QUESTION 1:
# In-sample RMSE of the AR model

# Out-of-sample RMSE of the AR model



##############################################################################
# 2. One-step ahead forecast with the Augmented Distributed Lag model: predicting GDP growth with its own lags and the term spread
##############################################################################
# Estimating the models
# ADL(1,2) model
adl12 <-dyn$lm(GDPGrowth ~ lag(GDPGrowth,-1) + + lag(term_spread,-1) + lag(term_spread,-2) ,data = data_insample)
predict(adl12,data_insample)
predict(adl12,data_outofsample)

# ADL(2,1) model
adl21 <-dyn$lm(GDPGrowth ~ lag(GDPGrowth,-1) + lag(GDPGrowth,-2) + lag(term_spread,-1),data = data_insample)
predict(adl21,data_insample)
predict(adl21,data_outofsample)

# ADL(2,2) model
adl22 <-dyn$lm(GDPGrowth ~ lag(GDPGrowth,-1) + lag(GDPGrowth,-2) + lag(term_spread,-1) + + lag(term_spread,-2),data = data_insample)
predict(adl22,data_insample)
predict(adl22,data_outofsample)

# ADL(1,1) model with the first difference of the term spread
adl11_diff <-dyn$lm(GDPGrowth ~ lag(GDPGrowth,-1) + + diff(lag(term_spread,-1)),data = data_insample)
predict(adl11_diff,data_insample)
predict(adl11_diff,data_outofsample)

# QUESTION 3:
# In-sample RMSE of the four models


# Out-of-sample RMSE of the four models


#####################################
# 3. Predicting recessions (2 quarters ahead) with the term spread
#####################################
# Estimating the models
# With the level of the term spread
fit1 <- glm(nber_lead2 ~ nber + term_spread,family=binomial(link='logit'), data=data_insample)
# With the first difference of the term spread
fit2 <- glm(nber_lead2 ~ nber + delta_term_spread,family=binomial(link='logit'), data=data_insample)
# With the dummy term spread
fit3 <- glm(nber_lead2 ~ nber + dummy_term_spread,family=binomial(link='logit'), data=data_insample)

# Visualize regression outputs
stargazer(fit1,fit2,fit3, title="Predicting Recessions", type="text",
          dep.var.caption  = "Model",
          dep.var.labels   = "",
          column.labels = c("level", "first diff","dummy"))

# Generate the values predicted by the three models
pr_recession_fit1 <- predict(fit1,data_insample, type = "response")
pr_recession_fit2 <- predict(fit1,data_insample, type = "response")
pr_recession_fit3 <- predict(fit1,data_insample, type = "response")

# QUESTION 6: (in-sample) misclassification rate of a classifier that predicts only non-recessionary events


# QUESTION 7: threshold that maximises the (in-sample) accuracy of the three models


