###############################
### DATA ANALYSIS PROJECT 2 ###
### ----------------------- ###
### PM10 Pollution Analysis ###
###############################


library(lubridate) # for datetimes
library(plyr) # for mapvalues() to rename factors
library(glmnet) # for glmnet()
set.seed(1) # set random seed


################################
## Data Access & Manipulation ##
################################

# read in data set
LINK = "http://lib.stat.cmu.edu/datasets/PM10.dat"
FILE_NAME = "PM10.txt"
download.file(LINK, FILE_NAME)
pm10 = read.table(FILE_NAME, col.names=c('log_pm10_concentration', 'log_cars_per_hour',
                                         'temp_at_2meters', 'wind_speed',
                                         'temp_diff', 'wind_direction',
                                         'hour_of_day','day_num'))

# date after which observations were taken
START_DATE = as.Date('2001/10/01')
# determine date of each observation
pm10$date = as_datetime(sprintf('%s %s:00:00',
                                as.character(START_DATE + pm10$day_num),
                                as.character(pm10$hour_of_day)))

# create variable for season
# {1: spring, 2: summer, 3: fall, 4: winter}
pm10$season = ifelse(month(pm10$date) %in% c(3,4,5), 1,
                     ifelse(month(pm10$date) %in% c(6,7,8), 2,
                            ifelse(month(pm10$date) %in% c(9,10,11), 3, 4)))

# create variable for period of day
# {1: early morning, 2: mid morning, 3: mid day,
#  4: early afternoon, 5: late afternoon, 6: night}
pm10$period = ifelse(pm10$hour_of_day %in% 1:4, 1,
                     ifelse(pm10$hour_of_day %in% 5:8, 2,
                            ifelse(pm10$hour_of_day %in% 9:12, 3,
                                   ifelse(pm10$hour_of_day %in% 13:16, 4,
                                          ifelse(pm10$hour_of_day %in% 17:20, 5, 6)))))

# convert to dummy variables
pm10$season = as.factor(pm10$season)
pm10$period = as.factor(pm10$period)

# rename factors
pm10$season = mapvalues(pm10$season,
                        from = as.character(1:4),
                        to = c('spring','summer','fall','winter'))
pm10$period = mapvalues(pm10$period,
                        from = as.character(1:6),
                        to = c('early morning','mid morning','mid day',
                               'early afternoon','late afternoon','night'))

# no need for date column
pm10 = subset(pm10, select=-c(date))
# there is just one 'summer' observation so get rid of it
pm10 = pm10[which(pm10$season != 'summer'), ]
# normalize numeric variables
temp = subset(pm10, select=-c(season, period))
temp = mapply(temp, FUN=function(x) (x-min(x))/(max(x)-min(x)))
temp = as.data.frame(temp)
temp$season = pm10$season
temp$period = pm10$period
pm10 = temp

# function to generate diagnostic plots
get_diagnostic_plots = function(lm)
{
  dev.new(width=600, height=400)
  par(mfrow=c(2,3))
  plot(lm, which=c(1:6))
  par(mfrow=c(1,1))
}

get_avg_error = function(s)
{
  ###################
  # Partition Data ##
  ###################
  
  index = c(1:nrow(pm10))
  train_index = sample(index, size=s*length(index), replace=FALSE)
  test_index = index[-train_index]
  train = pm10[train_index, ]
  test = pm10[test_index, ]
  
  
  ####################################
  ## Linear Model w/out Interaction ##
  ####################################
  
  # initial model (excluding variable 'period')
  pm10.lm.1 = lm(log_pm10_concentration ~ .,
               data=subset(train, select=-c(period)))
  # backward stepwise regression to reduce variable set
  pm10.lm.1 = step(pm10.lm.1, direction='backward', trace=0)
  
  
  ###################################
  ## Linear Models w/ Interactions ##
  ###################################
  
  # initial model containing interaction terms
  pm10.lm.2 = lm(log_pm10_concentration ~ (.)^2,
               data=train)
  # backward stepwise regression to reduce variable set
  pm10.lm.2 = step(pm10.lm.2, direction='backward', trace=0)
  
  # manually pruned model
  pm10.lm.3 = lm(log_pm10_concentration ~ log_cars_per_hour + temp_at_2meters + 
                   wind_speed + season + log_cars_per_hour:wind_speed + 
                   temp_at_2meters:wind_speed, data=train)
  
  
  ###################
  ## One-Way ANOVA ##
  ###################
  
  pm10.aov = aov(log_pm10_concentration ~ ., data=train)
  pm10.aov = step(pm10.aov, direction='backward', trace=0)
  
  #################
  ## Lasso Model ##
  #################
  
  x.train = as.matrix(subset(train, select=-c(log_pm10_concentration, season, period)))
  x.test = as.matrix(subset(test, select=-c(log_pm10_concentration, season, period)))
  y.train = train$log_pm10_concentration
  y.test = test$log_pm10_concentration
  grid = 10 ^ seq(10, -2, length=100)
  
  # generate lasso model
  pm10.lasso = glmnet(x.train, y.train, alpha=1, lambda=grid)
  lasso.cv = cv.glmnet(x.train, y.train, alpha=1)
  lasso.best_lambda = lasso.cv$lambda.min
  lasso.estimates = predict(pm10.lasso, type='coefficients', s=lasso.best_lambda)
  cat('LASSO - Number of non-zero coefficient estimates:', sum(lasso.estimates != 0))
  
  ############################
  ## Ridge Regression Model ##
  ############################
  
  # generate ridge regression model
  pm10.ridge = glmnet(x.train, y.train, alpha=0, lambda=grid)
  ridge.cv = cv.glmnet(x.train, y.train, alpha=0)
  ridge.best_lambda = ridge.cv$lambda.min
  ridge.estimates = predict(pm10.ridge, type='coefficients', s=ridge.best_lambda)
  cat('RIDGE REGRESSION - Number of non-zero coefficient estimates:', sum(ridge.estimates != 0))
  print(ridge.estimates)
  
  #####################
  ## Get Predictions ##
  #####################
  
  # test set of predictor variables for linear models & anova model
  x.test.lm = subset(test, select=-c(log_pm10_concentration))
  
  # linear models
  lm.1.pred = predict(pm10.lm.1, new=x.test.lm)
  lm.1.error = mean((lm.1.pred - y.test)^2)
  
  lm.2.pred = predict(pm10.lm.2, new=x.test.lm)
  lm.2.error = mean((lm.2.pred - y.test)^2)
  
  lm.3.pred = predict(pm10.lm.3, new=x.test.lm)
  lm.3.error = mean((lm.3.pred - y.test)^2)
  
  # anova model
  anova.pred = predict(pm10.aov, new=x.test.lm)
  anova.error = mean((anova.pred - y.test)^2)
  
  # lasso model
  lasso.pred = predict(pm10.lasso, s=lasso.best_lambda, newx=x.test)
  lasso.error = mean((lasso.pred - y.test)^2)
  
  # ridge regression model
  ridge.pred = predict(pm10.ridge, s=ridge.best_lambda, newx=x.test)
  ridge.error = mean((ridge.pred - y.test)^2)
  
  cat('\nLinear Model #1 Error:', lm.1.error)
  cat('\nLinear Model #2 Error:', lm.2.error)
  cat('\nLinear Model #3 Error:', lm.3.error)
  cat('\nAnova Model Error:', anova.error)
  cat('\nLasso Model Error:', lasso.error)
  cat('\nRidge Regression Model Error:', ridge.error)
}

get_avg_error(0.60)
get_avg_error(0.70)
get_avg_error(0.80)


##################################
## Principal Component Analysis ##
##################################

# for plot color coding
spring = c()
fall = c()
winter = c()
for (i in 1:nrow(pm10))
{t = pm10[i,'season'][1]
if (t == 'spring'){
  spring = append(spring, i)} else if (t == 'winter'){
    winter = append(winter, i)} else{
      fall = append(fall, i)}}

# for plot color coding
morning = c()
mid_day = c()
night = c()
for (i in 1:nrow(pm10))
{t = pm10[i,'period'][1]
if (t == 'early morning' | t == 'mid morning'){
  morning = append(morning, i)} else if (t == 'mid day' | t == 'early afternoon'){
    mid_day = append(mid_day, i)} else{
      night = append(night, i)}}

# implement PCA
pm10.pca = prcomp(subset(pm10,
                         select=-c(season, period)), scale=TRUE)
# grab score vectors of first two principal components
temp = as.data.frame(pm10.pca$x[,1:2])
# plot PCA for season
plot(temp, col=ifelse(attr(temp, "row.names") %in% fall, "blue",
                      ifelse(attr(temp, "row.names") %in% winter, "red", "green")),
     main='First Two Principal Component Score Vectors',
     xlab='Principal Component 1',
     ylab='Principal Component 2',
     cex.main=0.95, cex.lab=0.8)
legend('topleft', pch=c(1,1,1), legend=c('fall', 'winter', 'spring'),
       col=c('blue', 'red', 'green'), bty='n')

# plot PCA for time of day
plot(temp, col=ifelse(attr(temp, "row.names") %in% morning, "blue",
                      ifelse(attr(temp, "row.names") %in% mid_day, "red", "green")),
     main='First Two Principal Component Score Vectors',
     xlab='Principal Component 1',
     ylab='Principal Component 2',
     cex.main=0.95, cex.lab=0.8)
legend('topleft', pch=c(1,1,1), legend=c('morning', 'mid_day', 'night'),
       col=c('blue', 'red', 'green'), bty='n')

# get importance of PCA components
pm10.components = summary(pm10.pca)$importance
pca.props = pm10.components[2, ]
plot(pca.props, xlab='PCA Component',
     ylab='Proportion of Variance',
     main='PM10 Scree Plot')
lines(pca.props, type='l')
