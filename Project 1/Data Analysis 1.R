###############################
### DATA ANALYSIS PROJECT 1 ###
### ----------------------- ###
### PM10 Pollution Analysis ###
###############################


library(HH) # for vif() function
library(lubridate) # for datetimes
library(plyr) # for mapvalues() to rename factors


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

########################
## Data Visualization ##
########################

# The 'mar' argument of 'par' sets width of the margins in the order:
# 'bottom', 'left', 'top', 'right'
par(mar=c(5,6,4,2)) # allow room for y-label

# create boxplots according to factor levels
boxplot(log_pm10_concentration ~ season, data=pm10, las=2, col=c('gray'),
        main='Concentration of PM10 Particles\nAccording to Season',
        ylab='PM10 Concentration\n(micrograms/cubic meter of air)')
par(cex.axis=0.6) # make size of x-axis labels smaller
boxplot(log_pm10_concentration ~ period, data=pm10, las=2, col=c('gray'),
        main='Concentration of PM10 Particles\nAccording to Time of Day',
        ylab='PM10 Concentration\n(micrograms/cubic meter of air)')
par(cex.axis=1.0)

# there is just one 'summer' observation so get rid of it
pm10 = pm10[which(pm10$season != 'summer'), ]

# get scatter plot matrix
pairs(pm10)

# function to generate diagnostic plots
get_diagnostic_plots = function(lm)
{
  dev.new(width=600, height=400)
  par(mfrow=c(2,3))
  plot(lm, which=c(1:6))
}


####################################
## Linear Model w/out Interaction ##
####################################

# initial model (excluding variables 'date' and 'period')
pm10.lm = lm(log_pm10_concentration ~ ., data=subset(pm10, select=-c(date, period)))
# backward stepwise regression to reduce variable set
pm10.lm.stepwise = step(pm10.lm, direction='backward', trace=0)
# check results
summary(pm10.lm.stepwise)
# get diagnostic plots
get_diagnostic_plots(pm10.lm.stepwise)


#################################
## Linear Model w/ Interaction ##
#################################

# initial model containing interaction terms (excluding datetime variable 'date')
pm10.lm = lm(log_pm10_concentration ~ (.)^2,
             data=subset(pm10, select=-c(date)))
# backward stepwise regression to reduce variable set
pm10.lm.stepwise = step(pm10.lm, direction='backward', trace=0)
# check results
summary(pm10.lm.stepwise)
# get diagnostic plots
get_diagnostic_plots(pm10.lm.stepwise)

# manually pruned model
pm10.lm.manual = lm(log_pm10_concentration ~ log_cars_per_hour + temp_at_2meters + 
                    wind_speed + season + log_cars_per_hour:wind_speed + 
                    temp_at_2meters:wind_speed, data = pm10)
# check results
summary(pm10.lm.manual)
# get diagnostic plots
get_diagnostic_plots(pm10.lm.manual)


###################
## One-Way ANOVA ##
###################

# test to see if mean pm10 concentrations are equal for all periods of the day
pm10.aov.period = aov(log_pm10_concentration ~ period, data=pm10)
summary(pm10.aov.period)
# compare all pairs of means
pm10.mmc.period = mmc(pm10.aov.period)
pm10.mmc.period

# test to see if mean pm10 concentrations are equal for all seasons
pm10.aov.season = aov(log_pm10_concentration ~ season, data=pm10)
summary(pm10.aov.season)
# compare all pairs of means
pm10.mmc.season = mmc(pm10.aov.season)
pm10.mmc.season


###################
## Two-Way ANOVA ##
###################

pm10.aov.1 = aov(log_pm10_concentration ~ period + season, data=pm10)
summary(pm10.aov.1)

pm10.aov.2 = aov(log_pm10_concentration ~ period * season, data=pm10)
summary(pm10.aov.2)

# compare all pairs of means for each anova model
pm10.mmc.1 = mmc(pm10.aov.1, focus='period')
pm10.mmc.1
pm10.mmc.2 = mmc(pm10.aov.2, focus='period')
pm10.mmc.2

