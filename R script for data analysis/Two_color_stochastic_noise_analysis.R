#----- Install packages ---------
install.packages("dplyr")
install.packages("magrittr")
install.packages("tidyverse")
install.packages("ggplot2")                
install.packages("boot")
install.packages("devtools")
install.packages("colorspace")
install.packages("RColorBrewer")
install.packages("ggthemes")
install.packages("ggpubr")
install.packages("ggsci")
install.packages("scales")
install.packages("hexbin")
install.packages("viridis")
install.packages("viridisLite")
install.packages("cowplot")
install.packages("gridExtra")
install.packages("Matching")
install.packages("fitdistrplus")
install.packages("tictoc")

#----- Load packages  ---------
library(dplyr)
library(magrittr)
library(tidyverse)
library(ggplot2)                
library(boot)
library(devtools)
library(colorspace)
library(RColorBrewer)
library(ggthemes)
library(ggpubr)
library(ggsci)
library(scales)
library(hexbin)
library(viridis)
library(viridisLite)
library(cowplot)
library(gridExtra)
library(Matching)
library (fitdistrplus)
library(tictoc)

#----- Data structure -----
# Single cell data is present in a dataframe called 'df' with columns:
# [1] "x.coord"             x-coordinate of cell centroid          
# [2] "y.coord"             y-coordinate of cell centroid
# [3] "area"                nuclear area in pixels
# [4] "mean_green_channel"  mean sfGFP signal
# [5] "stand_dev_green"     std dev of sfGFP signal
# [6] "mean_DAPI_channel"   mean DAPI signal
# [7] "stand_dev_DAPI"      std dev of DAPI signal
# [8] "mean_red_channel"    mean mCherry signal
# [9] "stand_dev_red"       std dev of mCherry signal
# [10] "Disc_ID"            wing disc ID no. 
# [11] "Slice"              optical slice depth 
# [12] "Set"                experiment no. (unique) 
# [13] "Genotype"           cell genotype 
# [14] "Unique_ID"          unique ID number for each disc across experiments 

# Each row in the dataframe represents a single segmented nucleus

#----- Clean up segmentation artefacts -----
# remove segmented cells with area < 100 to ensure signal intensity was measured over at least n = 100 pixels
# remove cells with area > 4000 - cells fused incorrectly during segmentation
df <- subset(df, mean_green_channel > 0 & mean_red_channel > 0 & area > 100 & area < 4000)

#----- Check for experimental errors -------
# Check that all discs show flourescence signal in both channels 
ggplot()+ theme_bw() +
geom_point(data=df, aes(x=mean_green_channel, y=mean_red_channel), size = 1, col='darkblue') + 
facet_wrap(~Unique_ID, scales="fixed")

#----- Background Subtraction : prepare data to fit a normal curve to the background signal ------
# define Gaussian fitting function 'fitG' to optimally fit background signal 
fitG = function (x,y,mu,sig,scale) {
  f = function(p){
    d = p[3]*dnorm(x,mean=p[1], sd=p[2])
    sum((d-y)^2)
  }
  optim(c(mu,sig,scale),f)
}

# pick a suitable RFU cut-off before fitting background signal
# this is to ensure that extremely large true signal values don't skew optimization for background mean
sfGFP_cutoff = 3.5
mCherry_cutoff = 3
plot(density(df$mean_green_channel), col='green', xlim=c(0,10))
lines(density(df$mean_red_channel), col='red')
abline(v=c(sfGFP_cutoff, mCherry_cutoff), col =c("darkgreen","red3"), lty=2, lwd=2)

# background signal peak should be below chosen cut-off value for each disc
# check for sfGFP channel
ggplot()+theme_bw() + 
geom_density(data=df, aes(x=mean_green_channel), size = 1, col='green') + 
facet_wrap(~Unique_ID, scales="fixed") + 
scale_x_continuous(limits = c(0,5)) + 
geom_vline(xintercept = sfGFP_cutoff)

# check for mCherry channel
ggplot()+theme_bw() + 
geom_density(data=df, aes(x=mean_red_channel), size = 1, col='red3') + 
facet_wrap(~Unique_ID, scales="fixed") + 
scale_x_continuous(limits = c(0,5)) + 
geom_vline(xintercept = mCherry_cutoff)

#----- Background Subtraction : calculate Normal fit parameters for each disc ------
#define variable vectors to store mean and 84th percentile of Gaussian fit
p84_green = p84_red = c(NULL)
mu_green = mu_red = c(NULL)

# vector to store unique discID
discs = levels(df$Unique_ID) 

# loop through to fit both channels separately for each disc
for (n in seq_along(discs)){
  density <- density(collect(filter(df, Unique_ID == discs[n] , mean_green_channel < sfGFP_cutoff))$mean_green_channel)
  fitg = fitG(density$x, density$y,2,1,0.02)
  p = fitg$par
  mu_green[n] = p[1]                 # mean green channel background
  p84_green[n] = p[1] + (0.995*p[2]) # 84th percentile RFU for background green signal
   
  density <- density(collect(filter(df, Unique_ID == discs[n]  , mean_red_channel < mCherry_cutoff))$mean_red_channel)
  fitr = fitG(density$x, density$y,10,4,0.02)
  q = fitr$par
  mu_red[n] = q[1]                  # mean red channel background
  p84_red[n] = q[1] + (0.995*q[2])  # 84th percentile RFU for background red signal
}

# Check the distribution of mean background levels across discs - should be similar to each other
ggplot() + theme_bw() +
geom_jitter(aes(x=1, y = mu_red), col ='red3') +  # red background means
geom_jitter(aes(x=2, y = mu_green), col ='green') # green background means

# add a column to the dataframe which specifies the mean and 84th percentile background signal
# this is calculated separately for each disc and each channel
# thus cells belonging to the same disc are background corrected identically
pos <- list()
for(n in seq_along(discs)){
  # subset the data by Unique_ID for every disc to add columns with calculated fit values
  pos[[n]] <- subset(df, df$Unique_ID == discs[n])
  pos[[n]]$mu_red <- mu_red[n]  
  pos[[n]]$mu_green <- mu_green[n]
  pos[[n]]$p84_red <- p84_red[n]
  pos[[n]]$p84_green <- p84_green[n]
}
# bind data from multiple discs back into a single dataframe
df <- do.call(rbind.data.frame, pos)
df <- tbl_df(df)

# Assuming all cells below 84th percentile are Sens-negative
# Subtract 'background' signal (84th percentile of Gaussian fit) from raw signal
df$p84.red.subtracted <- df$mean_red_channel - df$p84_red
df$p84.green.subtracted <- df$mean_green_channel - df$p84_green

# discard cells that are below channel cut-offs for analysis
df = collect(filter(df, p84.green.subtracted > 0 , p84.red.subtracted > 0))

#----- Channel Normalization : Scale mCherry and sfGFP signal intensity units to be similar ------
# calculate linear fit coefficients for green ~ red signal for individual discs
# slope (alpha) and intercept 
alpha = intercept =  c()
for(i in seq_along(discs)) {
  disc_i = collect(filter(df, Unique_ID == discs[i]))
  fit_i = coef(lm(p84.green.subtracted ~ p84.red.subtracted, data = disc_i))
  alpha[i] <- fit_i[[2]]
  intercept[i] <- fit_i[[1]]
}

# add disc-specific estimates of slope and intercept to each cell in that disc
pos <- list()
for(n in seq_along(discs)){
  #adding columns with disc-wise values filled in
  pos[[n]] <- subset(df, df$Unique_ID == discs[n])
  pos[[n]]$alpha <- alpha[n]  
  pos[[n]]$intercept <- intercept[n]
}
df <- do.call(rbind.data.frame, pos)
df <- tbl_df(df)

# transform red signal values based on fit coefficients
df$red.prime <- (df$alpha * df$p84.red.subtracted) + df$intercept

# assign new columns which contain background and slope corrected values
df$normalized_green <- df$p99.green.subtracted
df$normalized_red <- df$red.prime 

# calculate total protein
df$total <- (df$normalized_green + df$normalized_red)
df <- df[df$total > 0, ] # safety check to ensure all non-negative values

# save transformed data
write.csv(file="yy-mm-dd 84% scaled formatted data", x=df, row.names = F)

#----- Fano factor calculation ------
# define function to calculate fano factor
fano_factor <- function (data, indices){
  d = data[indices,]
  x = d$normalized_green
  y = d$normalized_red
  sq_diff = (x - y)^2
  mean_sq_diff = mean(sq_diff)
  normalization = 2*mean(x)*mean(y)
  fractional_variance = mean_sq_diff/normalization
  total_protein = mean(x) + mean(y)
  fano_factor = fractional_variance*total_protein
  return(fano_factor)
}

# choose genotype for analysis
df$Genotype = as.factor(df$Genotype)
summary(df$Genotype)
selected = c("22A3-mutant") #choose dataset for analysis
df = collect(filter(df, Genotype %in% selected)) # cells of chosen genotype 

# sort cells into bins according to total protein level
n_bins = 600 # choose total number of bins
numseq <- seq(-2,4,length.out = n_bins) # vector of bin boundaries (log RFU)

df_bin = list(NULL)  # list of dataframes, each dataframe is one binned population of cells
sample_size_bin = c(NULL) # to store the population size of each bin

for(i in 1:(length(numseq)-1)){
  df_bin[[i]] <-  df[log(df$total) > numseq[i] & log(df$total) < numseq[i+1],]
  sample_size_bin[i] = nrow(df_bin[[i]])
}

# print out sample sizes to check n values for each bin
sample_size_bin 
# select suitable bins such that sample size > minimum number of cells
minimum_sample_size = 20
bin_range = which(sample_size_bin > minimum_sample_size)
# bin_range contains the bin indices which will be used for futher analysis

# calculate fano factor and mean protein level for each bin
protein_bin = fano_bin = c(NULL)

for (i in seq_along(bin_range)){
  fano_bin[i] = fano_factor(df_bin[[bin_range[i]]])
  protein_bin[i] = mean(df_bin[[bin_range[i]]]$total)
}
# This is a singular metric estimated from the binned population
# To understand the variance around this metric, we must resample cells and estimate it multiple times

#----- Boostrap to estimate confidence intervals for Fano factor ----------
numberOfResamples = 5000 # choose appropriate number
Fano = upr_bound = lwr_bound = Protein = Sample_size = medianFano = SEboot =c(NULL)

tic() # start timer to measure run time 
for(i in seq_along(bin_range)){
  print(i) # track loop progress
  bootobject = boot(data = df_bin[[bin_range[i]]], statistic = fano_factor, R = numberOfResamples, parallel = 'multicore')
  ci = boot.ci(bootobject, type="perc") #percentile method of estimating CI
  Fano[i] = ci$t0 #raw value from real dataset - without bootstrapping
  upr_bound[i] = ci$percent[1,5]
  lwr_bound[i] = ci$percent[1,4]
  Protein[i] = mean(df_bin[[bin_range[i]]]$total)
  Sample_size[i] = nrow(df_bin[[bin_range[i]]])
  medianFano[i] = median(bootobject$t) # calculates median of resamples
  SEboot[i] = sd(bootobject$t) # calculates S.E.M.
  # Check that bootstrap is sufficiently resampled - should yeild a normal distribution
  # skewed histograms indicate either sample sizes are too low or not enough re-sampling
  gghistogram(bootobject$t, bins=50)
}
toc() # end timer

# save bootstrapped Fano estimates
names  <- c("Protein","Fano","lwr_bound","upr_bound","Sample_size","medianFano","SEBoot","Genotype")
booted_bin_metrics = cbind(Protein, Fano, lwr_bound, upr_bound, Sample_size, medianFano, SEboot)
colnames(booted_bin_metrics) = names
booted_bin_metrics <- as.data.frame(booted_bin_metrics)
booted_bin_metrics$Genotype = selected # add a column to identify cell genotype
write.csv(file = 'yy-mm-dd Bootstrapped Bin Fano Factor for genotypeX.csv', x = booted_bin_metrics, row.names=F)

#----- Technical noise baseline subtraction ------
# load dataset for Tandem-tagged sfGFP-mCherry-Sens cells
technical_control = read.csv(file = 'yy-mm-dd Bootstrapped Bin Fano Factor for Tandem-tag-sens.csv', header = T)

# fit locally weighted polynomial to Fano factor vs Protein tandem data
loess_fit = loess(Fano ~ Protein, data = technical_control)

# predict the technical noise for every protein level in test genotype using previous fit coefficients
test_data$Predicted_baseline = predict(loess_fit, test_data$Protein)

# subtract technical noise to get corrected Fano factor estimates (Fano Factor is additive, so : total Fano factor = technical Fano factor + biological Fano factor)
test_data$Corrected_Fano = test_data$Fano - test_data$Predicted_baseline
test_data$Corrected_upr_bound = test_data$upr_bound - test_data$Predicted_baseline
test_data$Corrected_lwr_bound = test_data$lwr_bound - test_data$Predicted_baseline

#----- Convert RFU to molecule numbers and plot data -------------------------
# pick dataset and define conversion factor
plot_data = test_data # select dataset to plot
conv_factor = 138 # 1 RFU = 138 molecules

# scatter plot for sfGFP and mCherry signal correlation
  ggplot(data = plot_data, aes(x=normalized_green*conv_factor, y=normalized_red*conv_factor, col = as.factor(Genotype))) + 
  geom_point(size=0.5) + 
  theme_pubr(base_family = "Helvetica", base_size = 14)  + 
  scale_color_aaas() + 
  scale_x_continuous(limits = c(0,5000)) + xlab("sfGFP-Sens (molecules)") +
  scale_y_continuous(limits = c(0,5000)) + ylab("mCherry-Sens (molecules)") +
  coord_fixed()

# Fano factor plot
  ggplot(data = plot_data, aes(x = Protein*conv_factor, y = Corrected_Fano*conv_factor , col = as.factor(Genotype))) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'loess', span = 0.3, alpha = 0.15, lwd = 0.75, fill = 'grey') + 
  theme_pubr(base_family = "Helvetica", base_size = 14) +
  scale_color_aaas() +
  scale_x_continuous(limits = c(0,1500)) + xlab("Sens protein (molecules)") +
  scale_y_continuous(limits = c(0,150)) + ylab(expression("Fano factor (molecules)  " (sigma^{2}/mu)))
  
  