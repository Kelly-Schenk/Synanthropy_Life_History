#--------------------------------------------------------------------------------------------------------------------#
                                                  #DATA NOTES#
##Source for data: Ecke et al., 2022 (synanthropy categorization); Heldstab, 2021 (life history variables); and Shultz & Dunbar, 2010 (relative brain size).##
##Brain and body data are from separate data sets than the rest of the data, and thus will be analyzed together only.## 
##The generative artificial intelligence (GenAI) tool Microsoft Copilot was used to develop and correct code in this study; any code was reviewed for accuracy and correctness, as well as tested following its generation.##
##Code initially written for Penn State University's ANTH 571 course, Spring 2026 by Kelly Elizabeth Schenk - questions on code can be directed to kes6724@psu.edu. Please e-mail this address to report any issues with the code - I am still learning!##
##Code is subject to change.##
#--------------------------------------------------------------------------------------------------------------------#
                                                #INITIALIZING#

##Set libraries, includes standard packages used across my research##
  ###Core tidyverse packages (includes ggplot2, dplyr, tidyr, etc.)
  library(tidyverse)
  
  ###Data manipulation and performance
  library(data.table)
  
  ###Development tools
  library(devtools)
  
  ###Regression and statistical analysis
  library(olsrr)
  library(lmtest)
  library(nortest)
  library(car)
  library(rsq)
  library(groupcompare)
  library(MASS)
  library(quantreg)
  
  ###Visualization
  library(factoextra)
  library(plotly)
  library(ggiraph)
  library(patchwork)
  library(ggthemes)
  library(corrplot)
  
  ###Spatial data
  library(sf)
  
  ###Miscellaneous utilities
  library(Hmisc)
  
#--------------------------------------------------------------------------------------------------------------------#
                                      #DATA UPLOAD AND MANIPULATION#
##Set working directory##
  setwd("C:/Users/Kelly Elizabeth/Documents/PHD/Courses/Spring 2026/ANTH 571/R/Working Directory")
##Pulling data from computer into RStudio##
  data = read.table('Rodents_R_final.csv', sep=',', header=T)
  summary(data)
  data.frame(data)
  data$adult_mass <- as.numeric(data$adult_mass)

##Subset data by synanthropic and non-synanthropic##
  subset_anth <- data %>% group_by(Synanthropic)
  split_data <- group_split(subset_anth)
  nonsyn <- print(split_data[[1]])
  syn <- print(split_data[[2]])
  summary(nonsyn)
  summary(syn)

##Subset data by reservoir status##
  subset_path <- data %>% group_by(reservoir)
  split_data_2 <- group_split(subset_path)
  nonres <- print(split_data_2[[1]])
  res <- print(split_data_2[[2]])
  summary(nonres)
  summary(res)
  
#--------------------------------------------------------------------------------------------------------------------#
                                                #ANALYSIS#
#-------Assessing Data-------#
  boxplot(data$adult_mass,
          xlab = "Adult Mass (g)",
          col = "lightgreen",
          main = "Distribution of Adult Mass") #has outliers, but will not correct for them
  boxplot(data$gestation,
          xlab = "Gestation (in days)",
          col = "lightyellow",
          main = "Distribution of Gestation") #no outliers
  boxplot(data$litters_annually,
          xlab = "Litters Annually (in number of offspring per litter)",
          col = "blue",
          main = "Distribution of Number of Offspring") #has 1 outlier, but will not correct for it
  boxplot(data$litter_size,
          xlab = "Litter Size (in number of litters per year)",
          col = "deeppink",
          main = "Distribution of Number of Litters per Year") #has 1 outlier, but will not correct for it
  boxplot(data$brain,
          xlab = "Brain",
          col = "violet",
          main = "Distribution of Brain") #no outliers
  boxplot(data$body,
          xlab = "Body",
          col = "skyblue",
          main = "Distribution of Body") #no outliers

##Graphing Relationships##
  plot(data$adult_mass,data$gestation, xlab="Adult Body Mass (g)", ylab="Gestational Period (days)", main="Adult Body Mass and Gestational Period")
  plot(data$adult_mass,data$litter_size, xlab="Adult Body Mass (g)", ylab="Litter Size (Offspring Quantity)", main="Adult Body Mass and Litter Size")
  plot(data$adult_mass,data$litters_annually, xlab="Adult Body Mass (g)", ylab="Litters Produced Annually", main="Adult Body Mass and Litters Annually")
  plot(data$brain,data$body, xlab="Brain", ylab="Body", main="Brain-Body") #Brain and body data are from separate data sets than the rest of the data, and thus will be analyzed together only#
  
#-------Pearson Correlation-------#
  data_frame_set <- data.frame(data$adult_mass, data$gestation,data$litters_annually, data$litter_size)
  resultfull <- cor(data_frame_set, method = "pearson", use = "complete.obs")
    print(resultfull) 
    
  result1 <- cor(data$adult_mass, data$gestation, method = "pearson") #significant
    cat("Pearson correlation coefficient is:", result1)
    result_stat1 <- cor.test(data$adult_mass, data$gestation, method = "pearson")
    print(result_stat1) 
  
  result2 <- cor(data$adult_mass, data$litters_annually, method = "pearson") #significant
    cat("Pearson correlation coefficient is:", result2)
    result_stat2 <- cor.test(data$adult_mass, data$litters_annually, method = "pearson")
    print(result_stat2)
  
  result3 <- cor(data$adult_mass, data$litter_size, method = "pearson") #ns
    cat("Pearson correlation coefficient is:", result3)
    result_stat3 <- cor.test(data$adult_mass, data$litter_size, method = "pearson")
    print(result_stat3)
  
  resultx <- cor(data$gestation, data$litters_annually, method = "pearson") #ns
    cat("Pearson correlation coefficient is:", resultx)
    result_statx<- cor.test(data$gestation, data$litters_annually, method = "pearson")
    print(result_statx)
  
  resultx2 <- cor(data$gestation, data$litter_size, method = "pearson") #significant
    cat("Pearson correlation coefficient is:", resultx2)
    result_statx2 <- cor.test(data$gestation, data$litter_size, method = "pearson")
    print(result_statx2)
  
  resultx3 <- cor(data$litter_size, data$litters_annually, method = "pearson") #ns
    cat("Pearson correlation coefficient is:", resultx3)
    result_statx3 <- cor.test(data$litter_size, data$litters_annually, method = "pearson")
    print(result_statx3)
    
  resultxx <- cor(data$body, data$brain, method = "pearson") 
    cat("Pearson correlation coefficient is:", resultxx)
    result_statxx <- cor.test(data$body, data$brain, method = "pearson")
    print(result_statxx)
    
  ###Visuals###
    data_frame_pearson <- data.frame(data$adult_mass, data$gestation, data$litter_size, data$litters_annually)
      cor_results <- rcorr(as.matrix(data_frame_pearson), type = "pearson")
      cor_matrix <- cor_results$r   # Correlation coefficients
      p_matrix   <- cor_results$P   # P-values
    labelnames <- c("Adult Mass",
                    "Gestation Length",
                    "Litter Size",
                    "Litters Annually")
    corrplot(
      cor_matrix,
      method = "color",       # Color-coded heatmap
      order = "hclust",       # Cluster similar variables
      addCoef.col = "black",  # Add correlation coefficients
      tl.col = "black",       # Text label color
      tl.srt = 45,            # Rotate labels
      p.mat = p_matrix,       # Matrix of p-values
      sig.level = 0.05,       # Mark significant correlations
    )

#-------LINEAR REGRESSION-------#
##Running simple linear regressions on whole data set, piecewise##
  ###Adult Mass v Gestation Length###
  adult_mass_log <- log(data$adult_mass) #Log transformations to correct for non-normality#
  gestation_log <- log(data$gestation)
  model = lm(gestation_log~adult_mass_log, data=data)
    summary(model)
    anova(model) #Significant#
    
    ###Plot###
    ggplot(data, aes(x = exp(adult_mass_log), y = exp(gestation_log))) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Average Adult Mass vs. Average Gestation Length",
           x = "Average Adult Mass (in grams)",
           y = "Average Gestation Length (in days)") +
        theme(
        plot.title = element_text(hjust = 0.5)) +
      scale_fill_stata(scheme="s2color")
  
      ###Checking assumptions###
      #Normality# - Not normal
      hist(model$residuals, main="Histogram of Residuals") #data name precedes $
      qqnorm(model$residuals, main="QQ Plot") #normality plot residuals
      qqline(model$residuals) #want all or majority points fall near line
      ad.test(model$residuals) #Anderson-Darling normality test
      shapiro.test(model$residuals) #Shapiro-Wilk normality test
      lillie.test(model$residuals) #Kolmogorov-Smirnov normality test
      #Constant Residual Variance#
      plot(model$fitted,model$residuals, ylab="Residuals", xlab="Fitted Values")
      abline(0,0)
      bptest(model) #Breusch-Pagan for constant variance
      #Linearity#
      abline(model, col = "blue")
      cor.test(data$adult_mass,data$litters_annually)
      #Independence# - Not independent
      plot(model$residuals, main="Residuals vs Order")
      acf(model$residuals)
      dwtest(model)
  
  ###Adult Mass v Litter Size###
  litter_size_log <- log(data$litter_size)
  model2 = lm(litter_size_log~adult_mass_log, data=data) 
    summary(model2)
    anova(model2) #ns#
    
    ###Plot###
    ggplot(data, aes(x = adult_mass, y = litter_size)) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Average Adult Mass vs. Average Litter Size",
           x = "Average Adult Mass (in grams)",
           y = "Average Litter Size (number of offspring)") +
      theme(
        plot.title = element_text(hjust = 0.5))
    
  
  ###Adult Mass v Litters Annually###
  litters_annually_log <- log(data$litters_annually) #Log transformations to correct for non-normality#
  adult_mass_log <- log(data$adult_mass)
  model3 = lm(litters_annually_log~adult_mass_log, data=data)
    summary(model3)
    anova(model3) #Significant#
    
    ###Plot###
    ggplot(data, aes(x = exp(adult_mass_log), y = exp(litters_annually_log))) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Average Adult Mass vs. Average Litters Annually",
           x = "Average Adult Mass (in grams)",
           y = "Average Litters Annually (number of litters)") +
      theme(
        plot.title = element_text(hjust = 0.5))
  
      ###Checking assumptions###
      #Normality#
      hist(model3$residuals, main="Histogram of Residuals") #data name precedes $
      qqnorm(model3$residuals, main="QQ Plot") #normality plot residuals
      qqline(model3$residuals) #want all or majority points fall near line
      ad.test(model3$residuals) #Anderson-Darling normality test
      shapiro.test(model3$residuals) #Shapiro-Wilk normality test
      lillie.test(model3$residuals) #Kolmogorov-Smirnov normality test
      #Constant Residual Variance#
      plot(model3$fitted,model3$residuals, ylab="Residuals", xlab="Fitted Values")
      abline(0,0)
      bptest(model3) #Breusch-Pagan for constant variance
      #Linearity#
      abline(model3, col = "blue")
      cor.test(data$adult_mass,data$litters_annually)
      #Independence#
      plot(model3$residuals, main="Residuals vs Order")
      acf(model3$residuals)
      dwtest(model3)
      
  ###Gestation Length v Litters Annually###
  modelx1 = lm(litters_annually_log~gestation_log, data=data)
    summary(modelx1)
    anova(modelx1) #ns#
    
    ###Plot###
    ggplot(data, aes(x = exp(litters_annually_log), y = exp(gestation_log))) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Average Gestation Length vs. Average Litters Annually",
           x = "Average Gestation Length (in days)",
           y = "Average Litters Annually (number of litters)") +
      theme(
        plot.title = element_text(hjust = 0.5))
  
  ###Gestation Length v Litter Size###
  modelx2 = lm(litter_size_log~gestation_log, data=data)
    summary(modelx2)
    anova(modelx2) #Significant#
    
    ###Checking assumptions###
    #Normality# - Not normal
    hist(modelx2$residuals, main="Histogram of Residuals") #data name precedes $
    qqnorm(modelx2$residuals, main="QQ Plot") #normality plot residuals
    qqline(modelx2$residuals) #want all or majority points fall near line
    ad.test(modelx2$residuals) #Anderson-Darling normality test
    shapiro.test(modelx2$residuals) #Shapiro-Wilk normality test
    lillie.test(modelx2$residuals) #Kolmogorov-Smirnov normality test
    #Constant Residual Variance#
    plot(modelx2$fitted,modelx2$residuals, ylab="Residuals", xlab="Fitted Values")
    abline(0,0)
    bptest(modelx2) #Breusch-Pagan for constant variance
    #Linearity#
    abline(modelx2, col = "blue")
    cor.test(data$brain,data$body)
    #Independence#
    plot(modelx2$residuals, main="Residuals vs Order")
    acf(modelx2$residuals)
    dwtest(modelx2)
    
    ###Plot###
    ggplot(data, aes(x = exp(litter_size_log), y = exp(gestation_log))) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Average Gestation Length vs. Average Litter Size",
           x = "Average Gestation Length (in days)",
           y = "Average Litter Size (number of offspring)") +
      theme(
        plot.title = element_text(hjust = 0.5))
  
  ###Body v Brain###
  model4 = lm(body~brain, data=data)
    summary(model4)
    anova(model4) #Significant#
    
    ###Plot###
    ggplot(data, aes(x = body, y = brain)) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Relative Body Size vs. Relative Brain Size",
           x = "Relative Body Size",
           y = "Relative Brain Size") +
      theme(
        plot.title = element_text(hjust = 0.5))
  
      ###Checking assumptions###
      #Normality# - Not normal
      hist(model4$residuals, main="Histogram of Residuals") #data name precedes $
      qqnorm(model4$residuals, main="QQ Plot") #normality plot residuals
      qqline(model4$residuals) #want all or majority points fall near line
      ad.test(model4$residuals) #Anderson-Darling normality test
      shapiro.test(model4$residuals) #Shapiro-Wilk normality test
      lillie.test(model4$residuals) #Kolmogorov-Smirnov normality test
      #Constant Residual Variance#
      plot(model4$fitted,model4$residuals, ylab="Residuals", xlab="Fitted Values")
      abline(0,0)
      bptest(model4) #Breusch-Pagan for constant variance
      #Linearity#
      abline(model4, col = "blue")
      cor.test(data$brain,data$body)
      #Independence#
      plot(model4$residuals, main="Residuals vs Order")
      acf(model4$residuals)
      dwtest(model4)

##Running multiple linear regressions on significant data, primarily to determine if gestation and adult mass together can better predict reproductive variables##
  ###Adult Mass + Gestation vs Litter Size###
  litter_size_log <- log(data$litter_size) #Log transformations to correct for non-normality#
  mlmodel = lm(litter_size_log~adult_mass_log+gestation_log, data=data)
    summary(mlmodel)
    anova(mlmodel)
    Anova(mlmodel, type=3) #Significant#
    
    ###Checking assumptions###
    #Normality#
    hist(mlmodel$residuals, main="Histogram of Residuals") #data name precedes $
    qqnorm(mlmodel$residuals, main="QQ Plot") #normality plot residuals
    qqline(mlmodel$residuals) #want all or majority points fall near line
    ad.test(mlmodel$residuals) #Anderson-Darling normality test
    shapiro.test(mlmodel$residuals) #Shapiro-Wilk normality test
    lillie.test(mlmodel$residuals) #Kolmogorov-Smirnov normality test
    #Constant Residual Variance#
    plot(mlmodel$fitted,mlmodel$residuals, ylab="Residuals", xlab="Fitted Values")
    abline(0,0)
    bptest(mlmodel) #Breusch-Pagan for constant variance
    #Linearity#
    abline(mlmodel, col = "blue")
    #Independence#
    plot(mlmodel$residuals, main="Residuals vs Order")
    acf(mlmodel$residuals)
    dwtest(mlmodel)
    
    ###Checking for multicollinearity###
    vif(lm(litter_size_log~adult_mass_log+gestation_log, data=data))
    
  ###Adult Mass + Gestation vs Litters Annually###
  litters_annually_log <- log(data$litters_annually) #Log transformations to correct for non-normality#
    mlmodelx = lm(litters_annually_log~adult_mass_log+gestation_log, data=data)
    summary(mlmodelx)
    anova(mlmodelx)
    Anova(mlmodelx, type=3) #ns#
    
#--------------------------------------------------------------------------------------------------------------------#
                                  #SUB-ANALYSIS: BY SYNANTHROPY CATEGORIZATION#
    
#-------Pearson Correlation-------#
##By Synanthropy Categorization##
  result4 <- cor(syn$adult_mass, syn$gestation, method = "pearson") #ns
    cat("Pearson correlation coefficient is:", result4)
    result_stat4 <- cor.test(syn$adult_mass, syn$gestation, method = "pearson")
    print(result_stat4)
  
  result5 <- cor(nonsyn$adult_mass, nonsyn$gestation, method = "pearson") #significant
    cat("Pearson correlation coefficient is:", result5)
    result_stat5 <- cor.test(nonsyn$adult_mass, nonsyn$gestation, method = "pearson")
    print(result_stat5)
  
  result6 <- cor(syn$adult_mass, syn$litters_annually, method = "pearson") #significant
    cat("Pearson correlation coefficient is:", result6)
    result_stat6 <- cor.test(nonsyn$adult_mass, nonsyn$gestation, method = "pearson")
    print(result_stat6)
  
  result7 <- cor(nonsyn$adult_mass, nonsyn$litters_annually, method = "pearson") #ns
    cat("Pearson correlation coefficient is:", result7)
  
  result8 <- cor(syn$adult_mass, syn$litter_size, method = "pearson") #ns
    cat("Pearson correlation coefficient is:", result8)
  
  result9 <- cor(nonsyn$adult_mass, nonsyn$litter_size, method = "pearson") #ns
    cat("Pearson correlation coefficient is:", result9)
    
  resultx4 <- cor(nonsyn$gestation, nonsyn$litter_size, method = "pearson") #significant
    cat("Pearson correlation coefficient is:", resultx4)
    result_stat4 <- cor.test(nonsyn$gestation, nonsyn$litter_size, method = "pearson")
    print(result_stat4)
    
  resultx5 <- cor(syn$gestation, syn$litter_size, method = "pearson") #significant
    cat("Pearson correlation coefficient is:", resultx5)
    result_stat5 <- cor.test(syn$gestation, syn$litter_size, method = "pearson")
    print(result_stat5)
    
  resultx6 <- cor(nonsyn$gestation, nonsyn$litters_annually, method = "pearson") #na
    cat("Pearson correlation coefficient is:", resultx6)
  
  resultx7 <- cor(syn$gestation, syn$litters_annually, method = "pearson") #ns 
    cat("Pearson correlation coefficient is:", resultx7)
    result_statx7 <- cor.test(syn$gestation, syn$litters_annually, method = "pearson")
    print(result_statx7)
  
  resultx8 <- cor(nonsyn$litter_size, nonsyn$litters_annually, method = "pearson") #na
    cat("Pearson correlation coefficient is:", resultx8)
    
  resultx9 <- cor(syn$litter_size, syn$litters_annually, method = "pearson") #ns
    cat("Pearson correlation coefficient is:", resultx9)
    result_statx9 <- cor.test(syn$litter_size, syn$litters_annually, method = "pearson")
    print(result_statx9)
  
##Checkng conditions for T-tests and Whitney-Mann U as appropriate##
  summary(syn)
    shapiro.test(syn$adult_mass) #not normal
    shapiro.test(syn$litter_size) #not normal
    shapiro.test(syn$litters_annually) #not normal
    shapiro.test(syn$gestation) #not normal
    shapiro.test(syn$body)
    shapiro.test(syn$brain) 
  summary(nonsyn)
    shapiro.test(nonsyn$adult_mass) #not normal
    shapiro.test(nonsyn$litter_size) 
    shapiro.test(nonsyn$litters_annually) #not normal
    shapiro.test(nonsyn$gestation) #not normal
    shapiro.test(nonsyn$body) 
    shapiro.test(nonsyn$brain) #not normal
  ###Variance Tests###  
    var_test1 <- var.test(syn$adult_mass, nonsyn$adult_mass)
    print(var_test1)  #variance equality
    var_test2 <- var.test(syn$litter_size, nonsyn$litter_size)
    print(var_test2) #variance equality
    var_test3 <- var.test(syn$litters_annually, nonsyn$litters_annually)
    print(var_test3) #variance equality
    var_test4 <- var.test(syn$gestation, nonsyn$gestation)
    print(var_test4) #variance equality
    var_test5 <- var.test(syn$body, nonsyn$body)
    print(var_test5)
    var_test6 <- var.test(syn$brain, nonsyn$brain)
    print(var_test6) #no variance equality
  
##Mann-Whitney U Tests##
  wilcox_test_0 <- wilcox.test(syn$adult_mass, nonsyn$adult_mass)
    print(wilcox_test_0) #ns
  
  wilcox_1 <- wilcox.test(syn$litter_size, nonsyn$litter_size)
    print(wilcox_1) #ns
  
  wilcox_ <- wilcox.test(syn$litters_annually, nonsyn$litters_annually)
    print(wilcox_) #ns
  
  wilcox_2 <- wilcox.test(syn$gestation, nonsyn$gestation)
    print(wilcox_2) #ns
  
  wilcox_3 <- wilcox.test(syn$body, nonsyn$body)
    print(wilcox_3) #ns
  
  wilcox_4 <- wilcox.test(syn$brain, nonsyn$brain)
    print(wilcox_4) #ns

###Repeating plots above but subset###
  plot(syn$adult_mass,syn$gestation, xlab="Adult Body Mass (g)", ylab="Gestational Period (days)", main="Adult Body Mass and Gestational Period - Synanthropic")
  plot(nonsyn$adult_mass,nonsyn$gestation, xlab="Adult Body Mass (g)", ylab="Gestational Period (days)", main="Adult Body Mass and Gestational Period - Nonsynanthropic")
  plot(syn$adult_mass,syn$litter_size, xlab="Adult Body Mass (g)", ylab="Litter Size (Offspring Quantity)", main="Adult Body Mass and Litter Size - Synanthropic")
  plot(nonsyn$adult_mass,nonsyn$litter_size, xlab="Adult Body Mass (g)", ylab="Litter Size (Offspring Quantity)", main="Adult Body Mass and Litter Size - Nonsynanthropic")
  plot(syn$adult_mass,syn$litters_annually, xlab="Adult Body Mass (g)", ylab="Litters Annually (Litter Quantity)", main="Adult Body Mass and Litters Annually - Synanthropic")
  plot(nonsyn$adult_mass,nonsyn$litters_annually, xlab="Adult Body Mass (g)", ylab="Litters Annually (Litter Quantity)", main="Adult Body Mass and Litters Annually - Nonsynanthropic")
  plot(syn$body,syn$brain, xlab="Body", ylab="Brain", main="Body-Brain: Synanthropic")
  plot(nonsyn$body,nonsyn$brain, xlab="Body", ylab="Brain", main="Body-Brain: Non-synanthropic")

#-------LINEAR REGRESSION BY SYNANTHROPY CATEGORIZATION-------#
  ###Adult Mass v Litters Annually, Nonsynanthropic###
  model6 = lm(litters_annually~adult_mass, data=nonsyn) #ns#
    summary(model6)
    anova(model6) 
    
    ###Plot###
    ggplot(nonsyn, aes(x = adult_mass, y = litters_annually)) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Average Adult Mass vs. Average Litters Annually, Non-synanthropic Rodents",
           x = "Average Adult Mass (in grams)",
           y = "Average Litters Annually") +
      theme(
        plot.title = element_text(hjust = 0.5))
  
  ###Adult Mass v Litters Annually, Synanthropic###
  litters_annually_log_syn <- log(syn$litters_annually) #Log transformations to correct for non-normality#
  adult_mass_log_syn <- log(syn$adult_mass)
  model7 = lm(litters_annually_log_syn~adult_mass_log_syn, data=syn) #Significant#
    summary(model7)
    anova(model7) 
    
    ###Plot###
    ggplot(syn, aes(x = exp(adult_mass_log_syn), y = exp(litters_annually_log_syn))) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Average Adult Mass vs. Average Litters Annually, Synanthropic Rodents",
           x = "Average Adult Mass (in grams)",
           y = "Average Litters Annually") +
      theme(
        plot.title = element_text(hjust = 0.5))
    
      ###Checking assumptions###
      #Normality# - Not normal distribution
      hist(model7$residuals, main="Histogram of Residuals") #data name precedes $
      qqnorm(model7$residuals, main="QQ Plot") #normality plot residuals
      qqline(model7$residuals) #want all or majority points fall near line
      ad.test(model7$residuals) #Anderson-Darling normality test
      shapiro.test(model7$residuals) #Shapiro-Wilk normality test
      lillie.test(model7$residuals) #Kolmogorov-Smirnov normality test
      #Constant Residual Variance#
      plot(model7$fitted,model7$residuals, ylab="Residuals", xlab="Fitted Values")
      abline(0,0)
      bptest(model7) #Breusch-Pagan for constant variance
      #Linearity#
      abline(model7, col = "blue")
      cor.test(syn$adult_mass,syn$litters_annually)
      #Independence#
      plot(model7$residuals, main="Residuals vs Order")
      acf(model7$residuals)
      dwtest(model7)
  
  ###Adult Mass v Gestation Length, Nonsynanthropic###
  gestation_log_nonsyn <- log(nonsyn$gestation) #Log transformations to correct for non-normality#
  adult_mass_log_nonsyn <- log(nonsyn$adult_mass)
  model8 = lm(gestation_log_nonsyn~adult_mass_log_nonsyn, data=nonsyn) #Significant#
    summary(model8)
    anova(model8) 
    
    ###Plot###
    ggplot(nonsyn, aes(x = exp(adult_mass_log_nonsyn), y = exp(gestation_log_nonsyn))) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Average Adult Mass vs. Average Gestation Length, Non-synanthropic Rodents",
           x = "Average Adult Mass (in grams)",
           y = "Average Gestation Length (in days)") +
      theme(
        plot.title = element_text(hjust = 0.5))
  
      ###Checking assumptions###
      #Normality# - Not normal distribution
      hist(model8$residuals, main="Histogram of Residuals") #data name precedes $
      qqnorm(model8$residuals, main="QQ Plot") #normality plot residuals
      qqline(model8$residuals) #want all or majority points fall near line
      ad.test(model8$residuals) #Anderson-Darling normality test
      shapiro.test(model8$residuals) #Shapiro-Wilk normality test
      lillie.test(model8$residuals) #Kolmogorov-Smirnov normality test
      #Constant Residual Variance#
      plot(model8$fitted,model8$residuals, ylab="Residuals", xlab="Fitted Values")
      abline(0,0)
      bptest(model8) #Breusch-Pagan for constant variance
      #Linearity#
      abline(model8, col = "blue")
      cor.test(syn$adult_mass,syn$litters_annually)
      #Independence#
      plot(model8$residuals, main="Residuals vs Order")
      acf(model8$residuals)
      dwtest(model8)
    
    ###Adult Mass v Gestation Length, Synanthropic###
    model9 = lm(gestation~adult_mass, data=syn) #ns#
      summary(model9)
      anova(model9) 
      
      ###Plot###
      ggplot(syn, aes(x = adult_mass, y = gestation)) +
        geom_point(color = "black", size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
        labs(title = "Average Adult Mass vs. Average Gestation Length, Synanthropic Rodents",
             x = "Average Adult Mass (in grams)",
             y = "Average Gestation Length (in days)") +
        theme(
          plot.title = element_text(hjust = 0.5))
      
    ###Adult Mass v Litter Size, Nonsynanthropic###
    model10 = lm(litter_size~adult_mass, data=nonsyn) #ns
      summary(model10)
      anova(model10)
      
      ##Plot###
      ggplot(nonsyn, aes(x = adult_mass, y = litter_size)) +
        geom_point(color = "black", size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
        labs(title = "Average Adult Mass vs. Average Litter Size, Non-synanthropic Rodents",
             x = "Average Adult Mass (in grams)",
             y = "Average Litter Size (number of offspring)") +
        theme(
          plot.title = element_text(hjust = 0.5))
      
    ###Adult Mass v Litter Size, Synanthropic###
    model11 = lm(litter_size~adult_mass, data=syn) #ns
      summary(model11)
      anova(model11)
      
      ##Plot###
      ggplot(syn, aes(x = adult_mass, y = litter_size)) +
        geom_point(color = "black", size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
        labs(title = "Average Adult Mass vs. Average Litter Size, Synanthropic Rodents",
             x = "Average Adult Mass (in grams)",
             y = "Average Litter Size (number of offspring)") +
        theme(
          plot.title = element_text(hjust = 0.5))
      
    ###Gestation Length v Litter Size, Nonsynanthropic###
    litter_size_log_nonsyn <- log(nonsyn$litter_size)
    modelx8 = lm(litter_size_log_nonsyn~gestation_log_nonsyn, data=nonsyn) #Significant#
      summary(modelx8)
      anova(modelx8)
      
      ###Plot###
      ggplot(nonsyn, aes(x = exp(gestation_log_nonsyn), y = exp(litter_size_log_nonsyn))) +
        geom_point(color = "black", size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
        labs(title = "Average Gestation Length vs. Average Litter Size, Non-synanthropic Rodents",
             x = "Average Gestation Length (in days)",
             y = "Average Litter Size (number of offspring)") +
        theme(
          plot.title = element_text(hjust = 0.5))
      
    ###Gestation Length v Litter Size, Synanthropic###
    litter_size_log_syn <- log(syn$litter_size)
    gestation_log_syn <- log(syn$gestation)
    modelx9 = lm(litter_size_log_syn~gestation_log_syn, data=syn) #Significant#
      summary(modelx9)
      anova(modelx9) #Significant#
      
      ###Plot###
      ggplot(syn, aes(x = exp(gestation_log_syn), y = exp(litter_size_log_syn))) +
        geom_point(color = "black", size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
        labs(title = "Average Gestation Length vs. Average Litter Size, Synanthropic Rodents",
             x = "Average Gestation Length (in days)",
             y = "Average Litter Size (number of offspring)") +
        theme(
          plot.title = element_text(hjust = 0.5))
      
      ###Comparing slopes of the line###
        slope_nonsyn1<- coef(summary(modelx8))
        slope_nonsyn1
        slope_syn2<- coef(summary(modelx9))
        slope_syn2
        t_stat <- ((-4.51 - -4.93)/ sqrt(0.18^2 + 0.20^2))
        df <- ((0.18^2 + 0.20^2)^2) / (((0.18^2)^2 / (33 - 2)) + ((0.20^2)^2 / (24 - 2)))
        p_value <- 2 * pt(-abs(t_stat), df)
        t_stat
        df
        p_value
    
    ###Gestation Length v Litters Annually, Nonsynanthropic###
    modelx6 = lm(litters_annually~gestation, data=nonsyn) #ns#
      summary(modelx6)
      anova(modelx6)
      
      ###Plot###
      ggplot(nonsyn, aes(x = gestation, y = litters_annually)) +
        geom_point(color = "black", size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
        labs(title = "Average Gestation Length vs. Average Litters Annually, Non-synanthropic Rodents",
             x = "Average Gestation Length (in days)",
             y = "Average Litters Annually") +
        theme(
          plot.title = element_text(hjust = 0.5))
      
    ###Gestation Length v Litters Annually, Synanthropic### 
    modelx7 = lm(litters_annually~gestation, data=syn) #ns#
      summary(modelx7)
      anova(modelx7)
      
      ggplot(syn, aes(x = gestation, y = litters_annually)) +
        geom_point(color = "black", size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
        labs(title = "Average Gestation Length vs. Average Litters Annually, Synanthropic Rodents",
             x = "Average Gestation Length (in days)",
             y = "Average Litters Annually") +
        theme(
          plot.title = element_text(hjust = 0.5))
  
    ###Relative Brain v Body, Nonsynanthropic###
    model12 = lm(brain~body, data=nonsyn) #Significant
      summary(model12)
      anova(model12)
      
      ###Plot###
      ggplot(nonsyn, aes(x = body, y = brain)) +
        geom_point(color = "black", size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
        labs(title = "Relative Brain vs. Body Mass, Non-synanthropic Rodents",
             x = "Relative Body Mass",
             y = "Relative Brain Mass") +
        theme(
          plot.title = element_text(hjust = 0.5))
      
      ###Checking assumptions###
      #Normality# - Not normal distribution
      hist(model12$residuals, main="Histogram of Residuals") #data name precedes $
      qqnorm(model12$residuals, main="QQ Plot") #normality plot residuals
      qqline(model12$residuals) #want all or majority points fall near line
      ad.test(model12$residuals) #Anderson-Darling normality test
      shapiro.test(model12$residuals) #Shapiro-Wilk normality test
      lillie.test(model12$residuals) #Kolmogorov-Smirnov normality test
      #Constant Residual Variance#
      plot(model12$fitted,model12$residuals, ylab="Residuals", xlab="Fitted Values")
      abline(0,0)
      bptest(model12) #Breusch-Pagan for constant variance
      #Linearity#
      abline(model12, col = "blue")
      #Independence#
      plot(model12$residuals, main="Residuals vs Order")
      acf(model12$residuals)
      dwtest(model12)
  
  ###Relative Brain v Body, Synanthropic### 
  model13 = lm(brain~body, data=syn) 
    summary(model13)
    anova(model13) #Significant#
    
    ###Plot###
    ggplot(syn, aes(x = body, y = brain)) +
      geom_point(color = "black", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill="#69b3a2") +
      labs(title = "Relative Brain vs. Body Mass, Synanthropic Rodents",
           x = "Relative Body Mass",
           y = "Relative Brain Mass") +
      theme(
        plot.title = element_text(hjust = 0.5))
    
    ###Checking assumptions###
    #Normality# - Not normal distribution
    hist(model13$residuals, main="Histogram of Residuals") #data name precedes $
    qqnorm(model13$residuals, main="QQ Plot") #normality plot residuals
    qqline(model13$residuals) #want all or majority points fall near line
    ad.test(model13$residuals) #Anderson-Darling normality test
    shapiro.test(model13$residuals) #Shapiro-Wilk normality test
    lillie.test(model13$residuals) #Kolmogorov-Smirnov normality test
    #Constant Residual Variance#
    plot(model13$fitted,model13$residuals, ylab="Residuals", xlab="Fitted Values")
    abline(0,0)
    bptest(model13) #Breusch-Pagan for constant variance
    #Linearity#
    abline(model13, col = "blue")
    #Independence#
    plot(model13$residuals, main="Residuals vs Order")
    acf(model13$residuals)
    dwtest(model13)
    
      ###Comparing slopes of the line###
        slope_nonsyn <- coef(summary(model12))
          slope_nonsyn
        slope_syn <- coef(summary(model13))
          slope_syn
        t_stat <- ((0.84 - 0.64)/(sqrt((0.09^2+0.04^2))))
        df <- ((0.09^2 + 0.04^2)^2) / (((0.09^2)^2 / (33 - 2)) + ((0.04^2)^2 / (24 - 2)))
        p_value <- 2 * pt(-abs(t_stat), df)
        t_stat
        df
        p_value
        
    ###MLR - Litter Size vs Adult Mass + Gestation, Nonsynanthropic###
    litter_size_log_nonsyn <- log(nonsyn$litter_size)
    nonsynmlmodel = lm(litter_size_log_nonsyn~gestation_log_nonsyn+adult_mass_log_nonsyn, data=nonsyn) #Significant, but only for gestation#
      summary(nonsynmlmodel)
      anova(nonsynmlmodel)
        
      avPlots(nonsynmlmodel)
        
      ###Checking assumptions###
      #Normality#
      hist(nonsynmlmodel$residuals, main="Histogram of Residuals") #data name precedes $
      qqnorm(nonsynmlmodel$residuals, main="QQ Plot") #normality plot residuals
      qqline(nonsynmlmodel$residuals) #want all or majority points fall near line
      ad.test(nonsynmlmodel$residuals) #Anderson-Darling normality test
      shapiro.test(nonsynmlmodel$residuals) #Shapiro-Wilk normality test
      lillie.test(nonsynmlmodel$residuals) #Kolmogorov-Smirnov normality test
      #Constant Residual Variance#
      plot(nonsynmlmodel$fitted,nonsynmlmodel$residuals, ylab="Residuals", xlab="Fitted Values")
      abline(0,0)
      bptest(nonsynmlmodel) #Breusch-Pagan for constant variance
      #Linearity#
      abline(nonsynmlmodel, col = "blue")
      #Independence#
      plot(nonsynmlmodel$residuals, main="Residuals vs Order")
      acf(nonsynmlmodel$residuals)
      dwtest(nonsynmlmodel)
      
    ###MLR - Litter Size vs Adult Mass + Gestation, Synanthropic###
    litter_size_log_syn <- log(syn$litter_size)
    synmlmodel = lm(litter_size_log_syn~gestation_log_syn+adult_mass_log_syn+gestation_log_syn, data=syn) #Significant, but only for gestation#
      summary(synmlmodel)
      anova(synmlmodel)
      
      avPlots(synmlmodel)
      
      ###Checking assumptions###
      #Normality#
      hist(synmlmodel$residuals, main="Histogram of Residuals") #data name precedes $
      qqnorm(synmlmodel$residuals, main="QQ Plot") #normality plot residuals
      qqline(synmlmodel$residuals) #want all or majority points fall near line
      ad.test(synmlmodel$residuals) #Anderson-Darling normality test
      shapiro.test(synmlmodel$residuals) #Shapiro-Wilk normality test
      lillie.test(synmlmodel$residuals) #Kolmogorov-Smirnov normality test
      #Constant Residual Variance#
      plot(synmlmodel$fitted,synmlmodel$residuals, ylab="Residuals", xlab="Fitted Values")
      abline(0,0)
      bptest(synmlmodel) #Breusch-Pagan for constant variance
      #Linearity#
      abline(synmlmodel, col = "blue")
      #Independence#
      plot(synmlmodel$residuals, main="Residuals vs Order")
      acf(synmlmodel$residuals)
      dwtest(synmlmodel8)
