library(tidyr)
library(dplyr)
library(janitor)
library(ggplot2)
library(boot)
library(boot.pval)


#Analysis based on Meyer et al. 2019 https://www.frontiersin.org/articles/10.3389/fneur.2019.00556/full

########################
#step 1: tidy the data
########################

#extract lobe names order from first header
lobes <- as.vector(names(study1_mld_mibi_regression))
lobes <- lobes[!startsWith(lobes,'X')]

#define a data cleaning function
prana_clean <- function(dat, study_num=1) {
  #remove lobe names header
  names(dat) <- NULL
  
  #set 'subject', 'frame, 'GM', 'WM', etc. as the new header
  colnames(dat) <- dat[1,]
  dat <- dat[-1,]
  
  #make the header names unique and syntactically valid
  colnames(dat) <- make.names(colnames(dat), unique = TRUE)
  
  #remove blank columns
  dat <- remove_empty(dat, which = "cols")
  
  #add lobe names to the header names
  dat <- rename(dat, 
                GM.R.Frontal = GM,
                WM.R.Frontal = WM,
                CSF.R.Frontal = CSF,
                R2.R.Frontal = R2)
  
  colnames(dat) <- colnames(dat) %>% 
    stringr::str_replace_all(c("1" = "L.Frontal", 
                               "\\.2" = "\\.R.Temporal", # \\. prevents replacement in R2
                               "3" = "L.Temporal",
                               "4" = "R.Parietal",
                               "5" = "L.Parietal",
                               "6" = "R.Occipital",
                               "7" = "L.Occipital"))
  
  
  #replace empty subject names with "Average"
  dat$Subject[dat$Subject == ""] <- "Average"
  
  #replace empty cells with NA
  dat[dat == ""] <- NA
  
  #replace cells equal to <Min# with NA
  dat[dat == "<Min#"] <- NA
  
  #remove columns where everything is equal or NA
  dat <- remove_constant(dat, na.rm = TRUE)
  
  #remove rows where frame = NA
  dat<- dat[!is.na(dat$Frame), ]
  
  #add study number to dataset
  dat$Study.Number <- study_num
  
  #have R guess the class of each variable again, now that the data is clean
  dat <- type.convert(dat, as.is = TRUE)
  
  #pivot the data into wide format
  dat <- pivot_wider(dat, id_cols = c(Subject, Study.Number), names_from = Frame, values_from = colnames(dat)[3:26])
}

#import datasets for each study
study1_mld_mibi_regression <- read.delim("~/OneDrive - nyu.edu/Kirov Rotation Project/study1_mld_mibi_regression.txt")
study2_mld_mibi_regression <- read.delim("~/OneDrive - nyu.edu/Kirov Rotation Project/study2_mld_mibi_regression.txt")
study3_mld_mibi_regression <- read.delim("~/OneDrive - nyu.edu/Kirov Rotation Project/study3_mld_mibi_regression.txt")

#apply data cleaning function to each dataset
study1_mld_mibi_regression <- prana_clean(study1_mld_mibi_regression)
study2_mld_mibi_regression <- prana_clean(study2_mld_mibi_regression, study_num = 2)
study3_mld_mibi_regression <- prana_clean(study3_mld_mibi_regression, study_num = 3)

#combine all 3 studies into one dataframe
serial_mld_mibi_regression <- rbind(study1_mld_mibi_regression, study2_mld_mibi_regression, study3_mld_mibi_regression)

####
#average the R and L sides for each region
matter <- c('GM', 'WM')
measure <- c('Frontal_tCholine_Norm', 'Frontal_tCreatine_Norm', 'Frontal_Glx._Norm', 'Frontal_myo-Inositol_Norm', 'Frontal_NAcetylAspartate_Norm',
             'Temporal_tCholine_Norm', 'Temporal_tCreatine_Norm', 'Temporal_Glx._Norm', 'Temporal_myo-Inositol_Norm', 'Temporal_NAcetylAspartate_Norm',
             'Parietal_tCholine_Norm', 'Parietal_tCreatine_Norm', 'Parietal_Glx._Norm', 'Parietal_myo-Inositol_Norm', 'Parietal_NAcetylAspartate_Norm',
             'Occipital_tCholine_Norm', 'Occipital_tCreatine_Norm', 'Occipital_Glx._Norm', 'Occipital_myo-Inositol_Norm', 'Occipital_NAcetylAspartate_Norm')

#create a dataframe blank except for Subject and Study.Number to store new columns
serial_lobar_mld_mibi_regression <- serial_mld_mibi_regression %>% select(Subject, Study.Number)

for (mat in matter) {
  for (meas in measure) {
    #average the columns that start with matter and end with measure
    newcol <- rowMeans(serial_mld_mibi_regression %>% select(starts_with(mat) & ends_with(meas)))
    #append the averages to the end of the dataframe
    serial_lobar_mld_mibi_regression <- cbind(serial_lobar_mld_mibi_regression, newcol)
    #rename the new column matter_measure
    colnames(serial_lobar_mld_mibi_regression)[ncol(serial_lobar_mld_mibi_regression)] <- paste0(mat, "_", meas)
  }  
}

#replace dashes with underscores in the column names (myo-Inositol becomes myo_Inositol)
colnames(serial_lobar_mld_mibi_regression) = gsub("-", "_", colnames(serial_lobar_mld_mibi_regression))

###
#pull out the subject averages into separate dataframe
serial_averages_mld_mibi_regression <- serial_lobar_mld_mibi_regression %>% filter(Subject == 'Average')
serial_lobar_mld_mibi_regression <- serial_lobar_mld_mibi_regression %>% filter(Subject != 'Average')

#pull out subjects with 2 timepoints
serial_2_mld_mibi_regression <- serial_lobar_mld_mibi_regression %>%
  group_by(Subject) %>% filter(n() == 2)

#pull out subjects with all 3 timepoints
serial_3_mld_mibi_regression <- serial_lobar_mld_mibi_regression %>% 
  group_by(Subject) %>% filter(n() == 3)

###################
#Step 2: Linear Trends Over Time
###################

#pull the list of metabolites from the dataframe column names
metabolites <- serial_lobar_mld_mibi_regression %>% select(!c(Subject, Study.Number)) %>% colnames()
#create an empty list to store the looped lm models
storage <- list()
#create an empty data frame to store the looped pvalues
pvalues_linear <- data_frame(region_metabolite = NA, intercept_pvalue = NA, linear_pvalue = NA)

#generate lm model and calculate pvalue for each metabolite by Study.Number
for (i in metabolites) {
  storage[[i]] <- lm(get(i) ~ Study.Number, data = na.omit(serial_lobar_mld_mibi_regression))
  newpvalue <- c(i, summary(storage[[i]])$coefficients[,4])
  pvalues_linear <- rbind(pvalues_linear, newpvalue)
}

#get rid of intercept_pvalue column
pvalues_linear <- pvalues_linear %>% select(-intercept_pvalue)
#remove empty rows and sort the pvalues in ascending order
pvalues_linear <- pvalues_linear %>% remove_empty(which = 'rows') %>% arrange(linear_pvalue)

#######
#repeat for subset, subjects with 2 timepoints

storage_2 <- list()
pvalues_linear_2 <- data_frame(region_metabolite = NA, intercept_pvalue = NA, linear_pvalue = NA)

for (i in metabolites) {
  storage_2[[i]] <- lm(get(i) ~ Study.Number, data = na.omit(serial_2_mld_mibi_regression))
  newpvalue <- c(i, summary(storage_2[[i]])$coefficients[,4])
  pvalues_linear_2 <- rbind(pvalues_linear_2, newpvalue)
}

pvalues_linear_2 <- pvalues_linear_2 %>% select(-intercept_pvalue)
pvalues_linear_2 <- pvalues_linear_2 %>% remove_empty(which = 'rows') %>% arrange(linear_pvalue)

######
#repeat for subset, subjects with all 3 timepoints

storage_3 <- list()
pvalues_linear_3 <- data_frame(region_metabolite = NA, intercept_pvalue = NA, linear_pvalue = NA)

for (i in metabolites) {
  storage_3[[i]] <- lm(get(i) ~ Study.Number, data = na.omit(serial_3_mld_mibi_regression))
  newpvalue <- c(i, summary(storage_3[[i]])$coefficients[,4])
  pvalues_linear_3 <- rbind(pvalues_linear_3, newpvalue)
}

pvalues_linear_3 <- pvalues_linear_3 %>% select(-intercept_pvalue)
pvalues_linear_3 <- pvalues_linear_3 %>% remove_empty(which = 'rows') %>% arrange(linear_pvalue)


###########################
#step 3: bootstrap
###########################
#using storage and metabolites objects generated in step 2

#create dataframe of bootstrapped p values
pvalues_boot <- data_frame(region_metabolite = NA, bootstrap_pvalue = NA)
for (i in metabolites) {
  newpvalue <- c(i, boot_summary(storage[[i]])$p.value[2])
  pvalues_boot <- rbind(pvalues_boot, newpvalue)
}

#join linear pvalues and bootstrapped pvalues into one dataframe
pvalues_all <- inner_join(pvalues_linear, pvalues_boot)

#repeat for subset, subjects with 2 timepoints

pvalues_boot_2 <- data_frame(region_metabolite = NA, bootstrap_pvalue = NA)
for (i in metabolites) {
  newpvalue <- c(i, boot_summary(storage_2[[i]])$p.value[2])
  pvalues_boot_2 <- rbind(pvalues_boot_2, newpvalue)
}

pvalues_all_2 <- inner_join(pvalues_linear_2, pvalues_boot_2)


#repeat for subset, subjects with all 3 timepoints

pvalues_boot_3 <- data_frame(region_metabolite = NA, bootstrap_pvalue = NA)
for (i in metabolites) {
  newpvalue <- c(i, boot_summary(storage_3[[i]])$p.value[2])
  pvalues_boot_3 <- rbind(pvalues_boot_3, newpvalue)
}

pvalues_all_3 <- inner_join(pvalues_linear_3, pvalues_boot_3)

#####################
#step 4: plots
######################
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) +
    theme_bw()
}


WM_Frontal_tCreatine_Norm <- lm(WM_Frontal_tCreatine_Norm ~ Study.Number, data = na.omit(serial_2_mld_mibi_regression))
ggplotRegression(WM_Frontal_tCreatine_Norm)


WM_Parietal_Glx._Norm <- lm(WM_Parietal_Glx._Norm ~ Study.Number, data = na.omit(serial_2_mld_mibi_regression))
ggplotRegression(WM_Parietal_Glx._Norm)

WM_Frontal_Glx._Norm

WM_Frontal_Glx._Norm <- lm(WM_Frontal_Glx._Norm ~ Study.Number, data = na.omit(serial_2_mld_mibi_regression))
ggplotRegression(WM_Frontal_Glx._Norm)

WM_Occipital_tCholine_Norm <- lm(WM_Occipital_tCholine_Norm ~ Study.Number, data = na.omit(serial_2_mld_mibi_regression))
ggplotRegression(WM_Occipital_tCholine_Norm)

#r2 is stronger in recovered (transition from non-recovered to recovered, bigger decrease)
