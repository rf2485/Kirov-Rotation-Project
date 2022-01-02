#this function takes raw linear regression text files output by PRANA and puts them into a format workable in R
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
