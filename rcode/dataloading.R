###########################################################################
############### FYE: Lars Hertel ##########################################
###########################################################################

demo_outcome_data <- dget("http://www.ics.uci.edu/~dgillen/FYQualExam2016/demo_outcome_data")
audio_long <- dget("http://www.ics.uci.edu/~dgillen/FYQualExam2016/audio_long")
dcarb <- demo_outcome_data


## Demographic data
summary(demo_outcome_data)
levels(demo_outcome_data$ethnic)
# Recode sex variable
#demo_outcome_data$male <- ifelse(demo_outcome_data$sex=="Male", 1, 0) 
#demo_outcome_data$male[which(demo_outcome_data$sex=="")] <- NA
# Rename sex=="" as NA
demo_outcome_data$sex[which(demo_outcome_data$sex=="")] <- NA
demo_outcome_data$sex <- factor(demo_outcome_data$sex)
# Rest is okay!

###########################################################################
##### Binarize African American
demo_outcome_data$AfrAmeric <- ifelse(demo_outcome_data$ethnic=="Black", 1, 0)
sum(demo_outcome_data$AfrAmeric, na.rm=T) # Check, correct # of African Americans
names(demo_outcome_data)

## Audio data
head(audio_long)
# NA's?
sum(is.na(audio_long)) #...yes
# assume it's the missing hearing recordings as in description
# hope that mean will be robust to that, move on for now
#audio_long$subject <- study_no
audio_long$study_no <- as.numeric(audio_long$study_no)

# Make avg hearing threshold column
audio <- apply(as.matrix(audio_long[4:19]), 2, as.numeric)
audio_long$avg_db <- apply(audio, 1, mean, na.rm=T)

# Make time difference column
time_passed <- NULL
current_id <- 0
demographics <- demo_outcome_data[1, 2:9]
for(i in 1:nrow(audio_long)){
  if(as.numeric(audio_long[i, 1])!=current_id){
    patient_start_date <- as.Date(audio_long[i, 3])
    time_passed[i] <- 0
    current_id <- as.numeric(audio_long[i, 1])
    demographics <- rbind(demographics, demo_outcome_data[which(demo_outcome_data$study_no==current_id), 2:9])
  }else{
    time_passed[i] <- as.numeric(as.Date(audio_long[i, 3]) - patient_start_date)
    demographics <- rbind(demographics, demo_outcome_data[which(demo_outcome_data$study_no==current_id), 2:9])
  }
}
audio_long$time_passed <- time_passed
demographics <- demographics[2:nrow(demographics),]
audio_long <- cbind(audio_long, demographics)
audio_long <- cbind(audio_long, rep(c(1,2,3), 184))
names(audio_long)[30] <- "Test_no"
