###########################################################################
########################### EDA ###########################################
###########################################################################
library(survival)
library(xtable)
library(Hmisc)
source("http://www.ics.uci.edu/~dgillen/STAT211/Handouts/Stat211Functions.R")

###########################################################################
##### Missing Data
# NAs
# What is the structure?
describe(demo_outcome_data)
na.patterns <- naclus(demo_outcome_data)
pdf("./plots/na_patterns.pdf");
par(mar=c(3, 2, 2, 0) + 0.1); plot( na.patterns ) ; par(mar=c(5, 4, 4, 2) + 0.1);
dev.off();
# Ethnic, age and sex are missing together. Information lost for some people?

pdf("./plots/na_ethnic_regressed.pdf");
plot( summary(is.na(ethnic) ~ age + sex, data=demo_outcome_data))
dev.off()
# All that are missing in ethnic are also missing in age and sex
# But there are also Male, Female and individuals in higher two
# age brackets with missing ethnicity

###########################################################################
##### Table
summary(demo_outcome_data)

### Continuous
head(demo_outcome_data)
cont.names <- names(demo_outcome_data)[5]
cont.vars <- as.matrix(demo_outcome_data[,cont.names])
n <- nrow(cont.vars); p <- ncol(cont.vars)
cont.table <- matrix(0, p, 2)
for(i in 1:p){
  cont.table[i,2] <- sum(is.na(cont.vars[,i]))
  #cont.table[i, 1] <- cont.names[i]
  cont.table[i,1] <- paste(round(median(cont.vars[,i],na.rm=T),2), "(",round(sd(cont.vars[,i], na.rm=T),2),")", sep="")
}
cont.table <- as.matrix(cont.table)
rownames(cont.table) <- cont.names
colnames(cont.table) <- c("median (SD)", "Number of NAs")
xtable(cont.table)

### Categorical
df <- demo_outcome_data
cont.df <- subset(demo_outcome_data, demo_outcome_data$tx==1)
names(cont.df) <- c("Patient ID", "Site", "Treatment", "Aspirin use", "Age", "Sex", "Ethnicity", "Adenomas")
names(df) <- c("Patient ID", "Site", "Treatment", "Aspirin use", "Age", "Sex", "Ethnicity", "Adenomas")
# Categorical
cat.names <- names(cont.df)[c(2, 4, 6, 7, 8)]
cat.vars <- cont.df[,cat.names]
df.vars <- df[,cat.names]
n <- nrow(cat.vars); p <- ncol(cat.vars)
cat.table <- matrix(0, 50, 2)
cat.table.rownames <- NULL
levelcount <- 0
for(i in 1:p){
  var <- cat.vars[!is.na(cat.vars[,i]),i]
  if(class(var)=="factor")
  {
    var <- as.numeric(var)
  }
  no_na <- sum(is.na(cat.vars[,i]))
  if(max(var)>1){
    cat.table.rownames[i+levelcount] <- cat.names[i]
    cat.table[i+levelcount,1] <- "---"
    cat.table[i+levelcount,2] <- no_na
    for(j in 1:length(unique(df.vars[!is.na(cat.vars[,i]),i]))){
      q <- match(cat.names[i], names(cont.df))
      if(class(cont.df[,q])=="factor"){
        cat.table.rownames[i+levelcount+j] <- levels(df[,q])[j]
      }else{
        cat.table.rownames[i+levelcount+j] <- paste(cat.names[i], ": ", sort(unique(df[,q]))[j], sep="")
      }
      cat.table[i+levelcount+j,1] <- paste(round(sum(var==sort(unique(var))[j]),2), "(",round(mean(var==sort(unique(var))[j])*100,2),"%)", sep="")
      cat.table[i+levelcount+j,2] <- ""
    }
    levelcount <- levelcount + j
  }else{
    cat.table.rownames[i+levelcount] <- cat.names[i]
    cat.table[i+levelcount,1] <- paste(round(sum(var),2), "(",round(mean(var),2)*100,"%)", sep="")
    cat.table[i+levelcount,2] <- no_na
  }
}
cat.table <- as.matrix(cat.table)[1:(which(cat.table[,1]=="0")[1]-1),]
rownames(cat.table) <- cat.table.rownames
colnames(cat.table) <- c("Count (Proportion)", "NAs")
cat.table
xtable(cat.table)

cat.table.trt
cat.table.placebo 
cbind(cat.table.placebo, rbind(cat.table.trt, c("", "", ""), c("", "", "")))

###########################################################################
##### Bivariate

###########################################################################
##### Adenomas vs Trt
# First look at this - most interesting
table(demo_outcome_data$adenomas[which(demo_outcome_data$tx==1)])
table(demo_outcome_data$adenomas[which(demo_outcome_data$tx==0)])
# definitely very different, can we plot this?
counts <- table(demo_outcome_data$tx, demo_outcome_data$adenomas)
rownames(counts) <- c("Placebo", "Treatment")
xtable(counts)
pdf("./plots/bivariate_barplot_adenomas.pdf");
barplot(counts, main="Adenoma count by Treatment",
        xlab="Number of Adenomas", col=c("darkblue","red"),
        beside=TRUE)
legend(15, 150, c("Placebo", "Treatment"),
       fill = c("darkblue","red"), bty="n"); dev.off()
box(which = "plot", lty = "solid")

## Means across treatment groups
mean(as.numeric(as.matrix(subset(demo_outcome_data, tx==1, select="adenomas"))))
sd(as.numeric(as.matrix(subset(demo_outcome_data, tx==1, select="adenomas"))))
mean(as.numeric(as.matrix(subset(demo_outcome_data, tx==0, select="adenomas"))))
sd(as.numeric(as.matrix(subset(demo_outcome_data, tx==0, select="adenomas"))))


###########################################################################
##### Adenomas vs Trt by Aspirin use
trt.asa.inter <- as.factor(paste(demo_outcome_data$tx, demo_outcome_data$asa_use, sep=""))
inter.counts <- table(trt.asa.inter, demo_outcome_data$adenomas)/c(113,69,113,69)
pdf("./plots/trivariate_barplot_adenomas.pdf")
barplot(inter.counts, main="Adenoma count by Treatment and Aspirin",
        xlab="Number of Adenomas", col=c("darkred", "red", "darkblue", "blue"),
        beside=TRUE, ylab="Proportion")
legend(15, 0.8, c("Placebo, no Aspirin", "Placebo, Aspirin", "Treatment, no Aspirin", "Treatment, Aspirin"),
       fill = c("darkred", "red", "darkblue", "blue"), bty="n")
dev.off()

## Means across treatment*aspirin groups
mean(as.numeric(as.matrix(subset(demo_outcome_data, tx==0 & asa_use==0, select="adenomas"))))
sd(as.numeric(as.matrix(subset(demo_outcome_data, tx==0 & asa_use==0, select="adenomas"))))
mean(as.numeric(as.matrix(subset(demo_outcome_data, tx==0 & asa_use==1, select="adenomas"))))
sd(as.numeric(as.matrix(subset(demo_outcome_data, tx==0 & asa_use==1, select="adenomas"))))
mean(as.numeric(as.matrix(subset(demo_outcome_data, tx==1 & asa_use==0, select="adenomas"))))
sd(as.numeric(as.matrix(subset(demo_outcome_data, tx==1 & asa_use==0, select="adenomas"))))
mean(as.numeric(as.matrix(subset(demo_outcome_data, tx==1 & asa_use==1, select="adenomas"))))
sd(as.numeric(as.matrix(subset(demo_outcome_data, tx==1 & asa_use==1, select="adenomas"))))

###########################################################################
##### Longitudinal
head(audio_long)

## Plot 1: Panelled
p <- ggplot(data = audio_long, aes(x = time_passed, y = avg_db, group = study_no))
p + geom_line() + guides(colour=FALSE) + aes(colour = tx) + facet_wrap(~ tx)

## Plot 2: Mixed
pdf("./plots/longitudinal.pdf")
p <- ggplot(data = audio_long, aes(x = time_passed, y = avg_db, group = study_no))
p + geom_line() + guides(colour=FALSE, alpha=F) + aes(colour = tx)
dev.off()

## Plot 3: Boxplot stratified by treatment
pdf("./plots/longitudinal_box.pdf")
df <- data.frame(time_passed=cut(audio_long$time_passed, 3), avg_db=audio_long$avg_db, tx=audio_long$tx)
plt <- ggplot(data = df, aes(x=time_passed, y=avg_db))
plt + geom_boxplot(aes(fill = factor(tx)))
dev.off()

dim(audio_long)
## Plot 4: histogram of days since treatment stop
days_since_trt_stop <- ifelse(audio_long$Test_no==3, pmax(audio_long$time_passed-as.numeric(audio_long$days_on_tx),0), NA)
head(cbind(audio_long, days_since_trt_stop))
summary(days_since_trt_stop)
pdf("./plots/days_since_trt_stop.pdf", height=4, width=7)
par(mfrow=c(1,2))
hist(days_since_trt_stop[(audio_long$tx==0)], breaks=seq(0, 420, 20), col='darkred', xlab="Days from last treatment to hearing test", main="Placebo")
hist(days_since_trt_stop[(audio_long$tx==1)], breaks=seq(0, 420, 20), col='darkred', xlab="Days from last treatment to hearing test", main="Treatment")
par(mfrow=c(1,1))
dev.off()


## Mean of first / second / third test
db_table <- aggregate(audio_long$avg_db, list(Trt=audio_long$tx, Test_number=audio_long$Test_no), FUN=mean)
xtable(db_table)

## Correlation
fit.audio <- lm(avg_db ~ tx*time_passed + age + sex, data=audio_long)
summary(fit.audio)
R <- resid(fit.audio)
nR <- length(R)
R.minus.one <- R[1:(nR-1)]
R <- R[2:nR]
lm(R~0+R.minus.one)
# Or estimate manually
R%*%R.minus.one / R%*%R
sqrt(sum((R-R.minus.one)^2)/(length(R.minus.one)-1))
# This should be valid for estimating the overall covariance between observations
# But this is doing it all across subjects - need to do it per subject and then
# average