###########################################################################
####################### Model Fit and Inference ###########################
###########################################################################
library(lme4)
library(car)
library(nlme)



###########################################################################
##### Treatment
head(demo_outcome_data)
n <- dim(demo_outcome_data)[1]
fit.main <- glm(adenomas ~ tx + asa_use + age + sex + AfrAmeric,
                data=demo_outcome_data, family = poisson(link="log"))
summary(fit.main)
glmCI(fit.main)
## Note that due to overdispersion the variance estimates are off
## Also there is no clear form that we can adjust for. Using quasipoisson
## is not ideal because the line goes up and down. So phi is not constant.
## We could impose a structure with phi propto mu^2. But how do we do that?
## And do we really want to do that? Seems better to just use the robust
## variance estimator. Do we have enough observations?
results.main <- robust.se.glm(fit.main)
results.main.txt <- cbind(round(exp(results.main[,1]),3), round(exp(results.main[,3]),3), round(exp(results.main[,4]),3), ifelse(results.main[,6]<0.001, "<0.001", round(results.main[,6], 3)))
results.main.txt <- cbind(paste(results.main.txt[,1], " (", results.main.txt[,2], ", ", results.main.txt[,3], ")", sep=""), results.main.txt[,4])
rownames(results.main.txt) <- c("Intercept", "Treatment", "Aspirin", "Age (yrs)", "Male", "African American")
colnames(results.main.txt) <- c("Incidence Rate Ratio (95% CI)", "p-value")
xtable(results.main.txt)

###########################################################################
##### Treatment and Aspirin Interaction
fit.interaction <- glm(adenomas ~ tx*asa_use + age + sex + AfrAmeric,
             data=demo_outcome_data, family = poisson(link="log"))
summary(fit.interaction)
glmCI(fit.interaction)
## Same problem with the overdispersion (see Pearson residual plots). Use
## robust variance estimator again.
robust.se.glm(fit.interaction)
exp(robust.se.glm(fit.interaction)[,1])
point.estimates <- robust.se.glm(fit.interaction)[,1]
L <- rbind(c(1,0,0,0,0,0,0), c(0,0,0,1,0,0,0), c(0,0,0,0,1,0,0), c(0,0,0,0,0,1,0), c(0,0,1,0,0,0,0), c(0,1,0,0,0,0,0), c(0,1,1,0,0,0,1))
covs <- c("Intercept", "Age", "Male", "African American", "Aspirin, no D-Carb", "no Aspirin, D-Carb", "Aspirin and D-Carb")
CI.lo <- L%*%point.estimates - qnorm(0.975) * sqrt(diag(L%*%robust.vcov.glm(fit.interaction)%*%t(L)))
CI.hi <- L%*%point.estimates + qnorm(0.975) * sqrt(diag(L%*%robust.vcov.glm(fit.interaction)%*%t(L)))
z.scores <- L%*%point.estimates/sqrt(diag(L%*%robust.vcov.glm(fit.interaction)%*%t(L)))
p.vals <- 2*pnorm(abs(z.scores), lower.tail = F)

results.main.txt <- cbind(round(exp(L%*%point.estimates),3), round(exp(CI.lo),3), round(exp(CI.hi),3), ifelse(p.vals<0.001, "<0.001", round(p.vals, 3)))
results.main.txt <- cbind(paste(results.main.txt[,1], " (", results.main.txt[,2], ", ", results.main.txt[,3], ")", sep=""), results.main.txt[,4])
rownames(results.main.txt) <- c("Intercept", "Age", "Male", "African American", "Aspirin, no D-Carb", "no Aspirin, D-Carb", "Aspirin and D-Carb")
colnames(results.main.txt) <- c("Incidence Rate Ratio (95% CI)", "p-value")
xtable(results.main.txt)

lrtest(fit.main, fit.interaction)

### Adding all ethnicities
fit.interaction.ethnic <- glm(adenomas ~ tx*asa_use + age + sex+ ethnic,
                       data=demo_outcome_data, family = poisson(link="log"))
robust.se.glm(fit.interaction.ethnic)


###########################################################################
##### Treatment Hearing Loss
head(audio_long)
fit.audio <- lm(avg_db ~ tx*time_passed + age + sex, data=audio_long)
summary(fit.audio)

####### Reduced Model
fit.audio.reduced <- lme(avg_db ~ 1 + time_passed + age + sex, data=audio_long,
                        random= ~(1+time_passed)|study_no, na.action=na.exclude)
summary(fit.audio.reduced)

audio_long$days_on_tx <- as.numeric(audio_long$days_on_tx)

####### Main Model
fit.audio.full <- lme(avg_db ~ 1 + time_passed*tx + age + sex, data=audio_long,
                      random= ~(1+time_passed)|study_no, na.action=na.exclude)
res.long <- summary(fit.audio.full)

names(summary(fit.audio.full))


L <- rbind(c(1,0,0,0,0,0), c(0,0,0,1,0,0), c(0,0,0,0,1,0), c(0,300,0,0,0,0), c(0,600,0,0,0,0), c(0,900,0,0,0,0), c(0,300,0,0,0,300), c(0,600,0,0,0,600), c(0,1,0,0,0,1))
point.estimate <- L %*% res.long$tTable[,1]
CI.lo <- point.estimate - qnorm(.975) * sqrt(diag(L%*%res.long$varFix%*%t(L)))
CI.hi <- point.estimate + qnorm(.975) * sqrt(diag(L%*%res.long$varFix%*%t(L)))
z.scores <- point.estimate / sqrt(diag(L%*%res.long$varFix%*%t(L)))
p.vals <- 2*pnorm(abs(z.scores), lower.tail = F)
p.vals <- round(res.long$tTable[,5],3)[c(1, 4, 5, )]

results.main.txt <- cbind(round(point.estimate,3), round(CI.lo,3), round(CI.hi,3), ifelse(p.vals<0.001, "<0.001", round(p.vals, 3)))
results.main.txt <- cbind(paste(results.main.txt[,1], " (", results.main.txt[,2], ", ", results.main.txt[,3], ")", sep=""), results.main.txt[,4])
rownames(results.main.txt) <- c("Intercept", "Age", "Male", "300 days", "600 days", "900 days", "D-Carb, 300 days", "D-Carb, 600 days", "D-Carb, 900 days")
colnames(results.main.txt) <- c("Estimated Change in Mean (95% CI)", "p-value")
xtable(results.main.txt)



####### Outliers Removed
fit.audio.full <- lme(avg_db ~ 1 + time_passed*tx + age + sex, data=audio_long[-c(388,389,390),],
                        random= ~(1+time_passed)|study_no, na.action=na.exclude)
summary(fit.audio.full)

####### LRT
LRT <- -2*(fit.audio.reduced$logLik - fit.audio.full$logLik)
1 - pchisq(LRT, 2)
anova(fit.audio.reduced, fit.audio.full)

