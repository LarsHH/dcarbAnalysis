###########################################################################
########################### Diagnostics ###################################
###########################################################################

###########################################################################
##### Treatment
##### Pearson Residuals: Mean-variance relationship
pdf("./plots/model_i_pearson.pdf", height=4, width=7)
par(mfrow=c(1,3))
presids <- residuals( fit.main, type="pearson" )
muhat <- fitted( fit.main )
plot( muhat, presids^2, xlab="Fitted value", ylab="Pearson Residual Squared", main="(a)")
# There are two points that are extreme outliers

plot( muhat, presids^2, xlab="Fitted value", ylab="Pearson Residual Squared", ylim=c(0,5), main="(b)")
sfit <- loess(presids^2 ~ muhat )
lines( sfit$x[ order( sfit$x ) ], fitted(sfit)[ order( sfit$x ) ],, col="red", lwd=2 )
# Line is clearly pulled up by those outliers

plot( muhat[-c(55, 85)], presids[-c(55, 85)]^2, xlab="Fitted value", ylab="Pearson Residual Squared", ylim=c(0,5), main="(c)")
sfit <- loess(presids[-c(55, 85)]^2 ~ muhat[-c(55, 85)] )
lines( sfit$x[ order( sfit$x ) ], fitted(sfit)[ order( sfit$x ) ],, col="red", lwd=2 )
# Removing them, the line is still weird. Almost looks like a -mu^2 relationship in the variance
par(mfrow=c(1,1))
dev.off()

## Use Robust Variance Estimator
robust.se.glm(fit.main)


#####	Plot of deviance residuals
pdf("./plots/model_i_diagnostics.pdf", heigh=4, width=7)
par(mfrow=c(1,2))
dresid <- residuals( fit.main, type="deviance" )
plot( muhat, dresid, xlab="Fitted value", ylab="Deviance residual", main="(a)")
abline(h=c(-1.96, 1.96), col="darkred")

#####	Plot of df betas
dfbeta <- dfbeta(fit.main)
ids <- complete.cases(demo_outcome_data[,c("tx")])
i=1
plot( 1:length(dfbeta[,i]), dfbeta[,i], xlab="Patient ID",
      ylab="Delta Beta for Treatment", main="(b)" )
large.idxs <- match(sort(-abs(dfbeta[,i]))[1:3], -abs(dfbeta[,i]))
text( c(30,200,50), dfbeta[large.idxs,i], demo_outcome_data$study_no[ids][large.idxs], col="darkred" )
dev.off()
par(mfrow=c(1,1))



###########################################################################
##### Treatment and Aspirin
##### Pearson Residuals: Mean-variance relationship
pdf("./plots/model_ii_pearson.pdf", height=4, width=7)
par(mfrow=c(1,3))
presids <- residuals( fit.interaction, type="pearson" )
muhat <- fitted( fit.interaction )
plot( muhat, presids^2, xlab="Fitted value", ylab="Pearson Residual Squared", main="(a)")
# There are two points that are extreme outliers

plot( muhat, presids^2, xlab="Fitted value", ylab="Pearson Residual Squared", ylim=c(0,5), main="(b)")
sfit <- loess(presids^2 ~ muhat )
lines( sfit$x[ order( sfit$x ) ], fitted(sfit)[ order( sfit$x ) ],, col="red", lwd=2 )
# Line is clearly pulled up by those outliers

plot( muhat[-c(55, 85)], presids[-c(55, 85)]^2, xlab="Fitted value", ylab="Pearson Residual Squared", ylim=c(0,5), main="(c)")
sfit <- loess(presids[-c(55, 85)]^2 ~ muhat[-c(55, 85)] )
lines( sfit$x[ order( sfit$x ) ], fitted(sfit)[ order( sfit$x ) ],, col="red", lwd=2 )
# Removing them, the line is still weird. Almost looks like a -mu^2 relationship in the variance
par(mfrow=c(1,1))
dev.off()

## Use Robust Variance Estimator
robust.se.glm(fit.main)


#####	Plot of deviance residuals
pdf("./plots/model_ii_diagnostics.pdf", heigh=4, width=7)
par(mfrow=c(1,2))
dresid <- residuals( fit.interaction, type="deviance" )
plot( muhat, dresid, xlab="Fitted value", ylab="Deviance residual", main="(a)")
abline(h=c(-1.96, 1.96), col="darkred")

#####	Plot of df betas
dfbeta <- dfbeta(fit.interaction)
ids <- complete.cases(demo_outcome_data[,c("tx")])
i=1
plot( 1:length(dfbeta[,i]), dfbeta[,i], xlab="Patient ID",
      ylab="Delta Beta for Treatment", main="(b)" )
large.idxs <- match(sort(-abs(dfbeta[,i]))[1:3], -abs(dfbeta[,i]))
text( c(30,200,50), dfbeta[large.idxs,i], demo_outcome_data$study_no[ids][large.idxs], col="darkred" )
dev.off()
par(mfrow=c(1,1))


###########################################################################
##### Hearing effect
# Residuals vs fitted values
pdf("./plots/model_iii_resid.pdf")
plot(fit.audio.full)
dev.off()

pdf("./plots/model_iii_qq.pdf")
qqPlot(resid(fit.audio.full))
dev.off()

plot(fitted(fit.audio.full), resid(fit.audio.full)/2.960356487)

# ACF of residuals
acf(resid(fit.audio.full))
# lag 1 seems to be large
ACF(fit.audio.full)




