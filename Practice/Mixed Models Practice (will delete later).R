######### Mixed Models 1#########

remove(list = ls())

library(lme4)

Estuaries <- read.csv("Estuaries.csv", header = T)

#Fit a model for any dependent variables with a continuous distribution
ft.estu <- lmer(Total ~ Modification + (1|Estuary),data = Estuaries, REML=T)
#Total= dependent, Modificaiton = fixed, Estuary = random effect

#Assumption of normality reasonable if in a straight line
qqnorm(residuals(ft.estu))

#Diff in Variability between estuaries but variability doesn't increase w/ mean
scatter.smooth(residuals(ft.estu)~fitted(ft.estu))

#Fit full model by max likelihood & second model that lacks fixed effect oof modification
ft.estu <- lmer(Total ~ Modification + (1|Estuary), data=Estuaries, REML=F)
ft.estu.0 <- lmer(Total ~ (1|Estuary), data=Estuaries, REML=F)

#Compare both models w/ likelihood ratio test
anova(ft.estu.0,ft.estu) #shows effect of Modification p = 0.02385

#Calculate confidence intervals for each model 
confint(ft.estu) #effect of Modification has 95% confide4nce intervals that doo not overlap zero

nBoot=1000
lrStat=rep(NA,nBoot)
ft.null <- lm(Total~Modification,data=Estuaries)#null model
ft.alt <- lmer(Total~Modification+(1|Estuary),data=Estuaries,REML=F)#alternate model
lrObs <- 2*logLik(ft.alt) - 2*logLik(ft.null) #observed test stat
for(iBoot in 1:nBoot)
{
  Estuaries$TotalSim=unlist(simulate(ft.null)) #resampled data
  bNull <- lm(TotalSim~Modification,data=Estuaries)#null model
  bAlt <- lmer(TotalSim~Modification+(1|Estuary),data=Estuaries,REML=F)#alternate model
  lrStat[iBoot] <- 2*logLik(bAlt) - 2*logLik(bNull)#resampled test stat
}
mean(lrStat>lrObs) #P-value for test of Estuary effect

#Graph
ModEst <- unique(Estuaries[c("Estuary", "Modification")]) #find which Estuaries are modified
cols <- as.numeric(ModEst[order(ModEst[,1]),2])+3 #Assign colour by modification
boxplot(Total~ Estuary,data=Estuaries,col=cols,xlab="Estuary",ylab="Total invertebrates")
legend("bottomleft", inset=.02,
       c(" Modified "," Pristine "), fill=unique(cols), horiz=TRUE, cex=0.8)

is.mod <- as.numeric(ModEst[order(ModEst[,1]),2])-1 #0 if Modified, 1 if Pristine
Est.means <- coef(ft.estu)$Estuary[,1]+coef(ft.estu)$Estuary[,2]*is.mod #Model means
stripchart(Est.means~ sort(unique(Estuary)),data=Estuaries,pch=18,col="red",vertical = TRUE,add=TRUE)

##########################################################################################################################################
######### Mixed Models 2 #########

#Look at data and cross tabulation
Estuaries[1:10,]
xtabs(~ Estuary + Site, Estuaries, sparse = TRUE)

#Convert Site to a factor and create a new variable (SiteWithin= combo of Estuary and Site) ==> Create unique label for each site in the data
Estuaries$Site <- as.factor(Estuaries$Site)

Estuaries$SiteWithin <- with(Estuaries, factor(Estuary:Site))

#Check to see that each site is nested in only one Estuary
xtabs(~ Estuary + SiteWithin, Estuaries, sparse = TRUE)

#Fit a modoel for total abundance
fit.mod <- lmer(Total ~ Modification + (1|Estuary) + (1|SiteWithin), data = Estuaries)
summary(fit.mod)

#Fit wrong model to see if there's a difference of labels
fit.wrong <- lmer(Total ~ Modification + (1|Estuary) + (1|Site), data = Estuaries)
summary(fit.wrong )

#Look for straight line relationship on normal quantile plot & look for a fan and U-shape on the residual vs. fitted plot
par(mfrow=c(1,2))
qqnorm(residuals(fit.mod))
scatter.smooth(residuals(fit.mod)~fitted(fit.mod))

#Transform for better response
fit.mod <- lmer(log(Total) ~ Modification + (1|Estuary) + (1|SiteWithin), data = Estuaries)
par(mfrow=c(1,2))
qqnorm(residuals(fit.mod))
scatter.smooth(residuals(fit.mod)~fitted(fit.mod)) #residual plot

#Obtain approximate p-values for fixed effects
ft.mod <- lmer(log(Total)~Modification + (1|Estuary) + (1|SiteWithin), data=Estuaries, REML=F)
ft.mod.0 <- lmer(log(Total)~(1|Estuary) + (1|SiteWithin), data=Estuaries, REML=F)

#Parametric boootstrap
nBoot=1000
lrStat <- rep(NA,nBoot)
ft.null <- lmer(log(Total) ~ Modification + (1|Estuary) , Estuaries, REML=F)#null model
ft.alt <- lmer(log(Total) ~ Modification + (1|Estuary) + (1|SiteWithin), Estuaries, REML=F)#alternate model
lrObs <- 2*logLik(ft.alt) - 2*logLik(ft.null) #observed test stat
for(iBoot in 1:nBoot)
{
  Estuaries$logTotalSim <- unlist(simulate(ft.null)) #resampled data
  bNull <- lmer(logTotalSim~Modification + (1|Estuary) , Estuaries, REML=F)#null model
  bAlt <- lmer(logTotalSim~Modification+(1|Estuary)+ (1|SiteWithin), Estuaries, REML=F)#alternate model
  lrStat[iBoot] <- 2*logLik(bAlt) - 2*logLik(bNull)#resampled test stat
}
mean(lrStat>lrObs) #P-value for test of Estuary effect
anova(ft.mod.0,ft.mod)

##########################################################################################################################################
######### Mixed Models 3 #########

#Create variable to model presence/absence
Estuaries$HydroidPres <- Estuaries$Hydroid > 0

#Fit model for presence/absence of hydroids
fit.bin <- glmer(HydroidPres ~ Modification + (1|Estuary), family=binomial, data=Estuaries)

#Examine residual plots to check assumptions
par(mfrow=c(1,2))
plot(residuals(fit.bin)~fitted(fit.bin),main="residuals v.s. Fitted")
qqnorm(residuals(fit.bin))

#Use parametric bootstrap 
nBoot <- 1000
lrStat <- rep(NA,nBoot)
ft.null <- glmer(HydroidPres ~ 1 + (1|Estuary) ,family=binomial, data=Estuaries) #null model
ft.alt <- glmer(HydroidPres ~ Modification + (1|Estuary) ,family=binomial, data=Estuaries) #alternate model

lrObs <- 2*logLik(ft.alt) - 2*logLik(ft.null) #observed test stat

for(iBoot in 1:nBoot)
{
  Estuaries$HydroidPresSim <- unlist(simulate(ft.null)) #resampled data
  tryCatch({#sometimes the glmer code doesn't converge
    
    bNull <- glmer(HydroidPresSim ~ 1 + (1|Estuary) ,family=binomial, data=Estuaries)#null model
    bAlt <- glmer(HydroidPresSim ~ Modification + (1|Estuary) ,family=binomial, data=Estuaries)#alternate model
    lrStat[iBoot] <- 2*logLik(bAlt) - 2*logLik(bNull) #resampled test stat
  },warning=function(war) {lrStat[iBoot]=NA},error=function(err){lrStat[iBoot]=NA}) #if code doesn't converge skip sim
}
mean(lrStat>lrObs,na.rm=T) #P-value for test of Estuary effect

#Model the counts of hydroids
fit.pois <- glmer(Hydroid ~ Modification + (1|Estuary) ,family=poisson, data=Estuaries)

#Check the assumptions
par(mfrow=c(1,2))
plot(residuals(fit.pois)~fitted(fit.pois),main="Residuals vs. Fitted")
qqnorm(residuals(fit.pois))

#Parametric bootstrap to test for an effect of Modification
nBoot <- 1000
lrStat <- rep(NA,nBoot)
ft.null <- glmer(Hydroid ~ 1 + (1|Estuary) ,family=poisson, data=Estuaries) #null model
ft.alt <- glmer(Hydroid ~ Modification + (1|Estuary) ,family=poisson, data=Estuaries) #alternate model

lrObs <- 2*logLik(ft.alt) - 2*logLik(ft.null) #observed test stat
for(iBoot in 1:nBoot)
{
  Estuaries$HydroidSim <- unlist(simulate(ft.null)) #resampled data
  tryCatch({
    bNull <- glmer(HydroidSim ~ 1 + (1|Estuary)  ,family=poisson, data=Estuaries)#null model
    bAlt <- glmer(HydroidSim ~ Modification + (1|Estuary) ,family=poisson, data=Estuaries) #alternate model
    lrStat[iBoot] <- 2*logLik(bAlt) - 2*logLik(bNull) #resampled test stat
  },warning=function(war) {lrStat[iBoot]=NA},error=function(err){lrStat[iBoot]=NA})  #if code doesn't converge skip sim#   lrStat[iBoot]
}
mean(lrStat>lrObs,na.rm=TRUE) #P-value for test of Estuary effect

#Model counts of bryozoan, Schizoporella errata
fit.pois2 <- glmer(Schizoporella.errata ~ Modification + (1|Estuary), family=poisson,  data=Estuaries)
par(mfrow=c(1,2))
plot(residuals(fit.pois)~fitted(fit.pois),main="residuals vs. Fitted")
qqnorm(residuals(fit.pois))

#Communicate results
fit.pois <- glmer(Hydroid ~ Modification + (1|Estuary) ,family=poisson, data=Estuaries)
means <- fitted(fit.pois) #this will give the estimate at each data point
ModEst <- unique(Estuaries[c("Estuary", "Modification")])#find which Estuaries are modified
cols <- as.numeric(ModEst[order(ModEst[,1]),2])+3 #Assign colour by modification
boxplot(Hydroid~ Estuary,data=Estuaries,col=cols,xlab="Estuary",ylab="Count of hydroids")
legend("topleft", inset=.02,
       c("Modified","Pristine"), fill=unique(cols), horiz=TRUE, cex=0.8)

Est.means <- summarize(means,Estuaries$Estuary, mean)$means #extract means by Estuary
stripchart(Est.means~ sort(unique(Estuary)),data=Estuaries,pch=18,col="red",vertical = TRUE,add=TRUE) #plot means by estuary
