##### GLMM Analysis for individaul data ####
### Multilevel nested random effect model  Household Id nested with case control )
# required library
library(ggplot2)
library(tidyr)
library(dplyr)
library(MASS)
library(sjstats)
library(lme4)
library(glmulti)
library(broom)
## importing  data 
dd=read.csv("dd_cc_individual_data.csv",header=T,na.strings="")
### to see variable names 
colnames(dd)
#creating a case/control ID: 1/0
dd$case=NA
dd$case[dd$casexp=="Case"]=1
dd$case[dd$casexp=="Control"]=0
#creating a malaria positive ID: 1/0
dd$mal=NA
dd$mal[dd$malpos=="Positive"]=1
dd$mal[dd$malpos=="Negative"]=0
table(dd$case)
table(dd$mal)
#removing missing malpos data
dd2=dd[-which(dd$malpos=="NA"),]
dd$SITE=NA
dd$SITE[dd$site=="DD City"]=1
dd$SITE[dd$site=="DDU"]=0
### re-coding age category
dd$agecat=recode(dd$agecat,'Under 5'="0","5 - 15 Years"="1","Above 15 Years"="2")
dd$agecat <- factor(dd$agecat,
                    levels = c(0,1, 2),
                    labels = c("Under 5", "5 - 15 Years", "Above 15 Years"))

#GLMM Analysis 
##########################################################################
## univariate Analysis###
#########################################################################
## Stephensi presence
steph=glmer(mal~(1|idhh/case)+as.factor(stephpos),data=dd,family=binomial,nAGQ = 0)
summary(steph)
exp(1.3841)
table(dd$spray)
## confidence interval
se <- sqrt(diag(vcov(steph)))
tab <- cbind(Est = fixef(steph), LL = fixef(steph) - 1.96 * se, UL = fixef(steph) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## larvae presence
larv=glmer(mal~(1|idhh/case/SITE)+as.factor(larvaepos),data=dd,family=binomial,nAGQ = 0)
summary(larv)
se <- sqrt(diag(vcov(larv)))
tab <- cbind(Est = fixef(larv), LL = fixef(larv) - 1.96 * se, UL = fixef(larv) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
s=glmer(mal~(1|idhh/case)+site,data=dd,family=binomial,nAGQ = 1)
summary(s)
se <- sqrt(diag(vcov(s)))
tab <- cbind(Est = fixef(s), LL = fixef(s) - 1.96 * se, UL = fixef(s) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
s_sex=glmer(mal~(1|idhh/case)+site*sex,data=dd,family=binomial,nAGQ = 0)
summary(s_sex)
se <- sqrt(diag(vcov(s_sex)))
tab <- cbind(Est = fixef(s_sex), LL = fixef(s_sex) - 1.96 * se, UL = fixef(s_sex) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
s_age=glmer(mal~(1|idhh/case)+site*agecat,data=dd,family=binomial,nAGQ = 0)
summary(s_age)
se <- sqrt(diag(vcov(s_age)))
tab <- cbind(Est = fixef(s_age), LL = fixef(s_age) - 1.96 * se, UL = fixef(s_age) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)


## sex
sex=glmer(mal~(1|idhh/case/SITE)+as.factor(sex),data=dd,family=binomial,nAGQ = 0)
summary(sex)
se <- sqrt(diag(vcov(sex)))
tab <- cbind(Est = fixef(sex), LL = fixef(sex) - 1.96 * se, UL = fixef(sex) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## Age 
age= glmer(mal~(1|idhh/case/SITE) +agecat,data=dd,family=binomial,nAGQ = 0)
summary(age)
se <- sqrt(diag(vcov(age)))
tab <- cbind(Est = fixef(age), LL = fixef(age) - 1.96 * se, UL = fixef(age) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## adult presence
adult= glmer(mal~(1|idhh/case/SITE) +adultpos_2,data=dd,family=binomial,nAGQ = 0)
summary(adult)
se <- sqrt(diag(vcov(adult)))
tab <- cbind(Est = fixef(adult), LL = fixef(adult) - 1.96 * se, UL = fixef(adult) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## water body presence in the neighborhood
water= glmer(mal~(1|idhh/case/SITE) +as.factor(waterbody),data=dd,family=binomial,nAGQ = 0)
summary(water)
se <- sqrt(diag(vcov(water)))
tab <- cbind(Est = fixef(water), LL = fixef(water) - 1.96 * se, UL = fixef(water) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## use of spray 
spr= glmer(mal~(1|idhh/case/SITE) +as.factor(spray),data=dd,family=binomial,nAGQ = 0)
summary(spr)
se <- sqrt(diag(vcov(spr)))
tab <- cbind(Est = fixef(spr), LL = fixef(spr) - 1.96 * se, UL = fixef(spr) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## Multivariate Analysis 
ind_hh=glmer(mal~(1|idhh/case) +site+ as.factor(stephpos)+site*sex+site*agecat+as.factor(waterbody)+as.factor(spray),data=dd,family=binomial,nAGQ = 1)
summary(ind_hh)
performance::icc(ind_hh)
se <- sqrt(diag(vcov(ind_hh)))
tab <- cbind(OR = fixef(ind_hh), LL = fixef(ind_hh) - 1.96 * se, UL = fixef(ind_hh) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
All1=glmer(mal~(1|idhh/case/SITE)+stephpos+sex+relevel(as.factor(agecat),ref= "Under 5")+waterbody+scale(dist_river)+scale(dist_artificial_cont)+spray,data=dd,family=binomial)
summary(All1)
#### Household level analysis

hh=glmer(mal~(1|idhh/case/SITE)+factor(stephpos)+factor(waterbody)+factor(spray)+scale(dist_artificial_cont),data=dd,family=binomial,nAGQ = 0)
summary(hh)

se <- sqrt(diag(vcov(hh)))
tab <- cbind(Est = fixef(hh), LL = fixef(hh) - 1.96 * se, UL = fixef(hh) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
h=glmer(mal~(1|idhh)+stephpos+waterbody+as.factor(sp)+scale(dist_river)+scale(dist_artificial_cont),data=dd2,family=binomial)
summary(h)
exp(0.4744 )
A=glmer(mal~(1|idhh/case/SITE)+scale(dist_artificial_cont),data=dd,family=binomial,nAGQ = 0)
summary(A)
se <- sqrt(diag(vcov(A)))
tab <- cbind(Est = fixef(A), LL = fixef(A) - 1.96 * se, UL = fixef(A) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
d=glmer(mal~(1|idhh/case/SITE)+dist_river,data=dd,family=binomial,nAGQ = 0)
summary(d)
se <- sqrt(diag(vcov(d)))
tab <- cbind(Est = fixef(d), LL = fixef(d) - 1.96 * se, UL = fixef(d) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
DR=glmer(mal~(1|idhh:SITE)+scale(dist_river),data=dd,family=binomial,nAGQ = 0)
summary(DR)

HH_distance=glmer(mal~(1|idhh:case:SITE)+scale(dist_river)+scale(dist_artificial_cont),data=dd2,family=binomial,nAGQ = 0)
summary(HH_distance)

a=(glmer(mal ~ (1 | idhh/case/SITE) + stephpos + sex + relevel(as.factor(agecat),  ref = "Under 5") + waterbody + as.factor(spray) + scale(dist_river)+scale(dist_artificial_cont),data=dd,family=binomial,nAGQ=0))
summary(a)
### forest plot 
library(devtools)
devtools::install_github("strengejacke/sjPlot")
require(sjPlot)
set_theme(base = theme_classic())
plot_model(hh,title = "Forest plot HH level factors ")
plot_model(All)
plot_model(hh,show.values = T,title = "Forest plot HH level factors ")
plot_model(ind_hh,rm.terms=c("as.factor(waterbody)Yes","as.factor(agecat)Above 15 Years"),show.values = T,title = "Forest plot the full model")
plot_model(ind_hh,rm.terms=c("as.factor(waterbody)Yes","as.factor(agecat)Above 15 Years"),title = "Forest plot the full model")

## performance of the model 
performance::r2(ind_hh)
###############################################################################################
#Extracting the effect of random effects (so the background OR not explained by the fixed effects)
randoms=ranef(hh,condVar=TRUE)[[1]]
variances <- as.numeric(attr(randoms, "postVar"))
res=data.frame(WES = rownames(randoms), mean_effect = randoms$`(Intercept)`+sum(coef(summary(hh))[,"Estimate"]))
res$lower <- res$mean_effect - 2* sqrt(variances)
res$upper <- res$mean_effect + 2* sqrt(variances)
res$mean_effect <- exp(res$mean_effect)
res$lower <- exp(res$lower)
res$upper <- exp(res$upper)
res$WES <- reorder(res$WES, res$mean_effect, mean)
require(ggplot2)
ggplot(data=res,aes(x=WES,y=mean_effect))+geom_point(col="blue")+geom_errorbar(width=1,aes(ymin=lower,ymax=upper),
                                                                               col="black")+labs(title="Malaria Positiviy",y="Odds Ratio",x=NULL)

res$Case=sub('*:DD City','',as.character(res$WES))
res$Case=sub('*:DDU','',as.character(res$Case))
res$Case=substr(res$Case,nchar(res$Case),nchar(res$Case))

#ODDS between cases and controls (as random effects)
apply(res[res$Case==1,2:4],2,quantile,prob=0.5)
apply(res[res$Case==0,2:4],2,quantile,prob=0.5)

#ODDS between sites (as random effects)
apply(res[1:157,2:4],2,quantile,prob=0.5)#City
apply(res[158:290,2:4],2,quantile,prob=0.5)#DDU

#Finally odds city & cases VS city & controls (as random effects)
x=which(res$Case[1:157]==1)
apply(res[x,2:4],2,quantile,prob=0.5)
x=which(res$Case[1:157]==0)
apply(res[x,2:4],2,quantile,prob=0.5)

#Finally odds DDU & cases VS DDU & controls (as random effects)
x=which(res$Case[158:290]==1)
apply(res[c(158:290)[x],2:4],2,quantile,prob=0.5)
x=which(res$Case[158:290]==0)
apply(res[c(158:290)[x],2:4],2,quantile,prob=0.5)


