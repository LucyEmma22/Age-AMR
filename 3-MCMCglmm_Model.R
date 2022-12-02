setwd("~/OneDrive - University of Edinburgh/New mcmc model results")

library(MCMCglmm)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra) 
library(grid)  
library(HDInterval) 
library(modeest) 

##########################################################################################################
# MODELS 
##########################################################################################################

data<-read.csv("dataset_2.csv",header=TRUE)

# MODEL

nitt = 1100000
burnin = 100000
thin = 200

# FULL QUADRATIC MCMC MODEL
prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002),G2 = list(V = 1, nu = 0.002),G3 = list(V = diag(3), n = 3,alpha.mew=c(0,0,0),alpha.V=diag(3)*100),G4 = list(V = diag(3), n = 3,alpha.mew=c(0,0,0),alpha.V=diag(3)*100),G5 = list(V = diag(3), n = 3,alpha.mew=c(0,0,0),alpha.V=diag(3)*100)))
mcmc_model <- MCMCglmm(cbind(resistant, susceptible) ~ poly(age_scaled, 2, raw = TRUE),random = ~paper +paper_dataset +us(1 + poly(age_scaled, 2, raw = TRUE)):class +us(1 + poly(age_scaled, 2, raw = TRUE)):genus +us(1 + poly(age_scaled, 2, raw = TRUE)):class:genus, data = data, family = "multinomial2", prior = prior, nitt = nitt, burnin = burnin, thin = thin, pr = TRUE)

# FULL QUADRATIC GLMER MODEL
library(lme4)
model<-glmer(cbind(resistant, susceptible)~age_scaled+I(age_scaled^2)+(1|paper)+(1|paper_dataset)+(1+age_scaled+I(age_scaled^2)|genus)+(1+age_scaled+I(age_scaled^2)|class)+(1+age_scaled+I(age_scaled^2)|genus:class),data=data,family="binomial")

# FULL LINEAR MCMC MODEL
prior_linear <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002),G2 = list(V = 1, nu = 0.002),G3 = list(V = diag(2), n = 2,alpha.mew=c(0,0),alpha.V=diag(2)*100),G4 = list(V = diag(2), n = 2,alpha.mew=c(0,0),alpha.V=diag(2)*100),G5 = list(V = diag(2), n = 2,alpha.mew=c(0,0),alpha.V=diag(2)*100)))
mcmc_model_linear<- MCMCglmm(cbind(resistant, susceptible) ~ poly(age_scaled, 1,raw = TRUE),random = ~paper+paper_dataset+ us(1+poly(age_scaled, 1,raw = TRUE)):class+ us(1+poly(age_scaled, 1,raw = TRUE)):genus+ us(1+poly(age_scaled, 1,raw = TRUE)):class:genus, data=data,family="multinomial2",prior=prior_linear,nitt=nitt,burnin=burnin,thin=thin,pr=TRUE)

# QUADRATIC MCMC MODEL NO RANDOM EFFECTS
prior_no_random_effects <- list(R = list(V = 1, nu = 0.002))
mcmc_model_no_random_effects <- MCMCglmm(cbind(resistant, susceptible) ~ poly(age_scaled, 2, raw = TRUE), data = data,family = "multinomial2", prior = prior_no_random_effects, nitt = nitt, burnin = burnin, thin = thin, pr = TRUE)

save(list = c("mcmc_model", "mcmc_model_linear","mcmc_model_no_random_effects"),file="mcmc_models_new.Rdata")


# Gelman and Rubin

mcmc_model_2 <- MCMCglmm(cbind(resistant, susceptible) ~ poly(age_scaled, 2, raw = TRUE),random = ~paper +paper_dataset +us(1 + poly(age_scaled, 2, raw = TRUE)):class +us(1 + poly(age_scaled, 2, raw = TRUE)):genus +us(1 + poly(age_scaled, 2, raw = TRUE)):class:genus, data = data, start=list(QUASI=FALSE), family = "multinomial2", prior = prior, nitt = nitt, burnin = burnin, thin = thin, pr = TRUE)
mcmc_model_3 <- MCMCglmm(cbind(resistant, susceptible) ~ poly(age_scaled, 2, raw = TRUE),random = ~paper +paper_dataset +us(1 + poly(age_scaled, 2, raw = TRUE)):class +us(1 + poly(age_scaled, 2, raw = TRUE)):genus +us(1 + poly(age_scaled, 2, raw = TRUE)):class:genus, data = data, start=list(QUASI=FALSE), family = "multinomial2", prior = prior, nitt = nitt, burnin = burnin, thin = thin, pr = TRUE)

gelman.diag(mcmc.list(mcmc_model$Sol, mcmc_model_2$Sol, mcmc_model_3$Sol),multivariate=FALSE)
gelman.diag(mcmc.list(mcmc_model$Sol[,1, drop=TRUE], mcmc_model_2$Sol[,1, drop=TRUE], mcmc_model_3$Sol[,1, drop=TRUE]))
gelman.diag(mcmc.list(mcmc_model$Sol[,2, drop=TRUE], mcmc_model_2$Sol[,2, drop=TRUE], mcmc_model_3$Sol[,2, drop=TRUE]))
gelman.diag(mcmc.list(mcmc_model$Sol[,3, drop=TRUE], mcmc_model_2$Sol[,3, drop=TRUE], mcmc_model_3$Sol[,3, drop=TRUE]))

save(list = c("mcmc_model_2", "mcmc_model_3"),file="gelman_and_rubin_mcmc_models.Rdata")

##########################################################################################################
# POSTERIOR DISTRIBUTIONS
##########################################################################################################

# Calculating posterior distribution for each random effect by adding onto the fixed effect
all_results<-as.data.frame(mcmc_model$Sol) %>% select(matches("Intercept|age_scaled")) 

data$interaction<-paste0(data$class,":",data$genus)
druglist<-sort(unique(data$class))
buglist<-sort(unique(data$genus))
interactionlist<-sort(unique(paste0(data$class,":",data$genus)))

drugs<-c(paste0("Intercept_",druglist),paste0("Linear_",druglist),paste0("Quadratic_",druglist))
bugs<-c(paste0("Intercept_",buglist),paste0("Linear_",buglist),paste0("Quadratic_",buglist))
interaction<-c(paste0("Intercept_",interactionlist,"_Intercept"),paste0("Linear_",interactionlist,"_Linear"),paste0("Quadratic_",interactionlist,"_Quadratic"))

colnames(all_results)<-c("I_All","L_All","Q_All",drugs,bugs,interaction)

drugs2<-paste0(drugs,":")
bugs2<-c(paste0(":",buglist,"_Intercept"),paste0(":",buglist,"_Linear"),paste0(":",buglist,"_Quadratic"))

for (i in 1:length(drugs)){
  a<-select(all_results,drugs[i])
  all_results<-all_results %>% mutate_at(vars(contains(drugs2[i])), ~(.+a))
}

for (i in 1:length(bugs)){
  b<-select(all_results,bugs[i])
  all_results<-all_results %>% mutate_at(vars(contains(bugs2[i])), ~(.+b))
}

c<-select(all_results,I_All)
all_results<-all_results %>% mutate_at(vars(contains("Intercept")),( ~(.+c)))

d<-select(all_results,L_All)
all_results<-all_results %>% mutate_at(vars(contains("Linear")), ~(.+d))

e<-select(all_results,Q_All)
all_results<-all_results %>% mutate_at(vars(contains("Quadratic")), ~(.+e))

interaction2<-c(paste0("Intercept_",interactionlist),paste0("Linear_",interactionlist),paste0("Quadratic_",interactionlist))
colnames(all_results)<-c("Intercept_All","Linear_All","Quadratic_All",drugs,bugs,interaction2)

##########################################################################################################
# CALCULATE THE MEAN, CI AND VARIANCE FOR PARAMETERS FOR EVERY CLASS, GENUS AND COMBINATION
##########################################################################################################

average_parameters<-data.frame()
for (i in 1:ncol(all_results)){
  mean_and_CI<-data.frame(id=colnames(all_results)[i],
                          mean=mean(as.vector(unlist(all_results[,i]))),
                          mode=mfv(round(as.vector(unlist(all_results[,i])),2)),
                          lower=HPDinterval(mcmc(as.vector(unlist(all_results[,i]))))[1,1],
                          upper=HPDinterval(mcmc(as.vector(unlist(all_results[,i]))))[1,2],
                          variance=var(as.vector(unlist(all_results[,i]))))
  average_parameters<-rbind(average_parameters,mean_and_CI)
}

write.csv(average_parameters,file="average_parameters.csv",row.names=FALSE)

##########################################################################################################
# CALCULATE THE RESISTANCE FREQUENCY AT AGE 0-100 FOR EVERY EFFECT AT EVERY ITERATION (CURVES)
##########################################################################################################

# CURVES AND MINMAX+AGE: OVERALL, CLASS, GENUS AND INTERACTION

age<-seq(0,1,length=101)

all_drug_bug_interaction_list<-c("All",druglist,buglist,interactionlist)
all_average_curves<-data.frame()

for (j in 1:length(all_drug_bug_interaction_list)){
     
  # ITERATION CURVES: For each iteration, calculate resistance frequency for age 0-100 based on parameter estimates
  iteration_curves<-data.frame()
  
  for (i in 1:nrow(all_results)){
    c<-as.numeric(select(all_results,paste0("Intercept_",all_drug_bug_interaction_list[j]))[i,])
    m<-as.numeric(select(all_results,paste0("Linear_",all_drug_bug_interaction_list[j]))[i,])
    m2<-as.numeric(select(all_results,paste0("Quadratic_",all_drug_bug_interaction_list[j]))[i,])
    
    curves<- data.frame(age=age,resistance=exp(c+m*age+m2*I(age^2))/(1+exp(c+m*age+m2*I(age^2))))
    
    minmax<-data.frame(max=filter(curves,resistance==max(curves$resistance))$resistance,
      min=filter(curves,resistance==min(curves$resistance))$resistance,
      age_max=filter(curves,resistance==max(curves$resistance))$age,
      age_min=filter(curves,resistance==min(curves$resistance))$age,
      zero=filter(curves,age==0)$resistance,
      hundred=filter(curves,age==1)$resistance,
      max_difference=filter(curves,resistance==max(curves$resistance))$resistance-filter(curves,resistance==min(curves$resistance))$resistance,
      age_max_difference=filter(curves,resistance==max(curves$resistance))$age-filter(curves,resistance==min(curves$resistance))$age,
      zero_hundred_difference=filter(curves,age==1)$resistance-filter(curves,age==0)$resistance)
    
    iteration_curves<-rbind(iteration_curves,cbind(minmax,t(select(curves,resistance))))
  }
  rownames(iteration_curves)<-NULL
  
  # AVERAGE CURVES: For each age, calculate the mean, upper and lower CI and variance
  average_curves<-data.frame()
  for (i in 1:ncol(iteration_curves)){
    mean_and_CI<-data.frame(id=all_drug_bug_interaction_list[j],age=colnames(iteration_curves)[i],
                          mean=mean(iteration_curves[,i]),
                          mode=mfv(round(iteration_curves[,i],2)),
                          lower=HPDinterval(mcmc(iteration_curves[,i]))[1,1],
                          upper=HPDinterval(mcmc(iteration_curves[,i]))[1,2],
                          variance=var(iteration_curves[,i]))
    average_curves<-rbind(average_curves,mean_and_CI)
  }
  all_average_curves<-rbind(all_average_curves,average_curves)
  
}

write.csv(all_average_curves,file="all_average_curves_new.csv",row.names=FALSE)
