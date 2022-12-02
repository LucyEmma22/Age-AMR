setwd("~/OneDrive - University of Edinburgh/New mcmc model results")
library(lme4)
library(tidyr)
library(dplyr)
library(ggplot2)

# Import Processed Data
data<-read.csv("dataset_2.csv",header=T)

# Run regression for each paper to estimate of effect size and standard error
paper_list<-distinct(select(data,paper))
all_results<-data.frame()
for(i in 1:nrow(paper_list)){
  specific_data<-filter(data,paper==paper_list[i,1])
  model<-glm(cbind(resistant, susceptible)~age_scaled,data=specific_data,family="binomial")
  result<-data.frame(slope=summary(model)$coeff[2,1],error=summary(model)$coeff[2,2],paper=paper_list[i,1])
  all_results<-rbind(all_results,result)
}

# Funnel Plot
ggplot(all_results,aes(x=slope,y=log(1/error)))+
  geom_point(size=0.8)+
  geom_vline(xintercept=0,colour="black",linetype="dashed")+
  theme_light()+
  labs(title="Funnel Plot to Assess Publication Bias",x="Effect Size",y="1/Standard Error")+
  theme(plot.title = element_text(hjust = 0.5))
