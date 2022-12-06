setwd("~/OneDrive - University of Edinburgh/Age_AMR/Age_AMR_Github")
library(lme4)
library(ROCR)
library(cvAUC)
library(dplyr)

# Function to carry out cvAUC
iid_example <- function(data, V = 10){
  
  .cvFolds <- function(Y, V){  #Create CV folds (stratify by outcome)
    Y0 <- split(sample(which(Y==0)), rep(1:V, length = length(which(Y==0))))
    Y1 <- split(sample(which(Y==1)), rep(1:V, length = length(which(Y==1))))
    folds <- vector("list", length=V)
    for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}		
    return(folds)
  }
  .doFit <- function(v, folds, data){  #Train/test glm for each fold
    fit <- glm(Y~age_scaled+I(age_scaled^2), data = data[-folds[[v]],], family = binomial)
    pred <- predict(fit, newdata = data[folds[[v]],], type = "response")
    return(pred)
  }
  
  folds <- .cvFolds(Y = data$Y, V = V)  #Create folds
  predictions <- as.numeric(unlist(sapply(seq(V), .doFit, folds = folds, data = data)))  #CV train/predict
  predictions[unlist(folds)] <- predictions  #Re-order pred values
  # Get CV AUC and confidence interval
  out <- ci.cvAUC(predictions = predictions, labels = data$Y, 
                  folds = folds, confidence = 0.95)
  return(out)
}

# Group data by class and genus and filter out combinations with low sample size
data<-read.csv("dataset_2.csv",header=TRUE)

shapes<-c(1,3,14,4,5,7,8,9,10,11,12, 13, 0, 2, 6, 15, 16, 17, 18)
names(shapes) <- sort(unique(data$class))
# Give each genus a colour
colours<- c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", 
            "#1C8356", "#16FF32", "#F7E1A0", "#1CBE4F", 
            "#C4451C", "#DEA0FD", "#FE00FA", "#325A9B", "#FEAF16", 
            "#F8A19F", "#90AD1C", "#F6222E", "#1CFFCE", "#2ED9FF", 
            "#B10DA1", "#C075A6", "#FC1CBF", "#B00068", "#FBE426", 
            "#FA0087")
#colours<-c(polychrome()[1:25])
names(colours)<-sort(unique(data$genus))

data<-data.frame(data %>% dplyr::select(age_scaled,susceptible,resistant,class,genus) %>% group_by(age_scaled,class,genus) %>% summarise(resistant=sum(resistant),susceptible=sum(susceptible),age_scaled=age_scaled,class=class,genus=genus) %>% ungroup() %>% distinct())
drugbug<-data.frame(table(data$class,data$genus))%>% filter(Freq>=10) %>% rename("class"=Var1,"genus"=Var2) %>% dplyr::select(-Freq)

#View(data %>% mutate(count=1) %>% group_by(class,genus) %>% summarise(rfreq=sum(resistant)/(sum(susceptible)+sum(resistant)),total_s=sum(susceptible),total_r=sum(resistant),total=sum(resistant)+sum(susceptible),points=sum(count)))

# Run cvAUC function for each class-genus combination
cvAUC_results<-data.frame()
for (i in 1:nrow(drugbug)){
specific_data<-filter(data,class==drugbug$class[i], genus==drugbug$genus[i])

data_long<-data.frame()
for (j in 1:nrow(specific_data)){
  data_long<-rbind(data_long,data.frame(
    age_scaled=specific_data$age_scaled[j],
    Y=c(rep(1,specific_data$resistant[j]),rep(0,specific_data$susceptible[j]))))}

# Get performance
set.seed(1)
out <- iid_example(data = data_long, V = 10) 
cvAUC_results<-rbind(cvAUC_results,data.frame(class=unique(specific_data$class),genus=unique(specific_data$genus),cvAUC=out$cvAUC,se=out$se,ci_l=out$ci[1],ci_u=out$ci[2]))
}

library(ggplot2)
ggplot(cvAUC_results,aes(cvAUC,paste(genus,class)))+
  geom_point(aes(colour=genus,shape=class))+
  geom_segment(aes(x = ci_l, y = paste(genus,class), xend = ci_u, yend = paste(genus,class),colour=genus), data = cvAUC_results)+
  geom_vline(xintercept=0.5,linetype="dashed")+
  geom_vline(xintercept=0.7,linetype="dashed")+
  scale_colour_manual(values=colours)+
  scale_fill_manual(values=colours)+
  scale_shape_manual(values=shapes)+
  labs(y=NULL)+
  theme_light()+
  guides(color = "none", shape = "none")
