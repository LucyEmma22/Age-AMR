setwd("~/OneDrive - University of Edinburgh/Age_AMR/Age-AMR_Github")
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(MCMCglmm) 
library(ggnewscale) 
library(patchwork) 
library(pals)

data<-read.csv("dataset_2.csv")

# Calculate number of papers, datasets adn samples for each class and genus
data$interaction<-paste0(data$class,":",data$genus)
data1<-data
# Number papers
number_of_papers_datasets_and_samples <- data1 %>% group_by(class,genus) %>% summarise(paper_count = n_distinct(paper)) %>% 
  rbind(data1 %>% group_by(class) %>% summarise(paper_count = n_distinct(paper)) %>% mutate(genus="All")) %>% 
  rbind(data1 %>% group_by(genus) %>% summarise(paper_count = n_distinct(paper)) %>% mutate(class="All")) %>% 
  rbind(data.frame(class="All",genus="All",paper_count=length(unique(data1$paper))))%>% 
  # Number datasets
  full_join(data1 %>% group_by(class,genus) %>% summarise(dataset_count = n_distinct(paper_dataset)) %>% 
              rbind(data1 %>% group_by(class) %>% summarise(dataset_count = n_distinct(paper_dataset)) %>% mutate(genus="All")) %>% 
              rbind(data1 %>% group_by(genus) %>% summarise(dataset_count = n_distinct(paper_dataset)) %>% mutate(class="All")) %>% 
              rbind(data.frame(class="All",genus="All",dataset_count=length(unique(data1$paper_dataset))))) %>% 
  # Number samples 
  full_join(data1 %>% group_by(class,genus) %>% summarise(sample_count = sum(susceptible)+sum(resistant)) %>% 
              rbind(data1 %>% group_by(class) %>% summarise(sample_count = sum(susceptible)+sum(resistant)) %>% mutate(genus="All")) %>% 
              rbind(data1 %>% group_by(genus) %>% summarise(sample_count = sum(susceptible)+sum(resistant)) %>% mutate(class="All")) %>% 
              rbind(data.frame(class="All",genus="All",sample_count=sum(data1$susceptible)+sum(data1$resistant))))

number_of_papers_datasets_and_samples$class <- factor(number_of_papers_datasets_and_samples$class, levels = rev(sort(as.vector(unique(number_of_papers_datasets_and_samples$class)))))
number_of_papers_datasets_and_samples$genus <- factor(number_of_papers_datasets_and_samples$genus, levels = sort(as.vector(unique(number_of_papers_datasets_and_samples$genus))))

# Plot number of papers, datasets and samples for each class and genus
sample_size_plot<-ggplot(filter(number_of_papers_datasets_and_samples,class!="All",genus!="All"), aes(x = genus, y = fct_relabel(class,str_wrap,width = 16)))+
  geom_tile(aes(fill=log(sample_count)))+
  geom_text(aes(label=paper_count),size=3,vjust = 0, nudge_y = 0.2)+
  geom_text(aes(label=dataset_count),size=3)+
  geom_text(aes(label=sample_count),size=2,vjust = 0, nudge_y = -0.35)+
  theme_classic()+
  scale_fill_gradient(high = "mediumpurple",low="white")+
  theme(legend.position = "none")+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x=element_text(angle=20,hjust=1))+
  theme(text=element_text(family="Helvetica",size=12))

# Load model results
load("~/Downloads/New mcmc model results/mcmc_models_new.Rdata")

# Plot Trace and density
plot(mcmc_model$VCV[,c(1,2,3,7,11,12,16,20,21,25,29,30)])
plot(mcmc_model$Sol[,1:3])

#plot.MCMCglmm(mcmc_model)
#plot(mcmc_model$VCV[,1:4])
# plot(mcmc_model$VCV[,1],main="Paper")
# plot(mcmc_model$VCV[,2],main="Dataset")
# plot(mcmc_model$VCV[,3],main="Class (Intercept)")
# plot(mcmc_model$VCV[,7],main="Class (Linear)")
# plot(mcmc_model$VCV[,11],main="Class (Quadratic)")
# plot(mcmc_model$VCV[,12],main="Genus (Intercept)")
# plot(mcmc_model$VCV[,16],main="Genus (Linear)")
# plot(mcmc_model$VCV[,20],main="Genus (Quadratic)")
# plot(mcmc_model$VCV[,21],main="Class-Genus (Intercept)")
# plot(mcmc_model$VCV[,25],main="Class-Genus (Linear)")
# plot(mcmc_model$VCV[,29],main="Class-Genus (Quadratic)")
# plot(mcmc_model$VCV[,30],main="Units")
# 
# random_effects<-data.frame(mcmc_model$VCV)[,c(1:3,7,11,12,16,20,21,25,29,30)] %>% mutate(iteration=1:5000)
# colnames(random_effects)<-c("Paper","Dataset","Class_Intercept","Class_Linear", "Class_Quadratic", "Genus_Intercept","Genus_Linear", "Genus_Quadratic", "Class:Genus_Intercept","Class:Genus_Linear", "Class:Genus_Quadratic","Units", "Iteration")
# random_effects<-random_effects %>% gather("Random_Effect","Value",-Iteration)
# trace<-ggplot(data=random_effects,aes(x=Iteration,y=Value))+facet_wrap(~Random_Effect,scales = "free")+geom_line()+theme_light()
# density<-ggplot(data=random_effects,aes(x=Value))+facet_wrap(~Random_Effect,scales = "free")+geom_density()+theme_light()
# trace / density

# Plot random effects histograms
par(mfrow = c(3,4))
hist(mcmc(mcmc_model$VCV)[,"paper"],main="Paper",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"paper_dataset"],main="Dataset",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"(Intercept):(Intercept).class"],main="Class (Intercept)",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"poly(age_scaled, 2, raw = TRUE)1:poly(age_scaled, 2, raw = TRUE)1.class"],main="Class (Linear)",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"poly(age_scaled, 2, raw = TRUE)2:poly(age_scaled, 2, raw = TRUE)2.class"],main="Class (Quadratic)",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"(Intercept):(Intercept).genus"],main="Genus (Intercept)",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"poly(age_scaled, 2, raw = TRUE)1:poly(age_scaled, 2, raw = TRUE)1.genus"],main="Genus (Linear)",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"poly(age_scaled, 2, raw = TRUE)2:poly(age_scaled, 2, raw = TRUE)2.genus"],main="Genus (Quadratic)",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"(Intercept):(Intercept).class:genus"],main="Interaction (Intercept)",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"poly(age_scaled, 2, raw = TRUE)1:poly(age_scaled, 2, raw = TRUE)1.class:genus"],main="Interaction (Linear)",xlab="Variance")
hist(mcmc(mcmc_model$VCV)[,"poly(age_scaled, 2, raw = TRUE)2:poly(age_scaled, 2, raw = TRUE)2.class:genus"],main="Interaction (Quadratic)",xlab="Variance")
par(mfrow=c(1,1)) 

# Plot Random Variance from mcmc model
fixed_posterior<-as.data.frame(summary(mcmc_model)$solutions) %>% select(1,2,3,5) %>% rename("Mode"=post.mean) %>% rename("Lower 95% CI"='l-95% CI') %>% rename("Upper 95% CI"='u-95% CI') %>% rename("P Value"='pMCMC') %>% mutate(Term=c(paste("\U03B2\U2080"),paste("\U03B2\U2081"),paste("\U03B2\U2081\u00B2"))) %>% remove_rownames()
random_variance<-as.data.frame(summary(mcmc_model)$Gcovariances[c(1,2,3,7,11,12,16,20,21,25,29),]) %>% rename("Variance"=post.mean) %>% mutate(Effect=c("Paper","Dataset","Antibiotic","Antibiotic","Antibiotic","Bacteria","Bacteria","Bacteria","Anitbiotic:Bacteria","Anitbiotic:Bacteria","Anitbiotic:Bacteria")) %>% mutate(Term=c(paste("\U03B2\U2080"),paste("\U03B2\U2080"),paste("\U03B2\U2080"),paste("\U03B2\U2081"),paste("\U03B2\U2081\u00B2"),paste("\U03B2\U2080"),paste("\U03B2\U2081"),paste("\U03B2\U2081\u00B2"),paste("\U03B2\U2080"),paste("\U03B2\U2081"),paste("\U03B2\U2081\u00B2"))) %>% remove_rownames()
random_variance$Effect <- factor(random_variance$Effect, levels = unique(random_variance$Effect))

random_variance_plot<-ggplot(random_variance, aes(y=Variance,x=Term))+
  geom_bar(position=position_dodge(), stat="identity",colour="white", aes(fill=Effect))+
  scale_fill_manual(name="Random Effect",values=c("grey40","grey70","mediumpurple","mediumseagreen","goldenrod"))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`,group=Effect), width=0.2,position=position_dodge(0.9))+
  labs(x = "Coefficient", y = "Variance")+
  theme_light()+
  theme(text=element_text(family="Helvetica",size=12))+
  theme(plot.title = element_text(hjust=-0.2))

# Plot overall curve from mcmc model
all_average_curves<-read.csv("all_average_curves_new.csv",header=T) %>% 
  mutate(age=(as.numeric(age)-1) * (max(data$age)/100)) %>% 
  na.omit()
p_overall<-ggplot(filter(all_average_curves,id=="All"),aes(age,mean))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper), outline.type = "both", fill = "grey", alpha = 0.2,linetype="dashed",colour="black")+
  labs(x = "Age", y = "Resistance Frequency")+
  theme_light()+
  theme(text=element_text(family="Helvetica",size=12))+
  theme(plot.title = element_text(hjust = -0.1))

# Plot all curves from mcmc model
all_average_curves<-all_average_curves%>% filter(str_detect(id, ":")) %>% separate(id,c("class","genus"),sep=":",remove=F) 
# Give each class a shape
shapes<-c(1,3,14,4,5,7,8,9,10,11,12, 13, 0, 2, 6, 15, 16, 17, 18)
names(shapes) <- sort(unique(all_average_curves$class))
# Give each genus a colour
colours<- c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", 
              "#1C8356", "#16FF32", "#F7E1A0", "#1CBE4F", 
              "#C4451C", "#DEA0FD", "#FE00FA", "#325A9B", "#FEAF16", 
              "#F8A19F", "#90AD1C", "#F6222E", "#1CFFCE", "#2ED9FF", 
              "#B10DA1", "#C075A6", "#FC1CBF", "#B00068", "#FBE426", 
              "#FA0087")
#colours<-c(polychrome()[1:25])
names(colours)<-sort(unique(all_average_curves$genus))

all_curve_plots<-ggplot(all_average_curves,aes(age,mean))+
  geom_point(data=data,aes(x=age,y=resistant/(resistant+susceptible),colour=genus,shape=class,fill=genus),size=0.8)+
  geom_line(aes(colour=genus))+
  geom_ribbon(aes(ymin =lower, ymax = upper,fill=genus,colour=genus), outline.type = "both", alpha = 0.2,linetype="dashed",show.legend=F)+
  ylim(0,1)+
  xlim(0,100)+
  labs(x="Age", y="Resistance Frequency")+
  theme_minimal()+
  scale_colour_manual(values=colours)+
  scale_fill_manual(values=colours)+
  scale_shape_manual(values=shapes)+
  facet_grid(rows=vars(str_wrap(class,15)),cols=vars(genus),switch = "both",drop=T)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position="none", strip.text.y.left = element_text(angle = 0,hjust=1),strip.text.x.bottom = element_text(angle = 20))

# Plot 4 example curves
examples<-c("Penicillins:Proteus","Polymyxins:Pseudomonas","Quinolones:Acinetobacter","Cephalosporins:Streptococcus")
example_plot_data<- data %>% rename(id=interaction) %>% filter(id %in% examples) 
example_curves<-filter(all_average_curves,id %in% examples) 
example_plots<-ggplot(example_curves,aes(age,mean))+
  geom_point(data=example_plot_data,aes(x=age,y=resistant/(resistant+susceptible),colour=genus,shape=class,fill=genus),size=0.8)+
  geom_line(aes(colour=genus))+
  geom_ribbon(aes(ymin =lower, ymax = upper,fill=genus,colour=genus), outline.type = "both", alpha = 0.2,linetype="dashed",show.legend=F)+
  xlim(0,100)+
  labs(x="Age", y=str_wrap("Resistance Frequency",20))+
  theme_minimal()+
  scale_colour_manual(values=colours)+
  scale_fill_manual(values=colours)+
  scale_shape_manual(values=shapes)+
  facet_grid(rows=vars(factor(id, levels = examples)),scales="free_y")+
  theme(legend.position="none", strip.text.y=element_blank(),text=element_text(family="Helvetica",size=10))

# Calculate maximum fold change in resistance frequency for all curves
minmax<-all_average_curves %>% group_by(id) %>% 
  summarise(min=min(mean),max=max(mean),var=sum(variance),mean=mean,age=age) %>% distinct() %>% 
  separate(id,c("class","genus"),sep=":",remove=F) %>% filter(min==mean|max==mean) %>% mutate(x=ifelse(mean==min,"age_min","age_max")) %>% 
  select(!mean) %>% spread(x,age) %>% 
  mutate(fold_change=max/min) %>% mutate(abs_change=(max - min)) %>%
  mutate(directional_change=ifelse(age_min>age_max,min-max,max-min)) %>% 
  mutate(directional_fold_change=ifelse(age_min>age_max,min/max,max/min)) %>% ungroup() 

# Plot maximum fold change in resistance frequency for all curves
minmax$class <- factor(minmax$class, levels = rev(sort(as.vector(unique(minmax$class)))))
curve_bubbleplot<-ggplot(minmax, aes(x = genus, y = fct_relabel(class,str_wrap,width = 16)))+
  geom_point(aes(size=log(1/var),fill=log(directional_fold_change)),shape=21,stroke=0.5,colour="grey")+
  scale_size(name= "Log (1/Variance)",range = c(1,10))+
  theme_light()+
  scale_fill_gradient2(name="Log (Fold Change in Resistance Frequency)",low="mediumseagreen",mid="white",high = "mediumpurple",midpoint=0)+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x=element_text(angle=20,hjust=1),text=element_text(family="Helvetica",size=10),
        legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="grey"),legend.title.align=0.5,legend.position="bottom",legend.box="vertical")+
  guides(size = guide_legend(override.aes = list(colour = "grey",shape=21,stroke=1)))+
  new_scale("colour") +
  geom_point(data=filter(minmax,id %in% examples),aes(x = genus, y = fct_relabel(class,str_wrap,width = 16),colour=genus,size=log(1/var)+1,shape=class),stroke=1,show.legend=FALSE) +
  scale_colour_manual(values=colours)+
  scale_shape_manual(values=shapes)

# Plot linear and quadratic term for each class-genus combination 
average_parameters<-read.csv("average_parameters.csv",header=TRUE) %>% 
  separate(id,c("term","id"),sep="_") %>% 
  select(term,id,mean) %>% 
  distinct() %>% filter(str_detect(id, ":")) %>%
  spread(term,mean) %>% separate(id,c("class","genus"),sep=":",remove=F) 
point_graph<-ggplot(data = filter(average_parameters,!id %in% examples),aes(Linear,Quadratic))+
  geom_point(aes(col=genus,shape=class),size=1)+
  geom_point(data=filter(average_parameters,id %in% examples),aes(col=genus,shape=class),size=3,stroke=1,show.legend=F)+
  labs(x = paste("\U03B2\U1D62"), y = paste("\U03B2\U1D62\U00B2"))+
  theme_light()+
  geom_hline(yintercept=0, linetype="solid",colour="grey")+
  geom_vline(xintercept=0, linetype="solid",colour="grey")+
  geom_function(fun = function(x) (-x/2),colour="grey",linetype="solid")+
  scale_colour_manual(name="Bacteria",values=colours)+
  scale_shape_manual(name="Antibiotic",values=shapes)+
  guides(colour=guide_legend(nrow=4,byrow=TRUE,keyheight = 0.1,keywidth=0.1))+
  guides(shape=guide_legend(nrow=4,byrow=TRUE,keyheight = 0.1,keywidth=0.1))+
  theme(text=element_text(family="Helvetica",size=10),
        legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="grey"),legend.title.align=0.5,
        legend.position="bottom",legend.box="vertical")

# Arranging Plots
sample_size_plot + plot_annotation(title="B")
p_overall + random_variance_plot + plot_annotation(tag_levels = 'A') 
point_graph + example_plots + curve_bubbleplot + plot_layout(widths = c(5,1,4)) + plot_annotation(tag_levels = 'A')
all_curve_plots
