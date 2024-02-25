library(tidyr)
library(stringr)
library(dplyr)
library(gtools)

# # Read in CSV files and store as a single dataframe
# setwd("~/OneDrive - University of Edinburgh/Age_AMR/Age-AMR_Github/csv data files")
# file_list<-list.files()
# all_data<-data.frame()
# 
# for (i in 1:length(file_list)){
#   data<-data.frame(read.csv(file_list[i],header=T))[1:5] %>% drop_na(resistant)
#   colnames(data)<-c("min","max","age","susceptible","resistant")
#   data$paper_dataset<-gsub(".csv","",file_list[i])
#   all_data<-rbind(all_data,data)
# }
# all_data<-all_data %>% separate(paper_dataset,c("dataset","paper"),sep="_",remove=FALSE) %>% mutate(paper_dataset=gsub("_","",paper_dataset))
# setwd("~/OneDrive - University of Edinburgh/Age_AMR/Age-AMR_Github")
# write.csv(all_data,"all_data.csv",row.names=FALSE)

# "all_data.csv", "master_spreadsheet.csv", "Drug Classes.csv" and "genus.csv" are included as supporting files.
setwd("~/OneDrive - University of Edinburgh/Age_AMR/Age-AMR_Github")
all_data<-read.csv("all_data.csv")

# Calculate age where upper or lower limit was missing
all_data$age<-as.numeric(all_data$age)
data2<-na.omit(all_data)
minmodel<-lm(age~min,data=data2)
maxmodel<-lm(age~max,data=data2)

all_data$age<-ifelse(is.na(all_data$max)==TRUE,round(coef(minmodel)[1]+coef(minmodel)[2]*all_data$min,digits=1),all_data$age)
all_data$age<-ifelse(is.na(all_data$min)==TRUE,round(coef(maxmodel)[1]+coef(maxmodel)[2]*all_data$max,digits=1),all_data$age)
all_data$age<-ifelse(all_data$age<0,0,all_data$age)
all_data$age_scaled<-all_data$age/max(all_data$age)

# Taking out datasets with no resistance in any age class
no_resistance<-all_data %>% group_by(paper_dataset) %>% summarise(total_resistant = sum(resistant)) %>% filter(total_resistant==0) %>% select(paper_dataset)
all_data<-anti_join(all_data,no_resistance,by="paper_dataset")
all_data$paper<-as.numeric(all_data$paper)

# Define drug, bug, class and genus (ESBL and B-lactamase are tested using cephalosporins, drug for AmpC not defined) for each dataset
master_spreadsheet<-read.csv("master_spreadsheet.csv") %>% rename(paper=Paper) %>% separate_rows(Coder, sep="; ",convert = TRUE)  %>% separate(Coder,c("dataset","coder"),sep=" = ")  %>% separate(coder,c("drug","bug","coder"),sep="_")%>% 
  mutate(drug=ifelse(drug=="ESBL"|drug=="B-lactamase","cephalosporin",drug)) %>% filter(drug!="AmpC") %>%
  full_join(read.csv("Drug Classes.csv"),by="drug") %>% full_join(read.csv("genus.csv"),by="bug")

# Add drug/bug/class/genus, merge datasets which do not differ by drug or bug and remove datasets without a specific class and genus
all_data<- all_data %>% full_join(master_spreadsheet) %>% 
  filter(is.na(resistant)==FALSE) %>% 
  group_by(age,paper,drug,bug) %>% 
  mutate(susceptible = sum(susceptible), resistant = sum(resistant)) %>% 
  ungroup() %>%
  distinct(paper, drug, bug, age, .keep_all = TRUE)%>%
  filter(class!="Any" & genus!="Multiple" & class!="Multiple") 

# Additional filtering
all_data2<-filter(all_data,!class %in% c("Mupirocin", "Rifaximin")) %>% mutate(class=ifelse(class=="Fluoroquinolones","Quinolones",class)) # Remove classes which are not "J ANTIINFECTIVES FOR SYSTEMIC USE"
all_data2$paper<-as.factor(all_data2$paper)
all_data2<-all_data2[all_data2$resistant+all_data2$susceptible!=0,] # remove rows where there is no one in the age class
write.csv(all_data2,file="dataset_2.csv",row.names=FALSE)