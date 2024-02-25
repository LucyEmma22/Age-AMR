library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(deSolve)
library(grid)
library(gridtext)
library(patchwork)

population_y<-population_m<-population_o<-1
ct_y<-ct_m<-ct_o<-0.05
a_y<-a_m<-a_o<-1
cg_y<-cg_m<-cg_o<-0
g<-1/30
variable_parameters<-data.frame(
  
t_y=c(0.001, rep(c(0,     0.002, 0.0005, 0.00125, rep(0.001,8)),2)),
t_m=c(0.001, rep(c(0.001, 0.001, 0.002,  0.0005,  rep(0.001,8)),2)),
t_o=c(0.001, rep(c(0.002, 0,     0.0005, 0.00125, rep(0.001,8)),2)),

b_y=c(1/24, rep(c(rep(1/24,4),(1/(1-0.1))*(1/30), (1/(1-0.3))*(1/30), (1/(1-0.15))*(1/30), (1/(1-0.25))*(1/30),rep(1/24,4)),2)),
b_m=c(1/24,rep(c(rep(1/24,4),(1/(1-0.2))*(1/30), (1/(1-0.2))*(1/30), (1/(1-0.3))*(1/30), (1/(1-0.1))*(1/30),rep(1/24,4)),2)),
b_o=c(1/24,rep(c(rep(1/24,4),(1/(1-0.3))*(1/30), (1/(1-0.1))*(1/30), (1/(1-0.15))*(1/30), (1/(1-0.25))*(1/30),rep(1/24,4)),2)),

g_y=c(1/30,rep(c(rep(1/30,8),(1/24)/(1/(1-0.1)), (1/24)/(1/(1-0.3)), (1/24)/(1/(1-0.15)), (1/24)/(1/(1-0.25))),2)),
g_m=c(1/30,rep(c(rep(1/30,8),(1/24)/(1/(1-0.2)), (1/24)/(1/(1-0.2)), (1/24)/(1/(1-0.3)), (1/24)/(1/(1-0.1))),2)),
g_o=c(1/30,rep(c(rep(1/30,8),(1/24)/(1/(1-0.3)), (1/24)/(1/(1-0.1)), (1/24)/(1/(1-0.15)), (1/24)/(1/(1-0.25))),2)),

m=c(0, rep(0.01,12),rep(0.05,12)),
run=c(0, rep(1:4,6)),
varying=c("None", rep(c(rep("Treatment",4),rep("Transmission",4),rep("Clearance",4)),2))) %>% 
  mutate(prev_y=1-(g_y/b_y),prev_m=1-(g_m/b_m),prev_o=1-(g_o/b_o))


model<-function(times,state,pars){
  with(as.list(c(state,pars)),{
    
    b1_y<-b_y*((1-assortivity*2)*(s_y+0.5*sr_y)+assortivity*((s_m+0.5*sr_m)+(s_o+0.5*sr_o)))
    b2_y<-(1-ct_y)*b_y*((1-assortivity*2)*(r_y+0.5*sr_y)+assortivity*((r_m+0.5*sr_m)+(r_o+0.5*sr_o)))
    u_y<-1-(s_y+r_y+sr_y)
    ds_y<-b1_y*u_y+cg_y*sr_y-(b2_y+g_y+a_y*t_y)*s_y
    dr_y<-b2_y*u_y+(a_y*t_y)*sr_y-(b1_y+g_y)*r_y
    dsr_y<-b2_y*s_y+b1_y*r_y-(g_y+a_y*t_y+cg_y)*sr_y
    
    b1_m<-b_m*((1-assortivity*2)*(s_m+0.5*sr_m)+assortivity*((s_y+0.5*sr_y)+(s_o+0.5*sr_o)))
    b2_m<-(1-ct_m)*b_m*((1-assortivity*2)*(r_m+0.5*sr_m)+assortivity*((r_y+0.5*sr_y)+(r_o+0.5*sr_o)))
    u_m<-1-(s_m+r_m+sr_m)
    ds_m<-b1_m*u_m+(cg_m)*sr_m-(b2_m+g_m+a_m*t_m)*s_m
    dr_m<-b2_m*u_m+(a_m*t_m)*sr_m-(b1_m+g_m)*r_m
    dsr_m<-b2_m*s_m+b1_m*r_m-(g_m+a_m*t_m+cg_m)*sr_m
    
    b1_o<-b_o*((1-assortivity*2)*(s_o+0.5*sr_o)+assortivity*((s_m+0.5*sr_m)+(s_y+0.5*sr_y)))
    b2_o<-(1-ct_o)*b_o*((1-assortivity*2)*(r_o+0.5*sr_o)+assortivity*((r_m+0.5*sr_m)+(r_y+0.5*sr_y)))
    u_o<-1-(s_o+r_o+sr_o)
    ds_o<-b1_o*u_o+(cg_o)*sr_o-(b2_o+g_o+a_o*t_o)*s_o
    dr_o<-b2_o*u_o+(a_o*t_o)*sr_o-(b1_o+g_o)*r_o
    dsr_o<-b2_o*s_o+b1_o*r_o-(g_o+a_o*t_o+cg_o)*sr_o
    
    list(c(ds_y,dr_y,dsr_y,ds_m,dr_m,dsr_m,ds_o,dr_o,dsr_o))
  })
}

runtime<-3650
times<-seq(0,runtime,length.out=runtime+1) 

state<-c(s_y=0.2*0.9,
         r_y=0.2*0.1,
         sr_y=0,
         s_m=0.2*0.9,
         r_m=0.2*0.1,
         sr_m=0,
         s_o=0.2*0.9,
         r_o=0.2*0.1,
         sr_o=0)

#########################################################################################################################

results<-data.frame()
for (i in 1:nrow(variable_parameters)){
  t_y<-variable_parameters$t_y[i]
  t_m<-variable_parameters$t_m[i]
  t_o<-variable_parameters$t_o[i]
  b_y<-variable_parameters$b_y[i]
  b_m<-variable_parameters$b_m[i]
  b_o<-variable_parameters$b_o[i]
  g_y<-variable_parameters$g_y[i]
  g_m<-variable_parameters$g_m[i]
  g_o<-variable_parameters$g_o[i]
  assortivity<-variable_parameters$m[i]
  
  pars<-c(b_y,b_m,b_o,g_y,g_m,g_o,cg_y,cg_m,cg_o,a_y,a_m,a_o,t_y,t_m,t_o,ct_y,ct_m,ct_o)
  out<-as.data.frame(ode(y=state,times=times,func=model,parms=pars))
  
  finals<-mutate(out[nrow(out),], 
                 prev_y=s_y+r_y+sr_y, 
                 prev_m=s_m+r_m+sr_m, 
                 prev_o=s_o+r_o+sr_o, 
                 res_y=(r_y+0.5*sr_y)/(r_y+s_y+sr_y), 
                 res_m=(r_m+0.5*sr_m)/(r_m+s_m+sr_m), 
                 res_o=(r_o+0.5*sr_o)/(r_o+s_o+sr_o))
  
  results<-rbind(results,data.frame(run=variable_parameters$run[i],
                                    varying=variable_parameters$varying[i],
                                    age=c("Child","Adult","Elderly"),
                                    mixing=assortivity,
                                    Prevalence=c(finals$prev_y,finals$prev_m,finals$prev_o), 
                                    Resistance=c(finals$res_y, finals$res_m, finals$res_o), 
                                    Treatment=c(t_y,t_m,t_o),
                                    Transmission=c(b_y,b_m,b_o),
                                    Clearance=c(g_y,g_m,g_o)))
          }

results2<-results %>% mutate(mixing=paste0("Mixing = ",mixing))
results<-results %>% gather("parameter","value",5:9) %>% filter((varying==parameter&mixing==0.01)|parameter=="Resistance"|parameter=="Prevalence") %>% mutate(mixing=ifelse(!(parameter %in% c("Resistance", "Prevalence")),"",paste0("Mixing = ",mixing)))
results$age <- factor(c('Child', 'Adult', 'Elderly'), levels = c('Child', 'Adult', 'Elderly'))

# Treatment
treatment1<-ggplot(results2 %>% filter(varying=='Treatment', mixing=='Mixing = 0.01'), aes(x=age, y=Treatment, group=run))+
  geom_point(size=0.5, col='mediumpurple')+
  geom_line(col='mediumpurple')+
  theme_light()+
  labs(x=NULL, y='Treatment Rate')+
  facet_grid(rows=vars(run),cols=vars(mixing))+
  labs(tag = "A")+
  theme(strip.placement = "outside",
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank())

treatment2<-ggplot(results2 %>% filter(varying=='Treatment'), aes(x=age, y=Resistance, group=run))+
  geom_point(size=0.5)+
  geom_line()+
  theme_light()+
  labs(x=NULL, y='Resistance Probability')+
  facet_grid(rows=vars(run),cols=vars(mixing))+
  theme(strip.placement = "outside",
        strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black'),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank())

# Transmission
transmission1<-ggplot(results2 %>% filter(varying=='Transmission', mixing=='Mixing = 0.01'), aes(x=age, y=Transmission, group=run))+
  geom_point(size=0.5, col='mediumseagreen')+
  geom_line(col='mediumseagreen')+
  theme_light()+
  labs(x=NULL, y='Transmission Rate')+
  facet_grid(rows=vars(run),cols=vars(mixing))+
  labs(tag = "B")+
  theme(strip.placement = "outside",
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank())

transmission2<-ggplot(results2 %>% filter(varying=='Transmission'), aes(x=age, y=Resistance, group=run))+
  geom_point(size=0.5)+
  geom_line()+
  theme_light()+
  labs(x=NULL, y='Resistance Probability')+
  facet_grid(rows=vars(run),cols=vars(mixing))+
  theme(strip.placement = "outside",
        strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black'),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank())

# Clearance
clearance1<-ggplot(results2 %>% filter(varying=='Clearance', mixing=='Mixing = 0.01'), aes(x=age, y=Clearance, group=run))+
  geom_point(size=0.5, col='goldenrod')+
  geom_line(col='goldenrod')+
  theme_light()+
  labs(x=NULL, y='Clearance Rate')+
  facet_grid(rows=vars(run),cols=vars(mixing))+
  labs(tag = "C")+
  theme(strip.placement = "outside",
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank())

clearance2<-ggplot(results2 %>% filter(varying=='Clearance'), aes(x=age, y=Resistance, group=run))+
  geom_point(size=0.5)+
  geom_line()+
  theme_light()+
  labs(x=NULL, y='Resistance Probability')+
  facet_grid(rows=vars(run),cols=vars(mixing))+
  theme(strip.placement = "outside",
        strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black'),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank())

(treatment1 + treatment2 + plot_layout(widths = c(1,2))) / (transmission1 + transmission2 + plot_layout(widths = c(1,2))) / (clearance1 + clearance2 + plot_layout(widths = c(1,2)))
