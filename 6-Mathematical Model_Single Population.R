library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(deSolve)
library(grid)
library(gridtext)
library(patchwork)

population<-1
ct<-0.05
a<-1
cg<-0
t<-0.001
b<-1/24
g<-1/30

model<-function(times,state,pars){
  with(as.list(c(state,pars)),{
    
    b1<-b*(s+0.5*sr)
    b2<-(1-ct)*b*(r+0.5*sr)
    u<-1-(s+r+sr)
    ds<-b1*u+cg*sr-(b2+g+a*t)*s
    dr<-b2*u+(a*t)*sr-(b1+g)*r
    dsr<-b2*s+b1*r-(g+a*t+cg)*sr

    list(c(ds,dr,dsr))
  })
}

runtime<-3650
times<-seq(0,runtime,length.out=runtime+1) 

state<-c(s=0.2*0.9,
         r=0.2*0.1,
         sr=0)

pars<-c(b,g,cg,a,t,ct)
out<-as.data.frame(ode(y=state,times=times,func=model,parms=pars))
finals<-mutate(out[nrow(out),],prev=s+r+sr,res=(r+0.5*sr)/(r+s+sr))
print(paste0('Prevalence = ', finals$prev))
print(paste0('Resistance = ', finals$res))
