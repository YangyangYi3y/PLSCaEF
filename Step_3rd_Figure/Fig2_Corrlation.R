rm(list=ls())
library(R.matlab)
library(ggplot2)
library(lattice)
library(ggpubr)
library(ggpmisc)
library(patchwork)
library(grid)
library(tidyverse)
library(plyr)
library(dplyr)
library(PupillometryR)
library(hrbrthemes)

windowsFonts(arial = windowsFont(family = "Arial"))        

comp1_corlor = '#4472c4'
comp1_shade = '#b4c7e7'
comp2_corlor = '#ed7d31'
comp2_shade = '#f8cbad'
# draw for each component
component = 1
fold_num = 1
Folder = paste0('D:/projects/ABCD/EF/PLS/results/ef_task23/combat/2fold/comp',as.character(component));



df = as.data.frame(matrix(nrow=0,ncol=3))
names(df)=c('FC','EF','Fold')
for (i in 0:fold_num){
    corrmatt = readMat(paste0(Folder, '/Fold_',as.character(i),'_Score.mat'));
    BrainCat_temp = corrmatt$Brain.test.ca[,component];
    BehaviorCat_temp = corrmatt$Behavior.test.ca[,component];
    data_temp<-data.frame(BrainCat_temp,BehaviorCat_temp)
    num = length(BrainCat_temp)
    label <- as.character(i+1)
    data_temp$Fold <- matrix (label,num)
    names(data_temp)=c('FC','EF','Fold')
    df = rbind(df,data_temp)
}

  
range(df$FC)
range(df$EF)
dev.off()
dev.new()

#combine two folder
ggscatter(df, x = 'FC', y = 'EF',cex.axis=1.5,cex.lab=3,
          color = comp1_corlor, shape = 20,size = 1.5, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = comp1_corlor,fill=comp1_shade), # Customize reg. line
          conf.int = TRUE )+ # Add confidence interval
  labs(x = "Functional Connectivity", y = "Executive Function") +
  theme(axis.title.x = element_text(margin=margin(t=10),face="bold", color='black', size=20,family = "arial"),
        axis.text.x = element_text(margin=margin(t=0), size=12, angle=0))+ 
  theme(axis.title.y = element_text(margin=margin(r=20),angle=90, face="bold", size=20,color='black',family = "arial"),
        axis.text.y = element_text(margin=margin(r =0),size=12)) +
  scale_y_continuous(limits = c(-8,6), breaks = seq(-8,6,2))+
  scale_x_continuous(limits = c(-100, 70), breaks = seq(-100,70,20)) 
ggsave(paste0(Folder,'/EF_FC_corr.tiff'), width = 20, height = 15, dpi = 300, units = "cm");

#draw for each component
dev.off()
dev.new()
p <- ggplot(data=df,cex.axis=1.5,
            mapping = aes(
            x = FC,
            y = EF,
            color = Fold,
            fill = Fold,se = FALSE))+
            theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
            theme(panel.border = element_blank())+
            theme(axis.line = element_line(colour = "black")) 
p
p + geom_point() +
  geom_smooth(method="lm",se = FALSE)+
  scale_color_manual(values = c(comp2_corlor, comp2_shade))+
  labs(x = "Functional Connectivity", y = "Executive Function") +
  theme(axis.title.x = element_text(margin=margin(t=10),face="bold", color='black', size=20,family = "arial"),
        axis.text.x = element_text(margin=margin(t=0), size=12))+ 
  theme(axis.title.y = element_text(margin=margin(r=20),angle=90, face="bold", size=20,color='black',family = "arial"),
        axis.text.y = element_text(margin=margin(r =0),size=12)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,2))+
  scale_x_continuous(limits = c(-70, 60), breaks = seq(-70,60,20))

ggsave(paste0(Folder,'/EF_FC_corr_2.tiff'), width = 20, height = 15, dpi = 300, units = "cm")


#draw dense map
corrper <- readMat('D:/projects/ABCD/EF/PLS/results/ef_task23/combat/2fold/Comp_corr.mat')
comp1_act <-as.data.frame(corrper$Corr.Actual[,1])
names(comp1_act)<-'Value1_a'
comp2_act <- as.data.frame(corrper$Corr.Actual[,2])
names(comp2_act)<-'Value2_a'
comp1_per <-as.data.frame(corrper$Corr.Permutation[,1])
names(comp1_per)<-'Value1_p'
comp2_per = as.data.frame(corrper$Corr.Permutation[,2])
names(comp2_per)<-'Value2_p'

a_list = list(comp1_act = corrper$Corr.Actual[,1],
              comp2_act = corrper$Corr.Actual[,2],
              comp1_per = corrper$Corr.Permutation[,1] ,
              comp2_per = corrper$Corr.Permutation[,2])
a_df <- do.call(cbind, lapply(lapply(a_list, unlist), `length<-`, max(lengths(a_list))))
a_df <- data.frame(a_df)

comp1 = 0.209
comp2 = 0.178
dev.off()
dev.new()
p <- ggplot(a_df,aes(x=x)) +
  geom_density(aes(x = comp1_act, y = ..density..), fill=comp1_corlor) +
  geom_density(aes(x = comp1_per, y = ..density..), fill= comp1_shade)+
  labs(x = "R", y = "Density") +
  theme(axis.title.x = element_text(margin=margin(t=10),face="bold", color='black', 
                                    size=50,family = "arial"),
        axis.text.x = element_text(face="bold", color='black',family = "arial",size=2))+ 
  theme(axis.title.y = element_text(margin=margin(r=20),angle=90, face="bold", size=5,
                                    color='black',family = "arial"),
        axis.text.y = element_text(face="bold", color='black',family = "arial",size=5)) + 
  geom_text(aes(x=comp1, y=55, label="Actual correlation"), color=comp1_corlor, size=5) +
  geom_vline(aes(xintercept=comp1),
                color='#240101', linetype="dashed", size=0.5)+
  
  geom_text(aes(x=0, y=20, label="Permulation correlation"), color=comp1_shade, size=5)+ 
  geom_vline(aes(xintercept=median(comp1_per)),
            color='#240101', linetype="dashed", size=0.5)+
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,5),expand=c(0,0))+
  scale_x_continuous(limits = c(-0.15, 0.3), breaks = seq(-0.15,0.25,0.1))+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black")) 
p
ggsave(paste0(Folder,'/EF_FC_corr.tiff',as.character(component)), width = 20, height = 15, dpi = 300, units = "cm")

