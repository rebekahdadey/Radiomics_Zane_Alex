#Radiomics analysis for Zane and Alex
#Re-do analysis with melanoma patients updated
#9-21-23

#average First order Lung minimum, first order liver minimum 

library(ggplot2)
library(dplyr)
library(scales)
library(ctc)
library(broom)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(patchwork)
library(readxl)
library("openxlsx")
library(DescTools)
library(tidyr)
library("xlsx")
library(glmnet) #lasso software
library(tidyverse)
library("caret")
library(xgboost)
library(devtools)
#library(SHAPforxgboost)
library(e1071)
library(glmnet)
library(stringr)
library("rstatix")
#devtools::install_github("liuyanguu/SHAPforxgboost")

#ok note, the serafettin has more data in it.
setwd("/ix/jluke/Projects/becky/Radiomics_Zane_Alex/")

#moved all output to 9-26-22 and 2-10-23 data -----------------------------------------
#this is new analysis not in these folders
##-----------------------------------------------------------------------------------------



##FIX clinical stuff
#9-12-12
Data<-readxl::read_excel("Total_ICI_data_corrected_09162023.xlsx", sheet="Master")
str(Data)
############################################################################
######################################################################################
##re-run Figure for Zane with monootherapy vs all together 1-25-23, was put in 2-10-23 file folder
#now i am doing with the updated cohort for alex 8-11-23


#Filter samples to remove NE and N/A RECIST--------------------------------------
table(Data$`Best Response by RECIST`)
#CR  PD  PR  SD 
#20 115  60  47 

sum(is.na(Data$`Best Response by RECIST`))
#none that are NAs


#remove samples with NE
#rm(x)
#x = which(!is.na(Data$`Best Response by RECIST`) &
#            Data$`Best Response by RECIST`=='NE')
#print(paste0('index of NE row is:', x))

#Data$`Best Response by RECIST`[x] = NA

table(Data$`Best Response by RECIST`)
# CR  PD  PR  SD 
#20 115  60  47 

#no samples that have N/A as their RECIST score, but if did, use code below

#Data<-Data %>% filter(Data$`Best Response by RECIST`!="N/A")

##Collapse Responses into NR or R----------------------------------------
#according to Zane, keep SD as responder
Data$BRR = NA
Data$BRR[which(Data$`Best Response by RECIST` %in% c('PD'))] = 'NR'
Data$BRR[which(Data$`Best Response by RECIST` %in% c('PR','CR','SD'))] = 'R'

table(Data$BRR)
#NR   R 
#115 127 

table(Data$BRR)[2] / nrow(Data[!is.na(Data$BRR),])
#0.5247934 

sum(is.na(Data$BRR))
## 0 (removed)

saveRDS(Data, "Data.rds")

#add a new column with P001 as 1st entry, etc.
Data$PID = paste0('P',stringr::str_pad(1:nrow(Data), 3, pad = "0"))
tail(Data)
Data<-as.data.frame(Data)
my.levels<-c("PD","SD","PR","CR")
Data$`Best Response by RECIST`=factor(Data$`Best Response by RECIST`, levels=my.levels)
table(Data$`Target lesion % change from baseline`, Data$BRR)
Data$`Target lesion % change from baseline`<-as.numeric(Data$`Target lesion % change from baseline`)
str(Data)
Data<-Data %>% filter(!is.na(`Target lesion % change from baseline`))
table(Data$`Target lesion % change from baseline`, Data$BRR, useNA = "always")
#Data<-Data %>% mutate(Target.lesion.percent=Data$`Target lesion % change from baseline`*100) already in %

##----------------------------------------------------------------
#find shrinkage of tumors at indivudal target lesion. Already the Target lesion % change from baseline has been calculated

plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

plot.colors2<-c('CR' = '#0000CC','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#CC0000',
                'NA' = '#C0C0C0')
#'R' = '#CC0000', 'NR' = '#0000CC'

p1 = ggplot(Data, aes(x=reorder(PID, -`Target lesion % change from baseline`), y=`Target lesion % change from baseline`)) +
  geom_bar(aes(fill = `Best Response by RECIST`), color = '#000000', width = 0.6, stat = 'identity') +
  scale_fill_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90)+
  xlab("Patients")+
  geom_hline(yintercept = -30, linetype = 'dashed')
#scale_fill_discrete(drop=TRUE)

p1 

saveRDS(Data,"Data.rds")
Data<-readRDS("Data.rds")
pdf("Zane_and_Alex_Waterfall_Change_in_Target_Lesion_RECIST_9-21-23-update.pdf", width=50, height=15)
p1
dev.off()

my.levels<-c("PD","SD","PR","CR")
Data$`Best Response by RECIST`=factor(Data$`Best Response by RECIST`, levels=my.levels)

#Do the same but this time by NR vs R

p2 = ggplot(Data, aes(x=BRR, y=`Target lesion % change from baseline`)) +
  geom_boxplot(outlier.shape = NA) +
  theme_pubr(x.text.angle = 90)+
  xlab("Patients")+
  geom_hline(yintercept = -30, linetype = 'dashed')+
  geom_point(aes(fill = `Best Response by RECIST`), size = 5,  color = '#000000', shape = 21, position = "jitter")+
  scale_fill_manual(values = plot.colors) +
  stat_compare_means(comparisons = list(c('NR','R')),
                     method = 'wilcox.test')



p2

pdf("Zane_and_Alex_Boxplot_Change_in_Target_Lesion_RECIST_NR_v_R9-21-23-update.pdf")
p2
dev.off()


p3 = ggplot(Data, aes(x=`Best Response by RECIST`, y=`Target lesion % change from baseline`)) +
  geom_boxplot(outlier.shape = NA) +
  theme_pubr(x.text.angle = 90)+
  xlab("Patients")+
  geom_hline(yintercept = -30, linetype = 'dashed')+
  geom_point(aes(fill = BRR), size = 5, shape = 21, position = position_jitterdodge())


p3


pdf("Zane_and_Alex_Boxplot_Change_in_Target_Lesion_RECIST_all_RECIST_9-21-23-update.pdf")
p3
dev.off()


saveRDS(Data,"Data.rds")

my.site<-c("adrenal", "adrenal 2.adrenal.pct", "visceral other","visceral other 2", "peritoneum", "soft tissue",
           "soft tissue 2", "nodes", "nodes 2", "brain", "brain 2", "lung", "lung 2", "liver", "liver 2")
str(Data)
##convert all these to numeric
Data[, 45:74] <- lapply(45:74, function(x) as.numeric(Data[[x]]))
str(Data)
#for some reason, Sunny's loop wasn't working so I'm doing it by hand
adrenal.df<-Data[,c("adrenal (baseline)", "adrenal (best response)")]
adrenal.df<-as.data.frame(adrenal.df)
adrenal.df$subtract= adrenal.df[,2] - adrenal.df[,1]
adrenal.df$ratio = adrenal.df[,2] / adrenal.df[,1]
adrenal.df$pct = (adrenal.df[,2] / adrenal.df[,1] - 1) * 100
colnames(adrenal.df)[-c(1:2)] = paste0("adrenal",'.',colnames(adrenal.df)[-c(1:2)])
adrenal.df<-as.data.frame(adrenal.df)

Data = cbind(Data, adrenal.df[,-c(1:2)])
#adrenal 2
adrenal.2.df<-Data[,c("adrenal 2 (baseline)", "adrenal 2 (best response)")]
adrenal.2.df$subtract = adrenal.2.df[,2] - adrenal.2.df[,1]
adrenal.2.df$ratio = adrenal.2.df[,2] / adrenal.2.df[,1]
adrenal.2.df$pct = (adrenal.2.df[,2] / adrenal.2.df[,1] - 1) * 100
colnames(adrenal.2.df)[-c(1:2)] = paste0("adrenal 2",'.',colnames(adrenal.2.df)[-c(1:2)])
adrenal.2.df<-as.data.frame(adrenal.2.df)
Data = cbind(Data, adrenal.2.df[,-c(1:2)])

#visceral other
visceral.other.df<-Data[,c("visceral other (baseline)" , "visceral other (best response)")]
visceral.other.df$subtract = visceral.other.df[,2] - visceral.other.df[,1]
visceral.other.df$ratio = visceral.other.df[,2] / visceral.other.df[,1]
visceral.other.df$pct = (visceral.other.df[,2] / visceral.other.df[,1] - 1) * 100
colnames(visceral.other.df)[-c(1:2)] = paste0("visceral other",'.',colnames(visceral.other.df)[-c(1:2)])

Data = cbind(Data, visceral.other.df[,-c(1:2)])


#visceral other
visceral.other.2.df<-Data[,c("visceral other 2 (baseline)" , "visceral other 2 (best response)")]
visceral.other.2.df$subtract = visceral.other.2.df[,2] - visceral.other.2.df[,1]
visceral.other.2.df$ratio = visceral.other.2.df[,2] / visceral.other.2.df[,1]
visceral.other.2.df$pct = (visceral.other.2.df[,2] / visceral.other.2.df[,1] - 1) * 100
colnames(visceral.other.2.df)[-c(1:2)] = paste0("visceral other 2",'.',colnames(visceral.other.2.df)[-c(1:2)])

Data = cbind(Data, visceral.other.2.df[,-c(1:2)])

#peritoneum
peritoneum.df<-Data[,c("peritoneum (baseline)" , "peritoneum (best response)")]
peritoneum.df$subtract = peritoneum.df[,2] - peritoneum.df[,1]
peritoneum.df$ratio = peritoneum.df[,2] / peritoneum.df[,1]
peritoneum.df$pct = (peritoneum.df[,2] / peritoneum.df[,1] - 1) * 100
colnames(peritoneum.df)[-c(1:2)] = paste0("peritoneum",'.',colnames(peritoneum.df)[-c(1:2)])

Data = cbind(Data, peritoneum.df[,-c(1:2)])

#softTissue1

soft.tissue.df<-Data[,c("soft tissue (baseline)" , "soft tissue (best response)")]
soft.tissue.df$subtract = soft.tissue.df[,2] - soft.tissue.df[,1]
soft.tissue.df$ratio = soft.tissue.df[,2] / soft.tissue.df[,1]
soft.tissue.df$pct = (soft.tissue.df[,2] / soft.tissue.df[,1] - 1) * 100
colnames(soft.tissue.df)[-c(1:2)] = paste0("soft tissue",'.',colnames(soft.tissue.df)[-c(1:2)])

Data = cbind(Data, soft.tissue.df[,-c(1:2)])


#soft tissue 2
soft.tissue.2.df<-Data[,c("soft tissue 2 (baseline)" , "soft tissue 2 (best response)")]
soft.tissue.2.df$subtract = soft.tissue.2.df[,2] - soft.tissue.2.df[,1]
soft.tissue.2.df$ratio = soft.tissue.2.df[,2] / soft.tissue.2.df[,1]
soft.tissue.2.df$pct = (soft.tissue.2.df[,2] / soft.tissue.2.df[,1] - 1) * 100
colnames(soft.tissue.2.df)[-c(1:2)] = paste0("soft tissue 2",'.',colnames(soft.tissue.2.df)[-c(1:2)])

Data = cbind(Data, soft.tissue.2.df[,-c(1:2)])

#nodes
nodes.df<-Data[,c("nodes (baseline)" , "nodes (best response)")]
nodes.df$subtract = nodes.df[,2] - nodes.df[,1]
nodes.df$ratio = nodes.df[,2] / nodes.df[,1]
nodes.df$pct = (nodes.df[,2] / nodes.df[,1] - 1) * 100
colnames(nodes.df)[-c(1:2)] = paste0("nodes",'.',colnames(nodes.df)[-c(1:2)])

Data = cbind(Data, nodes.df[,-c(1:2)])


#nodes 2
nodes.2.df<-Data[,c("nodes 2 (baseline)" , "nodes 2 (best response)")]
nodes.2.df$subtract = nodes.2.df[,2] - nodes.2.df[,1]
nodes.2.df$ratio = nodes.2.df[,2] / nodes.2.df[,1]
nodes.2.df$pct = (nodes.2.df[,2] / nodes.2.df[,1] - 1) * 100
colnames(nodes.2.df)[-c(1:2)] = paste0("nodes 2",'.',colnames(nodes.2.df)[-c(1:2)])

Data = cbind(Data, nodes.2.df[,-c(1:2)])


#brain 

brain.df<-Data[,c("brain (baseline)" , "brain (best response)")]
brain.df$subtract = brain.df[,2] - brain.df[,1]
brain.df$ratio = brain.df[,2] / brain.df[,1]
brain.df$pct = (brain.df[,2] / brain.df[,1] - 1) * 100
colnames(brain.df)[-c(1:2)] = paste0("brain",'.',colnames(brain.df)[-c(1:2)])

Data = cbind(Data, brain.df[,-c(1:2)])


#brain 2
brain.2.df<-Data[,c("brain 2 (baseline)" , "brain 2 (best response)")]
brain.2.df$subtract = brain.2.df[,2] - brain.2.df[,1]
brain.2.df$ratio = brain.2.df[,2] / brain.2.df[,1]
brain.2.df$pct = (brain.2.df[,2] / brain.2.df[,1] - 1) * 100
colnames(brain.2.df)[-c(1:2)] = paste0("brain 2",'.',colnames(brain.2.df)[-c(1:2)])

Data = cbind(Data, brain.2.df[,-c(1:2)])

#lung

lung.df<-Data[,c("lung (baseline)" , "lung (best response)")]
lung.df$subtract = lung.df[,2] - lung.df[,1]
lung.df$ratio = lung.df[,2] / lung.df[,1]
lung.df$pct = (lung.df[,2] / lung.df[,1] - 1) * 100
colnames(lung.df)[-c(1:2)] = paste0("lung",'.',colnames(lung.df)[-c(1:2)])

Data = cbind(Data, lung.df[,-c(1:2)])



#lung 2

lung.2.df<-Data[,c("lung 2 (baseline)" , "lung 2 (best response)")]
lung.2.df$subtract = lung.2.df[,2] - lung.2.df[,1]
lung.2.df$ratio = lung.2.df[,2] / lung.2.df[,1]
lung.2.df$pct = (lung.2.df[,2] / lung.2.df[,1] - 1) * 100
colnames(lung.2.df)[-c(1:2)] = paste0("lung 2",'.',colnames(lung.2.df)[-c(1:2)])

Data = cbind(Data, lung.2.df[,-c(1:2)])

#liver
liver.df<-Data[,c("liver (baseline)" , "liver (best response)")]
liver.df$subtract = liver.df[,2] - liver.df[,1]
liver.df$ratio = liver.df[,2] / liver.df[,1]
liver.df$pct = (liver.df[,2] / liver.df[,1] - 1) * 100
colnames(liver.df)[-c(1:2)] = paste0("liver",'.',colnames(liver.df)[-c(1:2)])

Data = cbind(Data, liver.df[,-c(1:2)])

#liver 2

liver.2.df<-Data[,c("liver 2 (baseline)" , "liver 2 (best response)")]
liver.2.df$`liver 2 (best response)`<-as.numeric(liver.2.df$`liver 2 (best response)`)
liver.2.df$subtract = liver.2.df[,2] - liver.2.df[,1]
liver.2.df$ratio = liver.2.df[,2] / liver.2.df[,1]
liver.2.df$pct = (liver.2.df[,2] / liver.2.df[,1] - 1) * 100
colnames(liver.2.df)[-c(1:2)] = paste0("liver 2",'.',colnames(liver.2.df)[-c(1:2)])

Data = cbind(Data, liver.2.df[,-c(1:2)])




saveRDS(Data,"Data.rds") #Data has all the change in sizes among tissues
Data<-readRDS("Data.rds")
percent<-c( "adrenal.pct", "adrenal 2.pct", "visceral other.pct", "visceral other 2.pct", "peritoneum.pct", "soft tissue.pct","soft tissue 2.pct","nodes.pct","nodes 2.pct","brain.pct","brain 2.pct","lung.pct", "lung 2.pct","liver.pct","liver 2.pct" )     
Data1 = Data[,c("PID","Best Response by RECIST","BRR",
                percent)]


Data1 = reshape2::melt(Data1)
colnames(Data1)[(ncol(Data1)-1):ncol(Data1)] = c('Site','Change.pct')
saveRDS(Data1,"Data1")


## lesion 1 or 2
Data1$Lesion.number = "LS1"
Data1$Lesion.number[grep('2[.]pct$', Data1$Site)] = 'LS2'
table(Data1$Lesion.number)
# LS1  LS2 
#1744 1526

saveRDS(Data1,"Data1.rds")
## use line to connect lesions from the same pt instead of calling it site and site 2
#Data$Site = gsub('[.]2[.]pct$','.pct',Data$Site)
#print(table(Data$Site))

Data1$Group = paste0(Data1$BRR,'.',Data1$Lesion.number)


saveRDS(Data1,"Data1.rds")

##-------------------------------------------------
## for each site, take mean of multiple lesions 

#for some reason adrenal 2 has a weird name so change it to proper
#Data1<- lapply(Data1, function(x){gsub("adrenal 2.adrenal.pct", "adrenal 2.pct", x)})
Data1<-as.data.frame(Data1)
Data1$Lesion.number[grep('2[.]pct$', Data1$Site)] = 'LS2'
#need to fix data plot 2 as it isn't taking average across multiple lesions
Data1<-Data1 %>% group_by(PID)
saveRDS(Data1,"Data1.rds")
#i think the prolem is that theree are two sites. Change all sites to the same and rely on lesion for future separation
Data2<-data.frame(Data1)
Data2[Data2=="adrenal 2.pct"]<-"adrenal.pct"
Data2[Data2=="visceral other 2.pct"]<-"visceral other.pct" 
Data2[Data2=="soft tissue 2.pct"]<-"soft tissue.pct" 
Data2[Data2=="nodes 2.pct"]<-"nodes.pct" 
Data2[Data2=="brain 2.pct"]<-"brain.pct" 
Data2[Data2=="lung 2.pct"]<-"lung.pct" 
Data2[Data2=="liver 2.pct"]<-"liver.pct" 
saveRDS(Data2, "Data.2.rds")

Data2<-Data2 %>% group_by(Site)
Data2<-Data2[!is.na(Data2$Change.pct),] 
str(Data2)
Data2$Change.pct<-as.numeric(Data2$Change.pct)
#note two ways to perform mean expression and I verified each works (a and b below)
b<-Data2 %>% group_by(PID,Site) %>%
  summarise(mean = mean(Change.pct), na.rm=TRUE)
b

a<-aggregate(Data2$Change.pct, by=list(PID=Data2$PID, Site=Data2$Site),data=Data2,FUN=mean, na.rm=TRUE)


b1<-merge(Data2, b)
saveRDS(b1, "average_percent_change_by_patient_and_Site.rds")
b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
#b1 or average_percent_change_by_patient_and_Site is average while Data2 is individual
dim(b1) #[1][1] 236   9
my.levels<-c("PD","SD","PR","CR")
#b1<-b1 %>% rename(Best.Response.by.RECIST=Best.response.by.RECIST)
b1$Best.Response.by.RECIST<-factor(b1$Best.Response.by.RECIST, my.levels)
b1$Best.Response.by.RECIST<-na.omit(b1$Best.Response.by.RECIST)

#editing this because Zane wants the facet wrap not to say ".pct"
b1<-b1 %>% mutate(Site.renamed=b1$Site)
b1$Site.renamed<-str_remove(b1$Site.renamed, ".pct")
b2<-b1 %>% rename(`Lesion Number`=Lesion.number)
b2<-b2 %>% filter(!Site.renamed %in% c("visceral other", "peritoneum"))
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

plot.colors2 = c('CR' = '#0000CC','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' ='#CC0000' ,
                 'NA' = '#C0C0C0')

rm(p1, p2)
rm(p1, p2)
my_comparisons <- list( c("NR.LS1", "NR.LS2"), c("R.LS1", "R.LS2") )
p1 = ggplot(b2, aes(Group,
                    Change.pct)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_point(aes(color = Best.Response.by.RECIST,
                 shape = `Lesion Number`)) +
  geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  stat_compare_means(comparisons = my_comparisons, method="wilcox")

p1


pdf("Zane_and_Alex_boxplot_Change_percent_between_lesion_and_Group_9-21-23_update_colors.pdf", width=8)
p1
dev.off()

#edited for Zane's talk 9-26-22
b2<-b1 %>% rename(`Best Response by RECIST`=Best.Response.by.RECIST)
#edited for Zane's presentation
#•	Left and right figure, remove adrenal, visceral, peritoneum
b2<-b2 %>% rename(`Lesion Number`=Lesion.number)
b2<-b2 %>% filter(!Site.renamed %in% c("visceral other", "peritoneum"))
b2$`Best Response by RECIST`<-factor(b2$`Best Response by RECIST`, levels=c("CR","PR","SD","PD"))
p1 = ggplot(b2, aes(Group,
                    Change.pct)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_point(aes(color = `Best Response by RECIST`,
                 shape = `Lesion Number`)) +
  geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  stat_compare_means(comparisons = my_comparisons, method="wilcox")+
  ylab("Lesion Response (%)")

p1


pdf("Zane_and_Alex_combo_Edited_9-21-23_boxplot_Change_percent_between_lesion_and_Group_colors.pdf", width = 15)
p1
dev.off()


p2 = ggplot(b2, aes(BRR,
                    mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = `Best Response by RECIST`), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  stat_compare_means(comparisons = list(c('NR','R'))) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  xlab("Overall Best Response by RECIST")+
  ylab("Lesion Response (%)")

p2


pdf("Zane_and_Alex_combo_edited_boxplot_update_9-21-23_Change_percent_between_only_group_colors.pdf")
p2
dev.off()


a<-b2 %>% group_by(BRR, Site.renamed) %>% summarize (median(mean), total_count=n(), sd=sd(mean))
a<-as.data.frame(a)

write.xlsx(a, "Median_and_count.by.RECIST_location_with_sd.xlsx")


#separate into 4 groups
p2.2 = ggplot(b1, aes(Site.renamed,
                      mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  scale_color_manual(values = plot.colors) +
  stat_compare_means(comparisons = list(c('NR','R'))) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ BRR)

p2.2

pdf("Zane_and_Alex_boxplot_Change_percent_site_and_Site_NR_v_R_9-21-23_update.pdf")
p2.2
dev.off()

#Zane wanted median values for each Site so summarize that here
saveRDS(b1, "average_percent_change_by_patient_and_Site.rds")
d<-b1 %>% group_by(BRR, Site.renamed) %>% summarize(Median=median(Change.pct))
write.csv(d,"Median_for_Zane_and_Alex_boxplot_Change_percent_site_and_Site_NR_v_R_9-21-23_update.csv")

#mean values for eosinphil count and neutorphil:leukocyte ratio
#% eos pre-RX

Data<-Data %>%mutate(Eosinophil_dec=(as.numeric(Data$'% eos pre-RX'))/100)
Data$'WBC count pre-Rx (x10^3)'<-as.numeric(Data$'WBC count pre-Rx (x10^3)')
Data<-Data %>% mutate(Eosinophil_count=(Data$Eosinophil_dec*Data$`WBC count pre-Rx (x10^3)`))
#keep<-c("% eos pre-RX", "Eosinophil_dec","Eosinophil_count", "WBC count pre-Rx (x10^3)")
#d.1<-Data[,keep]
Data<-Data %>% rename("Eosinophil_count(x10^3)"="Eosinophil_count")
#to.remove<-c(78, 24)

#Data.3<-Data[-to.remove,]
e<-mean(Data$`Eosinophil_count(x10^3)`, na.rm=TRUE) #[1] 0.1810165

#neutrophil:leukocyte ratio 
saveRDS(Data, "Data.rds")
#%neutr pre-Rx/% lymph pre-Rx
Data$`%neutr pre-Rx`<-as.numeric(Data$`%neutr pre-Rx`)
Data$`% lymph pre-Rx`<-as.numeric(Data$`% lymph pre-Rx`)
Data<-Data %>% mutate(NLR=(Data$`%neutr pre-Rx`/Data$`% lymph pre-Rx`))
keep<-c("PID","%neutr pre-Rx", "% lymph pre-Rx","NLR")
d.2<-Data[,keep]

write.csv(d.2, "Zane_and_Alex_Neutrophil_to_Lymphocyte_9-21-23_update.csv")

b1$Best.Response.by.RECIST=factor(b1$Best.Response.by.RECIST, levels=c("CR","PR","SD","PD"))
p3 = ggplot(b1, aes(Site.renamed,
                    mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  scale_color_manual(values = plot.colors) +
  stat_compare_means(comparisons = list(c('NR','R'))) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Best.Response.by.RECIST)+
  geom_line(aes(group = PID))

p3

pdf("Zane_and_Alex_boxplot_9-21-23_Change_percent_site_and_Site_BRR.pdf")
p3
dev.off()
















b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
#Direction plots
b1$change.direction = NA
b1$change.direction[which(b1$Change.pct>0)] = 'up'
b1$change.direction[which(b1$Change.pct<0)] = 'down'
b1$change.direction[which(b1$Change.pct==0)] = 'nochange'
head(b1)

table(b1$change.direction)
dim(b1) #236, 10
sum(is.na(b1$Change.pct)) #0
length(unique(b1$Site)) #8 sites



data.plot.2.2 = NULL 
for(my.response in c('PD','SD','PR','CR')) {
  print(my.response)
  
  my.df = NULL 
  my.df = data.frame(table(
    b1[b1$Best.Response.by.RECIST==my.response,
       c('PID',
         'change.direction'),drop=F]))
  my.df = reshape2::dcast(PID ~ change.direction,
                          data = my.df,
                          value.var = 'Freq',
                          fun.aggregate = sum)
  my.df$Best.Response.by.RECIST = my.response
  rm(x)
  x = my.df[,c('PID','Best.Response.by.RECIST')]
  for( my.col in c('down','nochange','up')) {
    print(my.col)
    x$new = NA
    if(length(which(colnames(my.df)==my.col))>0) {
      x$new = my.df[,my.col]
      
    }
    colnames(x)[which(colnames(x)=='new')] = my.col
  }
  my.df = x
  data.plot.2.2 = rbind(data.plot.2.2,
                        my.df)
}
data.plot.2.2[is.na(data.plot.2.2)] = 0

## bar graphs
data.plot.2.2.melt = NULL 
data.plot.2.2.melt = reshape2::melt(data.plot.2.2)
colnames(data.plot.2.2.melt) = c('PID',
                                 'Best.Response.by.RECIST',
                                 'change.direction',
                                 'sites')
data.plot.2.2.melt$Best.Response.by.RECIST = 
  factor(data.plot.2.2.melt$Best.Response.by.RECIST,
         levels = c('PD','SD','PR','CR'))
p4.2 = NULL 
p4.2 = ggplot(data.plot.2.2.melt,aes(PID,sites)) +
  geom_bar(aes(fill = change.direction), width = 0.6,
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values = c('up' = '#0000CC',
                               'nochange' = '#C0C0C0',
                               'down' = '#CC0000')) + 
  theme_pubr(x.text.angle = 90) +
  facet_grid(~ Best.Response.by.RECIST, scale = 'free_x', space = 'free_x')
p4.2


pdf("Zane_and_Alex_change_direction_by_patient.9-21-23_update.pdf")
p4.2
dev.off()


## tables 
data.plot.2.3 = data.plot.2.2
data.plot.2.3$total = apply(data.plot.2.3[,-c(1:2)], 1, sum)
data.plot.2.3[,3:5] = data.plot.2.3[,3:5] / data.plot.2.3[,6]
data.plot.2.3$category = NA
data.plot.2.3$category[which(data.plot.2.3$down==0 &
                               data.plot.2.3$nochange==0 &
                               data.plot.2.3$up == 1)] = 'uniform progression' #change 9-26 to uniform progression
data.plot.2.3$category[which(data.plot.2.3$down==1 &
                               data.plot.2.3$nochange==0 &
                               data.plot.2.3$up == 0)] = 'uniform regression'#change 9-26 to uniform regression
data.plot.2.3$category[which(data.plot.2.3$down==0 &
                               data.plot.2.3$nochange==1 &
                               data.plot.2.3$up == 0)] = 'all.nochange'
data.plot.2.3$category[is.na(data.plot.2.3$category)] = 'mixed response'

table(data.plot.2.3$category)
#  mixed response uniform progression  uniform regression 
###20                  30                  32 


data.plot.2.3 = data.frame(
  table(data.plot.2.3[,c('Best.Response.by.RECIST','category')]),
  stringsAsFactors = F)
data.plot.2.3$Best.Response.by.RECIST = factor(
  data.plot.2.3$Best.Response.by.RECIST,
  levels = c('PD','SD','PR','CR')
)
p4.3 = NULL 
p4.3 = ggplot(data.plot.2.3, aes(Best.Response.by.RECIST,
                                 Freq)) +
  geom_bar(aes(fill = category),width = 0.6,
           position = 'stack', stat = 'identity',
           color = '#000000') +
  theme_pubr() +
  scale_fill_manual(values = c('uniform regression' = '#0000CC',
                               #'all.nochange' = '#C0C0C0',
                               'mixed response' = '#00CC00',
                               'uniform progression' = '#CC0000')) 
p4.3





pdf("Zane_and_Alex_change_direction_by_response_all_9-21-23.pdf")
p4.3
dev.off()


write.csv(data.plot.2.3, "Change.direction.by.response.n.group_9-21-23_all.csv")



#Zane wants to know heterogenity betwen % change in each lesion
#try abs(l1-l2)/(mean(l1+l2)))
#and try abs(l1-l1))

#probably a simple way to do this but cna't figure out. going to export to excel and do it there
write.csv(b1, "average_percent_change_by_patient_and_Site.csv")

B1.1<-b1 %>% filter(BRR=="NR")
B1.1<-B1.1 %>% group_by(PID, Site) %>% filter(n()>1)

B1.1<-B1.1 %>% select (PID, Site, Change.pct, mean, Lesion.number)
B1.1<-B1.1 %>% group_by(PID) %>% mutate(delta.over.mean="N/A")
B1.1<-B1.1[,-6]
B1.2<-inner_join(filter(B1.1, Lesion.number=="LS1"), filter(B1.1, Lesion.number=="LS2"), by=c("PID","Site")) %>%
  mutate(
    homog.over.mean=abs(Change.pct.y-Change.pct.x)/(mean(Change.pct.y+Change.pct.x)),
    homog.only=abs(Change.pct.y-Change.pct.x))


saveRDS(B1.2, "average_percent_change_by_patient_and_site_with_heterogenous_change.rds")

p1<-ggplot(B1.2, aes(x=Site, y=homog.over.mean))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  ylab("Difference over mean in percent change Lesion 1 and Lesion 2 per tissue")+
  stat_compare_means(method="wilcox.test", ref.group = "liver.pct", label = "p.signif", label.y = 1.6)                                       


pdf("Zane_and_Alex_Difference_in_percent_change_between_lesions_Over_Mean_9-21-23_update.pdf")
p1
dev.off()


p2<-ggplot(B1.2, aes(x=Site, y=homog.only))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  ylab("Difference in percent change Lesion 1 and Lesion 2 per tissue")+
  stat_compare_means(method="wilcox.test", ref.group = "liver.pct", label = "p.signif") 
p2

pdf("Zane_and_Alex_Difference_in_percent_change_between_lesions_9-21-23-update.pdf")
p2
dev.off()


##Re-run all plots for split in clinical factors
#2-7-23
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
##Mutational status#################################################################
##start with mutational status (BRAF, NF1, NRAS)---------------------------------
#dataplot has all clinical, b1 has averaged change.pct
Data<-readRDS("Data.rds")
mutation<- Data %>% select("PID","BRAF status", "RAS (NRAS/KRAS)", "NF1", "Other (specify)")

Mutational.df<-merge(b1, mutation)
saveRDS(Mutational.df, "Mutational.df.contains.all.replicates.rds")
Mutational.df.1<-filter(Mutational.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(Mutational.df.1, "Mutational.df.contains.only.mean.rds")
Mutational.df.1<-readRDS("Mutational.df.contains.only.mean.rds")
table(Mutational.df.1$`BRAF status`)
table(Mutational.df.1$`RAS (NRAS/KRAS)`)
table(Mutational.df.1$NF1)
table(Mutational.df.1$`Other (specify)`)

Mutational.df.1$Mutational_group<-paste0(Mutational.df.1$`BRAF status`, "_",
                                         Mutational.df.1$`RAS (NRAS/KRAS)`, "_",
                                         Mutational.df.1$NF1)

table(Mutational.df.1$Mutational_group)
Mutational.df.1<-Mutational.df.1 %>% filter(Mutational_group!="WT_NRAS_NF1")
table(Mutational.df.1$Mutational_group)



WT.NA<-c("N/A_N/A_N/A", "NA_N/A_N/A","WT_N/A_N/A","WT_WT_N/A","WT_WT_WT","WT_WT_NA")
BRAF<-c("V600_N/A_N/A",  "V600_WT_N/A",   "V600_WT_WT","V600_N/A_WT")
NRAS<-c("WT_NRAS_N/A",   "WT_NRAS_WT")
NF1<-c("WT_WT_NF1")


Mutational.df.1$Simplified_mutational_Group=NA
Mutational.df.1$Simplified_mutational_Group[which(Mutational.df.1$Mutational_group %in% WT.NA)] = 'NAorWT' 
Mutational.df.1$Simplified_mutational_Group[which(Mutational.df.1$Mutational_group %in% BRAF)] = 'BRAF' 
Mutational.df.1$Simplified_mutational_Group[which(Mutational.df.1$Mutational_group %in% NRAS)] = 'NRAS' 
Mutational.df.1$Simplified_mutational_Group[which(Mutational.df.1$Mutational_group %in% NF1)] = 'NF1' 
table(Mutational.df.1$Simplified_mutational_Group, Mutational.df.1$PID)
Mutational.df.1$Site.renamed<-"N/A"
Mutational.df.1<-mutate(Mutational.df.1, Site.renamed=str_remove(Mutational.df.1$Site, ".pct"))
#Mutational.df.1<-Mutational.df.1 %>% rename(Best.Response.by.RECIST=Best.response.by.RECIST)

saveRDS(Mutational.df.1, "Mutational.df.contains.only.mean.rds")
Mutational.df.1<-readRDS("Mutational.df.contains.only.mean.rds")

to.remove<-c("peritoneum", "visceral other")
Mutational.df.1<-Mutational.df.1 %>% filter(!Site.renamed %in% to.remove)


p10 = ggplot(Mutational.df.1, aes(Simplified_mutational_Group,
                                  mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  stat_compare_means(comparisons = list(c('NR','R'))) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Mutational Status")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p10

pdf("Zane_and_Alex_edited_9-21-23_Mutational_status_vs_response.pdf")
p10
dev.off()

#perform pairwise wilcox test (need to remove one sample since it only has one of the mutational groups (need more than 2 to compare))

a<-Mutational.df.1 %>% group_by(Site.renamed, BRR) %>% count(Simplified_mutational_Group)
a.1<-a %>% count(Site.renamed)
a.2<-merge(Mutational.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one LDH per group

#pairwise wilcox is giving me so much trouble. I'm going to try to break it down to see which samples are the issue
nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Simplified_mutational_Group, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R")
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Simplified_mutational_Group, p.adjust.method="bonf")
#& Site.renamed %in% c("brain", "lung","nodes","soft tissue", "visceral other"

b<-rbind(nr.1, respond.1) 

write.csv(b, "Zane_and_Alex_Pairwise_edited_9-21-23-wilcox_test_for_mutational_status_vs_response.csv")
saveRDS(Data,"Data.rds")
saveRDS(b1,"average_percent_change_by_patient_and_Site.rds")
######################################################################
##LDH elevated vs low---------------------------------
#Text message from Zane 9-12-22, LDH 171 and above is elevated, below is within normal limit (WNL)
Data<-readRDS("Data.rds")
#Data<-Data %>% rename(Best.Response.by.RECIST=Best.response.by.RECIST)
LDH<- Data %>% select("PID","LDH pre-rx")

LDH$LDH.levels[LDH$`LDH pre-rx` >170] <- "elevated"
LDH$LDH.levels[LDH$`LDH pre-rx` =="N/A"] <- "N/A"
#for some reason there is one that just has "n" so replace that with N/A
LDH$LDH.levels[LDH$`LDH pre-rx` =="n"] <- "N/A"
LDH$LDH.levels[LDH$`LDH pre-rx` <171] <- "WNL"

b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
b1<-b1 %>% mutate(Site.renamed=b1$Site)
b1$Site.renamed<-str_remove(b1$Site.renamed, ".pct")
#b1<-b1 %>% rename(Best.Response.by.RECIST=Best.response.by.RECIST)
saveRDS(b1,"average_percent_change_by_patient_and_Site.rds")
View(b1)

LDH.df<-merge(b1, LDH)
saveRDS(LDH.df, "LDH.df.contains.all.replicates.rds")
LDH.df<-readRDS("LDH.df.contains.all.replicates.rds")
LDH.df.1<-filter(LDH.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(LDH.df.1, "LDH.df.contains.only.mean.rds")
LDH.df.1<-readRDS("LDH.df.contains.only.mean.rds")
to.remove<-c("peritoneum", "visceral other")
LDH.df.1<-LDH.df.1 %>% filter(!Site.renamed %in% to.remove)
p11 = ggplot(LDH.df.1, aes(LDH.levels,
                           mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  stat_compare_means(comparisons = list(c('NR','R'))) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("LDH levels")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p11

pdf("Zane_and_Alex_9-21-23_update_LDH_vs_response.pdf")
p11
dev.off()

a<-LDH.df.1 %>% group_by(Site.renamed, BRR) %>% count(LDH.levels)
a.1<-a %>% count(Site.renamed)
a.2<-merge(LDH.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one LDH per group

#pairwise wilcox is giving me so much trouble. I'm going to try to break it down to see which samples are the issue
nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ LDH.levels, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R" )
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ LDH.levels, p.adjust.method="bonf")
#& Site.renamed %in% c("lung","nodes","soft tissue")

b<-rbind(nr.1, respond.1) 


write.csv(b, "Zane_and_Alex_9-21-23_Pairwise_wilcox_test_for_LDH_vs_response.csv")
saveRDS(Data, "Data.rds")

######################################################################
##BMI-----------------------------------------------------------------
Data<-Data %>% mutate(`height (m)`=`height (cm)`*0.01)
Data<-Data %>% mutate(BMI=`weight (Kg)`/(`height (m)`)^2)
saveRDS(Data,"Data.rds")

BMI<- Data %>% select("PID","BMI")
BMI$BMI.25[BMI$BMI >25] <- "high"
BMI$BMI.25[BMI$BMI <=25] <- "low"


BMI$BMI.30[BMI$BMI >30] <- "high"
BMI$BMI.30[BMI$BMI <=30] <- "low"

table(Data$BMI)
b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
BMI.df<-merge(b1, BMI)

saveRDS(BMI.df, "BMI.df.contains.all.replicates.rds")
BMI.df<-readRDS("BMI.df.contains.all.replicates.rds")
BMI.df.1<-filter(BMI.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(BMI.df.1, "BMI.df.contains.only.mean.rds")
BMI.df.1<-readRDS("BMI.df.contains.only.mean.rds")

to.remove<-c("peritoneum", "visceral other")
BMI.df.1<-BMI.df.1 %>% filter(!Site.renamed %in% to.remove)


#BMI>25
p12 = ggplot(BMI.df.1, aes(BMI.25,
                           mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("BMI>25")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p12

pdf("Zane_and_Alex_9-21-23_update_BMI_over_25_vs_response.pdf")
p12
dev.off()

#calculate pairwise significance levels
a<-BMI.df.1 %>% group_by(Site.renamed, BRR) %>% count(BMI.25) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(BMI.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one LDH per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ BMI.25, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_BMI_greater_than_25_vs_response.csv")




#BMI>30
p13 = ggplot(BMI.df.1, aes(BMI.30,
                           mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("BMI>30")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p13

pdf("Zane_and_Alex_9-21-23_update_BMI_over_30_vs_response.pdf")
p13
dev.off()


a<-BMI.df.1 %>% group_by(Site.renamed, BRR) %>% count(BMI.30) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(BMI.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one LDH per group


nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ BMI.30, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R")
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ BMI.30, p.adjust.method="bonf")
#& Site.renamed %in% c("adrenal", "liver","lung","nodes","soft tissue", "visceral other"

b<-rbind(nr.1, respond.1) 

write.csv(b, "Zane_and_Alex_9-21-23_Pairwise_wilcox_test_for_BMI_greater_than_30_vs_response.csv")

#make density plot of BMI
d <- density(BMI.df.1$BMI) # returns the density data 
plot(d)

pdf("Zane_and_Alex_edited_9-21-23_BMI_density.pdf")
plot(d)
dev.off()

length(unique(BMI.df.1$PID)) #number of unique patients(not including multiple tissues per patient)
unique.patient.df<- BMI.df.1[!duplicated(BMI.df.1$PID),]

table(unique.patient.df$BMI.25)
table(unique.patient.df$BMI.30)

#skipping Metabolic syndrome question-------------------------------

#Albumin levels Albumin <3.5
Data<-readRDS("Data.rds")
Data$`Albumin pre-rx`<-as.numeric(Data$`Albumin pre-rx`)
Data<-Data %>% mutate(Albumin.level=ifelse(Data$`Albumin pre-rx`< 3.5, "low", "WNL"))
saveRDS(Data, "Data.rds")
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

Albumin<- Data %>% select("PID","Albumin.level")


b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
Albumin.df<-merge(b1, Albumin)

saveRDS(Albumin.df, "Albumin.df.contains.all.replicates.rds")
Albumin.df.1<-filter(Albumin.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(Albumin.df.1, "Albumin.df.contains.only.mean.rds")
Albumin.df.1<-readRDS("Albumin.df.contains.only.mean.rds")

#Albumin
p14 = ggplot(Albumin.df.1, aes(Albumin.level,
                               mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Albumin Level, low is <3.5")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p14

pdf("Zane_and_Alex_9-21-23_update_Albumin_vs_response.pdf")
p14
dev.off()


a<-Albumin.df.1 %>% group_by(Site.renamed, BRR) %>% count(Albumin.level) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Albumin.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one LDH per group

nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Albumin.level, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R")
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Albumin.level, p.adjust.method="bonf")
#& Site.renamed %in% c("adrenal", "brain","lung","nodes","soft tissue", "visceral other")

b<-rbind(nr.1, respond.1) 



write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_Albumin_vs_response.csv")
#############################################################
#melanoma subtype
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
Data<-readRDS("Data.rds")
table(Data$`Melanoma subtype`)
Data<-Data %>% mutate(Melanoma_subtype_simplified=ifelse(Data$`Melanoma subtype`=="Cutaneous", "Cutaneous", "Others"))
table(Data$Melanoma_subtype_simplified, Data$`Melanoma subtype`)
#Data$Melanoma_subtype_simplified[Data$`Melanoma subtype`=="Cutaneous/Acral"]<-"Cutaneous"
saveRDS(Data, "Data.rds")

Subtype<- Data %>% select("PID","Melanoma_subtype_simplified")


b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
subtype.df<-merge(b1, Subtype)

#subtype.df<-subtype.df %>% filter (`Melanoma subtype` !="unknown primary")
table(subtype.df$Melanoma_subtype_simplified)

saveRDS(subtype.df, "subtype.df.contains.all.replicates.rds")
subtype.df.1<-filter(subtype.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(subtype.df.1, "subtype.df.contains.only.mean.rds")
subtype.df.1<-readRDS("subtype.df.contains.only.mean.rds")
#Melanoma subtype
p15 = ggplot(subtype.df.1, aes(Melanoma_subtype_simplified,
                               mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Melanoma Subtype")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p15

pdf("Zane_and_Alex_edited_9-21-23_Melanoma_subtype_vs_response.pdf")
p15
dev.off()



a<-subtype.df.1 %>% group_by(Site.renamed, BRR) %>% count(Melanoma_subtype_simplified) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(subtype.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1)
#a.2<- a.2 %>% filter(Melanoma.subtype !="unknown primary")


#remove out samples that only have one  per group
#pairwise wilcox is giving me so much trouble. I'm going to try to break it down to see which samples are the issue
nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Melanoma_subtype_simplified, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R")
respond<-respond %>% filter(n>1) 
respond <-respond %>% filter(Site.renamed !="peritoneum")
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Melanoma_subtype_simplified, p.adjust.method="bonf")


b<-rbind(nr.1, respond.1) #note, had to do it differently wheere we only look at nodes since that is the only 
#tissue that has R with more than just one sample in comparison group.

write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_melanoma_subtype_vs_response.csv")
saveRDS(Data, "Data.rds")
####################--------------------------------------------------
#GENDER
Data<-readRDS("Data.rds")
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
table(Data$Sex)
#Female   Male 
#85    133
Gender<- Data %>% select("PID","Sex")


b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
gender.df<-merge(b1, Gender)


table(gender.df$Sex)

saveRDS(gender.df, "gender.df.contains.all.replicates.rds")
gender.df.1<-filter(gender.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(gender.df.1, "gender.df.contains.only.mean.rds")

#Melanoma subtype
p17 = ggplot(gender.df.1, aes(Sex,
                              mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Sex")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p17

pdf("Zane_and_Alex_9-21-23_update_Melanoma_Gender_vs_response.pdf")
p17
dev.off()

a<-gender.df.1 %>% group_by(Site.renamed, BRR) %>% count(Sex) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(gender.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ Sex, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_gender_vs_response.csv")
####################--------------------------------
#Baseline tumor size, assuming this is in cm

plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

Data<-readRDS("Data.rds")
table(Data$`Baseline lesion measurement`)

Tumor.size<- Data %>% select("PID","Baseline lesion measurement")
Tumor.size$`Baseline lesion measurement`<-as.numeric(Tumor.size$`Baseline lesion measurement`)


Tumor.size$size[Tumor.size$`Baseline lesion measurement`>5]<-"large"
Tumor.size$size[Tumor.size$`Baseline lesion measurement`<=5]<-"small"

table(Tumor.size$size, useNA="always")

b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
baseline.size.df<-merge(b1, Tumor.size)

table(baseline.size.df$`Baseline lesion measurement`)

saveRDS(baseline.size.df, "baseline.size.df.contains.all.replicates.rds")
baseline.size.df.1<-filter(baseline.size.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(baseline.size.df.1, "baseline.size.df.contains.only.mean.rds")

#size
p18 = ggplot(baseline.size.df.1, aes(size,
                                     mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Baseline Tumor size")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p18

pdf("Zane_and_Alex_9-21-23_update_Melanoma_Baseline_tumor_size_vs_response.pdf")
p18
dev.off()

a<-baseline.size.df.1 %>% group_by(Site.renamed, BRR) %>% count(size) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(baseline.size.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ size, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_baseline_tumor_size_vs_response.csv")

table(baseline.size.df.1$size)
####################---------------------------------------------------

#total WBC
#Baseline tumor size, assuming this is in cm
Data<-readRDS("Data.rds")
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
table(Data$`WBC count pre-Rx (x10^3)`)
str(Data$`WBC count pre-Rx (x10^3)`)
WBC.count<- Data %>% select("PID","WBC count pre-Rx (x10^3)")

#4.5-11 is normal range 
WBC.count$WBC_group[WBC.count$`WBC count pre-Rx (x10^3)`<4.5]<-"low"
WBC.count$WBC_group[WBC.count$`WBC count pre-Rx (x10^3)`>11]<-"high"
WBC.count$WBC_group[WBC.count$`WBC count pre-Rx (x10^3)`>=4.5 &WBC.count$`WBC count pre-Rx (x10^3)`<=11 ]<-"normal"


table(WBC.count$WBC_group, useNA="always") #4 #NAs
b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
WBC.count.df<-merge(b1, WBC.count)

table(WBC.count.df$WBC_group)

saveRDS(WBC.count.df, "WBC.count.df.contains.all.replicates.rds")
WBC.count.df.1<-filter(WBC.count.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(WBC.count.df.1, "WBC.count.df.contains.only.mean.rds")

#WBC
p19 = ggplot(WBC.count.df.1, aes(WBC_group,
                                 mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("WBC Group, <4.5 is low, >11 is high")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p19

pdf("Zane_and_Alex_9-21-23_update_Melanoma_WBC_group_vs_response.pdf")
p19
dev.off()

a<-WBC.count.df.1 %>% group_by(Site.renamed, BRR) %>% count(WBC_group) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(WBC.count.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ WBC_group, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R")
respond<-respond %>% filter(n>1) 
respond <-respond %>% filter(!Site.renamed  %in% c("peritoneum","visceral other"))
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ WBC_group, p.adjust.method="bonf")


b<-rbind(nr.1, respond.1) #note, had to do it differently wheere we only look at nodes since that is the only 
write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_WBC_group_vs_response.csv")

####################--------------------------------
#neutrophil lymphocyte ratio (from above)
Data<-readRDS("Data.rds")
Data$`%neutr pre-Rx`<-as.numeric(Data$`%neutr pre-Rx`)
Data$`% lymph pre-Rx`<-as.numeric(Data$`% lymph pre-Rx`)

Data<-Data %>% mutate (NLR=Data$`%neutr pre-Rx`/Data$`% lymph pre-Rx`)
Data$NLR_group[Data$NLR>3]<-"high"
Data$NLR_group[Data$NLR<=3]<-"normal"
b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
NLR.table<-merge(b1, Data, by="PID")
saveRDS(Data, "Data.rds")

plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
table(NLR.table$NLR)

NLR<- NLR.table %>% select("PID","NLR_group")

table(NLR$NLR_group, useNA="always")

NLR.df<-merge(b1, NLR)

table(NLR.df$NLR_group) #high normal 
###########################1448    692 

saveRDS(NLR.df, "NLR.df.contains.all.replicates.rds")
NLR.df.1<-filter(NLR.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(NLR.df.1, "NLR.df.contains.only.mean.rds")
NLR.df.1<-readRDS("NLR.df.contains.only.mean.rds")
NLR.df.1<-NLR.df.1 %>% filter(!Site.renamed %in% c( "visceral other", "peritoneum"))

#NLR
p20 = ggplot(NLR.df.1, aes(NLR_group,
                           mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Neutrophil:Lymphocyte Ratio")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p20

pdf("Zane_and_Alex_9-21-23_update_Melanoma_NLR_vs_response_colors.pdf")
p20
dev.off()

a<-NLR.df.1 %>% group_by(Site.renamed, BRR) %>% count(NLR_group) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(NLR.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group


nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ NLR_group, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R")
respond<-respond %>% filter(n>1) 
respond <-respond %>% filter(!Site.renamed  %in% c("peritoneum","visceral other"))
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ NLR_group, p.adjust.method="bonf")

b<-rbind(nr.1, respond.1)


write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_NLR_vs_response.csv")
###--------------------------------------------------------------------
##NLR split by median
#neutrophil lymphocyte ratio (from above)
Data<-readRDS("Data.rds")
Data$NLR<-as.numeric(Data$NLR)
a<-Data$NLR
a<-na.omit(a)
median(a)

Data$NLR_group_median[Data$NLR>median(a)]<-"high"
Data$NLR_group_median[Data$NLR<=median(a)]<-"low"
head(Data)

NLR_median<-Data %>% select("PID", "NLR_group_median")

b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
NLR_median.df<-merge(b1, NLR_median)

table(NLR_median.df$NLR_group) 

saveRDS(NLR_median.df, "NLR_median.df.contains.all.replicates.rds")
NLR_median.df.1<-filter(NLR_median.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(NLR_median.df.1, "NLR_median.df.contains.only.mean.rds")

#NLR by median
p21 = ggplot(NLR_median.df.1, aes(NLR_group_median,
                                  mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Neutrophil:Lymphocyte Ratio, split by median")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p21

pdf("Zane_and_Alex_9-21-23_update_Melanoma_NLR_split_by_median_vs_response.pdf")
p21
dev.off()

a<-NLR_median.df.1 %>% group_by(Site.renamed, BRR) %>% count(NLR_group_median) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(NLR_median.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ NLR_group_median, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R")
respond<-respond %>% filter(n>1) 
respond <-respond %>% filter(!Site.renamed  %in% c("peritoneum","visceral other"))
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ NLR_group_median, p.adjust.method="bonf")

b<-rbind(nr.1, respond.1)

write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_NLR_split_by_median_vs_response.csv")



####################--------------------------------
#eosinophil count 
Data<-readRDS("Data.rds")
Data$"% eos pre-RX"<-as.numeric(Data$"% eos pre-RX")
Data$"WBC count pre-Rx (x10^3)"<-as.numeric(Data$"WBC count pre-Rx (x10^3)")
Data<-Data %>% mutate(Eosinophil_count=Data$`WBC count pre-Rx (x10^3)`*(Data$`% eos pre-RX`/100))
Data<-Data %>%mutate(Eosinophil_count=Eosinophil_count*1000)
a<-Data %>% select("PID", "Eosinophil_count") %>%na.omit()
median(a$Eosinophil_count) #124.9



plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
table(Data$Eosinophil_count)
saveRDS(Data, "Data.rds")

a$Eosinophil_group[a$Eosinophil_count>median(a$Eosinophil_count)]<-"high"
a$Eosinophil_group[a$Eosinophil_count<=median(a$Eosinophil_count)]<-"low"

table(a$Eosinophil_group, useNA="always")
b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
Eosinophil.df<-merge(b1, a)

table(Eosinophil.df$Eosinophil_group) #110 high 120 normal

saveRDS(Eosinophil.df, "Eosinophil.df.contains.all.replicates.rds")
Eosinophil.df.1<-filter(Eosinophil.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(Eosinophil.df.1, "Eosinophil.df.contains.only.mean.rds")
Eosinophil.df.1<-readRDS("Eosinophil.df.contains.only.mean.rds")
Eosinophil.df.1<-Eosinophil.df.1 %>% filter(!Site.renamed %in% c( "visceral other", "peritoneum"))


#eosinophil
p21 = ggplot(Eosinophil.df.1, aes(Eosinophil_group,
                                  mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Eosinophil (split by median)")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p21

pdf("Zane_and_Alex_9-21-23_update_Melanoma_Eosinophil_count_vs_response_color.pdf")
p21
dev.off()

a<-Eosinophil.df.1 %>% group_by(Site.renamed, BRR) %>% count(Eosinophil_group) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Eosinophil.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ Eosinophil_group, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_9-21-23_update_wilcox_test_for_Eosinophil_count_vs_response.csv")


#Zane wants to re-run again but not separating by BRR 9-15-22-----------------------------
p21.1 = ggplot(Eosinophil.df.1, aes(Eosinophil_group,
                                    mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  xlab("Eosinophil (split by median)")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p21.1

pdf("Zane_and_Alex_9-21-23_update_Melanoma_Eosinophil_count_not_separated_by_BRR_vs_response.pdf")
p21.1
dev.off()

a<-Eosinophil.df.1 %>% group_by(Site.renamed) %>% count(Eosinophil_group) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Eosinophil.df.1, a.1, by=c("Site.renamed"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed)%>% pairwise_wilcox_test(mean ~ Eosinophil_group, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_Eosinophil_count_not_separated_by_BRR_vs_response.csv")






##sun exposure---------------------------------------
#per Zane's text message, we do not need to look at sun exposure

###################################################################
##progression on adjuvant therapy (pd1 or braf)

Data<-readRDS("Data.rds")


#
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

# 

Prog.on.adjuvant<- Data %>% select("PID","Prior adjuvant immunotherapy  treatment")

b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
Prog.on.adjuvant.df<-merge(b1, Prog.on.adjuvant)

#subtype.df<-subtype.df %>% filter (`Melanoma subtype` !="unknown primary")
table(Prog.on.adjuvant.df$`Prior adjuvant immunotherapy  treatment`)
#lots of variables with Yes

saveRDS(Prog.on.adjuvant.df, "Prog.on.adjuvant.df.contains.all.replicates.rds")
Prog.on.adjuvant.df.1<-filter(Prog.on.adjuvant.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(Prog.on.adjuvant.df.1, "Prog.on.adjuvant.df.contains.only.mean.rds")
table(Prog.on.adjuvant.df.1$`Prior adjuvant immunotherapy  treatment`, useNA="always")

Prog.on.adjuvant.df.1<-Prog.on.adjuvant.df.1 %>% mutate(Prog.on.therapy=ifelse(Prog.on.adjuvant.df.1$`Prior adjuvant immunotherapy  treatment`=="N", "no", "yes"))
saveRDS(Prog.on.adjuvant.df.1, "Prog.on.adjuvant.df.contains.only.mean.rds")


#prog on adjuvant
p23 = ggplot(Prog.on.adjuvant.df.1, aes(Prog.on.therapy,
                                        mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Progression on adjuvant therapy")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p23

pdf("Zane_and_Alex_9-21-23_update_Melanoma_prog_on_adjuv_therapy_vs_response.pdf")
p23
dev.off()


a<-Prog.on.adjuvant.df.1 %>% group_by(Site.renamed, BRR) %>% count(Prog.on.therapy) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Prog.on.adjuvant.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ Prog.on.therapy, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_prog_on_adjuvant_count_vs_response.csv")

########################################################################
#line of metastatic immunotherapy (need to upload from other excel sheet.)

Data<-readRDS("Data.rds")
#read in line of immunotherapy. this was fixed on the form from Alex on 9-21-23 but I didn't want to re-load so i'm just 
#adding it manually here

df<-Data %>% select("PID", "Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)")
df$`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`<-as.character(df$`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`)
str(df)

#df<-df[,-2]
saveRDS(Data, "Data.rds")

#Line.of.therapy<- Data %>% select("PID","Line.of.therapy.combined")

b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
Line.of.therapy.df<-merge(b1, df)

#subtype.df<-subtype.df %>% filter (`Melanoma subtype` !="unknown primary")
table(Line.of.therapy.df$`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`)
#     1 2_and_3 
#####449     117

saveRDS(Line.of.therapy.df, "Line.of.therapy.df.contains.all.replicates.rds")
Line.of.therapy.df.1<-filter(Line.of.therapy.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(Line.of.therapy.df.1, "Line.of.therapy.df.contains.only.mean.rds")
table(Line.of.therapy.df.1$Line.of.therapy.combiend, useNA="always")

saveRDS(Line.of.therapy.df.1, "Line.of.therapy.df.contains.only.mean.rds")
Line.of.therapy.df.1<-readRDS("Line.of.therapy.df.contains.only.mean.rds")
Line.of.therapy.df.1<-na.omit(Line.of.therapy.df.1)
#line of therapy
p24 = ggplot(Line.of.therapy.df.1, aes(`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`,
                                       mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Line of therapy 0=first line 1=combined 2 and 3 line")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p24

pdf("Zane_and_Alex_9-21-23_update_Melanoma_line_of_therapy_vs_response.pdf")
p24
dev.off()



a<-Line.of.therapy.df.1 %>% group_by(Site.renamed, BRR) %>% count(`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Line.of.therapy.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group
a.2<-a.2 %>% rename(line.tx="Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)")
b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ line.tx, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_line_on_therapy_vs_response.csv")


######edited for zane's presentaiton

plot.colors2 = c('CR' = '#0000CC','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' ='#CC0000' ,
                 'NA' = '#C0C0C0')

Line.of.therapy.df.1<-readRDS("Line.of.therapy.df.contains.only.mean.rds")
Line.of.therapy.df.1<-Line.of.therapy.df.1 %>% rename(`Best Response by RECIST`=Best.Response.by.RECIST)
#Line.of.therapy.df.1$Line.of.therapy.combined[Line.of.therapy.df.1$Line.of.therapy.combined=="2_and_3"]<-"≥2"

#Line.of.therapy.df.1$Line.of.therapy.combined<-factor(Line.of.therapy.df.1$Line.of.therapy.combined, levels=c("1", "≥2"))
Line.of.therapy.df.1<-Line.of.therapy.df.1 %>% filter(!Site.renamed %in% c("adrenal", "visceral other", "peritoneum"))
#line of therapy
Line.of.therapy.df.1<-na.omit(Line.of.therapy.df.1)
p24 = ggplot(Line.of.therapy.df.1, aes(`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`,
                                       mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = `Best Response by RECIST`), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Line of therapy")+
  ylab("Lesion Response (%)")

p24

pdf("Line_of_therapy_9-21-23_update_by_response")
p24
dev.off()

#split into response
Response<-Line.of.therapy.df.1 %>% filter(BRR=="R")


p25 = ggplot(Response, aes(`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`,
                           mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = `Best Response by RECIST`), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  xlab("Line of therapy")+
  ylab("Lesion Response (%)")

p25


pdf("Zane_and_Alex_9-21-23_update_Melanoma_RESPONDERS.line_of_therapy_vs_response.pdf")
p25
dev.off()
b<-table(Response$Site.renamed, Response$`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`)
b<-as.data.frame(b)
b$Var2<-as.character(b$Var2)
write.csv(b, "Responders.Line.of.response.numbers.csv")

#NR
N.R<-Line.of.therapy.df.1 %>% filter(BRR=="NR")


p26 = ggplot(N.R, aes(`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`,
                      mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = `Best Response by RECIST`), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  xlab("Line of therapy")+
  ylab("Lesion Response (%)")

p26


pdf("Zane_and_Alex_9-21-23_update_Melanoma_NR.line_of_therapy_vs_response.pdf")
p26
dev.off()

a<-table(N.R$Site.renamed, N.R$`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`)
a<-as.data.frame(a)
a$Var2<-as.character(a$Var2)
write.csv(a, "N.R.Line.of.response.numbers.csv")

#Plot again, Zane wants without stratified by response 9-15-22--------------------------------


#line of therapy
p24.1 = ggplot(Line.of.therapy.df.1, aes(`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`,
                                         mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = `Best Response by RECIST`), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  xlab("Line of therapy (combined 2 and 3)")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p24.1

pdf("Zane_and_Alex_9-21-23_update_Melanoma_line_of_therapy_not_separated_by_BRR_vs_response_color.pdf")
p24.1
dev.off()


a<-Line.of.therapy.df.1 %>% group_by(Site.renamed) %>% count(`Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)`) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Line.of.therapy.df.1, a.1, by=c("Site.renamed"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group
a.2<-a.2 %>% rename(line.tx="Line of immunotherapy (0 = first line, 1 = 2nd or 3rd line)")

b<-a.2 %>% group_by(Site.renamed)%>% pairwise_wilcox_test(mean ~ line.tx, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_9-21-23_update_Pairwise_wilcox_test_for_line_on_therapy_not_separated_by_BRR_vs_response.csv")




























































































































































































######################################################################################
##repeat everyting for monotherapy
Data<-readRDS("Data.rds")
colnames(Data)
Mono<-Data %>% filter(`Ipi/nivo = 1, PD-1 = 0`=="0")




percent<-c( "adrenal.pct", "adrenal 2.pct", "visceral other.pct", "visceral other 2.pct", "peritoneum.pct", "soft tissue.pct","soft tissue 2.pct","nodes.pct","nodes 2.pct","brain.pct","brain 2.pct","lung.pct", "lung 2.pct","liver.pct","liver 2.pct" )     
Mono1 = Mono[,c("PID","Best Response by RECIST","BRR",
                percent)]


Mono1 = reshape2::melt(Mono1)
colnames(Mono1)[(ncol(Mono1)-1):ncol(Mono1)] = c('Site','Change.pct')
saveRDS(Mono1,"Mono1")


## lesion 1 or 2
Mono1$Lesion.number = "LS1"
Mono1$Lesion.number[grep('2[.]pct$', Mono1$Site)] = 'LS2'
table(Mono1$Lesion.number)
#LS1  LS2 
#1088  952 

saveRDS(Mono1,"Mono1")
## use line to connect lesions from the same pt instead of calling it site and site 2
#Mono$Site = gsub('[.]2[.]pct$','.pct',Mono$Site)
#print(table(Mono$Site))

Mono1$Group = paste0(Mono1$BRR,'.',Mono1$Lesion.number)


saveRDS(Mono1,"Mono1.rds")

##-------------------------------------------------
## for each site, take mean of multiple lesions 

#for some reason adrenal 2 has a weird name so change it to proper
#Mono1<- lapply(Mono1, function(x){gsub("adrenal 2.adrenal.pct", "adrenal 2.pct", x)})
Mono1<-as.data.frame(Mono1)
Mono1$Lesion.number[grep('2[.]pct$', Mono1$Site)] = 'LS2'
#need to fix Mono plot 2 as it isn't taking average across multiple lesions
Mono1<-Mono1 %>% group_by(PID)
saveRDS(Mono1,"Mono1.rds")
#i think the prolem is that theree are two sites. Change all sites to the same and rely on lesion for future separation
Mono2<-data.frame(Mono1)
Mono2[Mono2=="adrenal 2.pct"]<-"adrenal.pct"
Mono2[Mono2=="visceral other 2.pct"]<-"visceral other.pct" 
Mono2[Mono2=="soft tissue 2.pct"]<-"soft tissue.pct" 
Mono2[Mono2=="nodes 2.pct"]<-"nodes.pct" 
Mono2[Mono2=="brain 2.pct"]<-"brain.pct" 
Mono2[Mono2=="lung 2.pct"]<-"lung.pct" 
Mono2[Mono2=="liver 2.pct"]<-"liver.pct" 
saveRDS(Mono2, "Mono.2.rds")

Mono2<-Mono2 %>% group_by(Site)
Mono2<-Mono2[!is.na(Mono2$Change.pct),] 
str(Mono2)
Mono2$Change.pct<-as.numeric(Mono2$Change.pct)
#note two ways to perform mean expression and I verified each works (a and b below)
b<-Mono2 %>% group_by(PID,Site) %>%
  summarise(mean = mean(Change.pct), na.rm=TRUE)
b

a<-aggregate(Mono2$Change.pct, by=list(PID=Mono2$PID, Site=Mono2$Site),Mono=Mono2,FUN=mean, na.rm=TRUE)


b1<-merge(Mono2, b)
saveRDS(b1, "average_percent_change_by_patient_and_Site_monotherapy.rds")
b1<-readRDS("average_percent_change_by_patient_and_Site_monotherapy.rds")
#b1 or average_percent_change_by_patient_and_Site is average while Mono2 is individual
dim(b1) #[1][1] 236   9
my.levels<-c("PD","SD","PR","CR")
#b1<-b1 %>% rename(Best.Response.by.RECIST=Best.response.by.RECIST)
b1$Best.Response.by.RECIST<-factor(b1$Best.Response.by.RECIST, my.levels)
b1$Best.Response.by.RECIST<-na.omit(b1$Best.Response.by.RECIST)

#editing this because Zane wants the facet wrap not to say ".pct"
b1<-b1 %>% mutate(Site.renamed=b1$Site)
b1$Site.renamed<-str_remove(b1$Site.renamed, ".pct")
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

plot.colors2 = c('CR' = '#0000CC','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' ='#CC0000' ,
                 'NA' = '#C0C0C0')

rm(p1, p2)
rm(p1, p2)
my_comparisons <- list( c("NR.LS1", "NR.LS2"), c("R.LS1", "R.LS2") )
p1 = ggplot(b1, aes(Group,
                    Change.pct)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_point(aes(color = Best.Response.by.RECIST,
                 shape = Lesion.number)) +
  geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  stat_compare_means(comparisons = my_comparisons, method="wilcox")+
  ggtitle("Monotherapy_Plots")

p1

saveRDS(b1, "average_percent_change_by_patient_and_Site_monotherapy.rds")
#tiff("Zane_and_Alex_boxplot_Change_percent_between_lesion_and_Group.tiff", res=300, width=8, height=6, units='in')
#p1
#dev.off()

#edited for Zane's talk 9-26-22, and 9-21-23
b2<-b1 %>% rename(`Best Response by RECIST`=Best.Response.by.RECIST)
#edited for Zane's presentation
#•	Left and right figure, remove adrenal, visceral, peritoneum
b2<-b2 %>% rename(`Lesion Number`=Lesion.number)
b2<-b2 %>% filter(!Site.renamed %in% c( "visceral other", "peritoneum"))
b2$`Best Response by RECIST`<-factor(b2$`Best Response by RECIST`, levels=c("CR","PR","SD","PD"))
p1 = ggplot(b2, aes(Group,
                    Change.pct)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_point(aes(color = `Best Response by RECIST`,
                 shape = `Lesion Number`)) +
  geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  stat_compare_means(comparisons = my_comparisons, method="wilcox")+
  ylab("Lesion Response (%)")

p1


pdf("Zane_and_Alex_9-21-23_update_boxplot_Monotherapy_Change_percent_between_lesion_and_Group.pdf")
p1
dev.off()


p2 = ggplot(b2, aes(BRR,
                    mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = `Best Response by RECIST`), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  stat_compare_means(comparisons = list(c('NR','R'))) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  xlab("Overall Best Response by RECIST")+
  ylab("Lesion Response (%)")

p2


pdf("Zane_and_Alex_edited_boxplot_Monotherapy_update_9-21-23_Change_percent_between_only_group.pdf")
p2
dev.off()


a<-b2 %>% group_by(BRR, Site.renamed) %>% summarize (median(mean), total_count=n(), sd=sd(mean))
a<-as.data.frame(a)

write.xlsx(a, "Monotherapy_Median_and_count.by.RECIST_location_with_sd.xlsx")
saveRDS(b2, "average_percent_change_by_patient_and_Site_reduced_sites_monotherapy.rds")


b1<-readRDS("average_percent_change_by_patient_and_Site_monotherapy.rds")
#Direction plots
b1$change.direction = NA
b1$change.direction[which(b1$Change.pct>0)] = 'up'
b1$change.direction[which(b1$Change.pct<0)] = 'down'
b1$change.direction[which(b1$Change.pct==0)] = 'nochange'
head(b1)

table(b1$change.direction)
dim(b1) #236, 10
sum(is.na(b1$Change.pct)) #0
length(unique(b1$Site)) #8 sites



data.plot.2.2 = NULL 
for(my.response in c('PD','SD','PR','CR')) {
  print(my.response)
  
  my.df = NULL 
  my.df = data.frame(table(
    b1[b1$Best.Response.by.RECIST==my.response,
       c('PID',
         'change.direction'),drop=F]))
  my.df = reshape2::dcast(PID ~ change.direction,
                          data = my.df,
                          value.var = 'Freq',
                          fun.aggregate = sum)
  my.df$Best.Response.by.RECIST = my.response
  rm(x)
  x = my.df[,c('PID','Best.Response.by.RECIST')]
  for( my.col in c('down','nochange','up')) {
    print(my.col)
    x$new = NA
    if(length(which(colnames(my.df)==my.col))>0) {
      x$new = my.df[,my.col]
      
    }
    colnames(x)[which(colnames(x)=='new')] = my.col
  }
  my.df = x
  data.plot.2.2 = rbind(data.plot.2.2,
                        my.df)
}
data.plot.2.2[is.na(data.plot.2.2)] = 0

## bar graphs
data.plot.2.2.melt = NULL 
data.plot.2.2.melt = reshape2::melt(data.plot.2.2)
colnames(data.plot.2.2.melt) = c('PID',
                                 'Best.Response.by.RECIST',
                                 'change.direction',
                                 'sites')
data.plot.2.2.melt$Best.Response.by.RECIST = 
  factor(data.plot.2.2.melt$Best.Response.by.RECIST,
         levels = c('PD','SD','PR','CR'))
p4.2 = NULL 
p4.2 = ggplot(data.plot.2.2.melt,aes(PID,sites)) +
  geom_bar(aes(fill = change.direction), width = 0.6,
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values = c('up' = '#0000CC',
                               'nochange' = '#C0C0C0',
                               'down' = '#CC0000')) + 
  theme_pubr(x.text.angle = 90) +
  facet_grid(~ Best.Response.by.RECIST, scale = 'free_x', space = 'free_x')
p4.2

pdf("Monotherapy_update_9-21-23_update_change_direction_plot.pdf")
p4.2
dev.off()

## tables 
data.plot.2.3 = data.plot.2.2
data.plot.2.3$total = apply(data.plot.2.3[,-c(1:2)], 1, sum)
data.plot.2.3[,3:5] = data.plot.2.3[,3:5] / data.plot.2.3[,6]
data.plot.2.3$category = NA
data.plot.2.3$category[which(data.plot.2.3$down==0 &
                               data.plot.2.3$nochange==0 &
                               data.plot.2.3$up == 1)] = 'uniform progression' #change 9-26 to uniform progression
data.plot.2.3$category[which(data.plot.2.3$down==1 &
                               data.plot.2.3$nochange==0 &
                               data.plot.2.3$up == 0)] = 'uniform regression'#change 9-26 to uniform regression
data.plot.2.3$category[which(data.plot.2.3$down==0 &
                               data.plot.2.3$nochange==1 &
                               data.plot.2.3$up == 0)] = 'all.nochange'
data.plot.2.3$category[is.na(data.plot.2.3$category)] = 'mixed response'

table(data.plot.2.3$category)
#  mixed response uniform progression  uniform regression 
###20                  30                  32 


data.plot.2.3 = data.frame(
  table(data.plot.2.3[,c('Best.Response.by.RECIST','category')]),
  stringsAsFactors = F)
data.plot.2.3$Best.Response.by.RECIST = factor(
  data.plot.2.3$Best.Response.by.RECIST,
  levels = c('PD','SD','PR','CR')
)
p4.3 = NULL 
p4.3 = ggplot(data.plot.2.3, aes(Best.Response.by.RECIST,
                                 Freq)) +
  geom_bar(aes(fill = category),width = 0.6,
           position = 'stack', stat = 'identity',
           color = '#000000') +
  theme_pubr() +
  scale_fill_manual(values = c('uniform regression' = '#0000CC',
                               #'all.nochange' = '#C0C0C0',
                               'mixed response' = '#00CC00',
                               'uniform progression' = '#CC0000')) 
p4.3



pdf("Zane_and_Alex_change_direction_by_response_MONOTHERAPY_9-21-23.pdf")
p4.3
dev.off()


write.csv(data.plot.2.3, "Change.direction.by.response.n.group_9-21-23-MONOTHERAPY.csv")













































###############################################################################
###repeat everything for combo
colnames(Data)
Combo.only<-Data %>% filter(`Mono (0) vs Combo(1)`=="1")



percent<-c( "adrenal.pct", "adrenal 2.pct", "visceral other.pct", "visceral other 2.pct", "peritoneum.pct", "soft tissue.pct","soft tissue 2.pct","nodes.pct","nodes 2.pct","brain.pct","brain 2.pct","lung.pct", "lung 2.pct","liver.pct","liver 2.pct" )     
Combo.only1 = Combo.only[,c("PID","Best Response by RECIST","BRR",
                            percent)]


Combo.only1 = reshape2::melt(Combo.only1)
colnames(Combo.only1)[(ncol(Combo.only1)-1):ncol(Combo.only1)] = c('Site','Change.pct')
saveRDS(Combo.only1,"Combo.only1")


## lesion 1 or 2
Combo.only1$Lesion.number = "LS1"
Combo.only1$Lesion.number[grep('2[.]pct$', Combo.only1$Site)] = 'LS2'
table(Combo.only1$Lesion.number)
#LS1 LS2 
#656 574

saveRDS(Combo.only1,"Combo.only1")
## use line to connect lesions from the same pt instead of calling it site and site 2
#Mono$Site = gsub('[.]2[.]pct$','.pct',Mono$Site)
#print(table(Mono$Site))

Combo.only1$Group = paste0(Combo.only1$BRR,'.',Combo.only1$Lesion.number)


saveRDS(Combo.only1,"Combo.only1.rds")

##-------------------------------------------------
## for each site, take mean of multiple lesions 

#for some reason adrenal 2 has a weird name so change it to proper
#Combo.only1<- lapply(Combo.only1, function(x){gsub("adrenal 2.adrenal.pct", "adrenal 2.pct", x)})
Combo.only1<-as.Mono.frame(Combo.only1)
Combo.only1$Lesion.number[grep('2[.]pct$', Combo.only1$Site)] = 'LS2'
#need to fix Mono plot 2 as it isn't taking average across multiple lesions
Combo.only1<-Combo.only1 %>% group_by(PID)
saveRDS(Combo.only1,"Combo.only1.rds")
#i think the prolem is that theree are two sites. Change all sites to the same and rely on lesion for future separation
Combo.only2<-data.frame(Combo.only1)
Combo.only2[Combo.only2=="adrenal 2.pct"]<-"adrenal.pct"
Combo.only2[Combo.only2=="visceral other 2.pct"]<-"visceral other.pct" 
Combo.only2[Combo.only2=="soft tissue 2.pct"]<-"soft tissue.pct" 
Combo.only2[Combo.only2=="nodes 2.pct"]<-"nodes.pct" 
Combo.only2[Combo.only2=="brain 2.pct"]<-"brain.pct" 
Combo.only2[Combo.only2=="lung 2.pct"]<-"lung.pct" 
Combo.only2[Combo.only2=="liver 2.pct"]<-"liver.pct" 
saveRDS(Combo.only2, "Combo.only2.rds")

Combo.only2<-Combo.only2 %>% group_by(Site)
Combo.only2<-Combo.only2[!is.na(Combo.only2$Change.pct),] 
str(Combo.only2)
Combo.only2$Change.pct<-as.numeric(Combo.only2$Change.pct)
#note two ways to perform mean expression and I verified each works (a and b below)
b<-Combo.only2 %>% group_by(PID,Site) %>%
  summarise(mean = mean(Change.pct), na.rm=TRUE)
b

a<-aggregate(Combo.only2$Change.pct, by=list(PID=Combo.only2$PID, Site=Combo.only2$Site),Mono=Combo.only2,FUN=mean, na.rm=TRUE)


b1<-merge(Combo.only2, b)
saveRDS(b1, "average_percent_change_by_patient_and_Site_combo_only.rds")
b1<-readRDS("average_percent_change_by_patient_and_Site_combo_only.rds")
#b1 or average_percent_change_by_patient_and_Site is average while Combo.only2 is individual
dim(b1) #[1][1] 236   9
my.levels<-c("PD","SD","PR","CR")
b1<-b1 %>% rename(Best.Response.by.RECIST=Best.response.by.RECIST)
b1$Best.Response.by.RECIST<-factor(b1$Best.Response.by.RECIST, my.levels)
b1$Best.Response.by.RECIST<-na.omit(b1$Best.Response.by.RECIST)

#editing this because Zane wants the facet wrap not to say ".pct"
b1<-b1 %>% mutate(Site.renamed=b1$Site)
b1$Site.renamed<-str_remove(b1$Site.renamed, ".pct")
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

plot.colors2 = c('CR' = '#0000CC','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' ='#CC0000' ,
                 'NA' = '#C0C0C0')

rm(p1, p2)
rm(p1, p2)
my_comparisons <- list( c("NR.LS1", "NR.LS2"), c("R.LS1", "R.LS2") )
p1 = ggplot(b1, aes(Group,
                    Change.pct)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_point(aes(color = Best.Response.by.RECIST,
                 shape = Lesion.number)) +
  geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  stat_compare_means(comparisons = my_comparisons, method="wilcox")

p1

saveRDS(b1, "average_percent_change_by_patient_and_Site_combo_only.rds")
#tiff("Zane_and_Alex_boxplot_Change_percent_between_lesion_and_Group.tiff", res=300, width=8, height=6, units='in')
#p1
#dev.off()

#edited for Zane's talk 9-26-22
b2<-b1 %>% rename(`Best Response by RECIST`=Best.Response.by.RECIST)
#edited for Zane's presentation
#•	Left and right figure, remove adrenal, visceral, peritoneum
b2<-b2 %>% rename(`Lesion Number`=Lesion.number)
b2<-b2 %>% filter(!Site.renamed %in% c( "visceral other", "peritoneum"))
b2$`Best Response by RECIST`<-factor(b2$`Best Response by RECIST`, levels=c("CR","PR","SD","PD"))
p1 = ggplot(b2, aes(Group,
                    Change.pct)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_point(aes(color = `Best Response by RECIST`,
                 shape = `Lesion Number`)) +
  geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  stat_compare_means(comparisons = my_comparisons, method="wilcox")+
  ylab("Lesion Response (%)")

p1


pdf("Zane_and_Alex_Edited_02-07-23_boxplot_combo_only_Change_percent_between_lesion_and_Group.pdf")
p1
dev.off()


p2 = ggplot(b2, aes(BRR,
                    mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = `Best Response by RECIST`), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  stat_compare_means(comparisons = list(c('NR','R'))) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  xlab("Overall Best Response by RECIST")+
  ylab("Lesion Response (%)")

p2


pdf("Zane_and_Alex_edited_boxplot_combo_only_update_02-07-23_Change_percent_between_only_group.pdf")
p2
dev.off()


a<-b2 %>% group_by(BRR, Site.renamed) %>% summarize (median(mean), total_count=n(), sd=sd(mean))
a<-as.data.frame(a)

write.xlsx(a, "Combo_only_Median_and_count.by.RECIST_location_with_sd.xlsx")
saveRDS(b2, "average_percent_change_by_patient_and_Site_reduced_sites_combo_only.rds")

