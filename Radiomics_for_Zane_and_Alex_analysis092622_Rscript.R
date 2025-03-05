
#Radiomics analysis for Zane and Alex
#Becky OLD analysis


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
library(SHAPforxgboost)
library(e1071)
library(glmnet)
library(stringr)
library("rstatix")
#devtools::install_github("liuyanguu/SHAPforxgboost")




Clinical<-readxl::read_excel('Melanoma_IpiNivo_THIS_IS_THE_MOST_UPTODATE.xlsx', 
                             sheet = "Master")
str(Clinical) 
Clinical<-as.data.frame(Clinical)


to.remove<-c("...104","...105","...106","...107", "...108","...109")
#it added on a few extra columns since there were 2 extra columns with random stuff in it. Remove that
Clinical<-Clinical[, !(names(Clinical) %in% to.remove)]

#need to edit some of the columns to be the right data type
Clinical$`Date of diagnosis (mm/dd/yyyy)`<-as.numeric(Clinical$`Date of diagnosis (mm/dd/yyyy)`)
Clinical$`Date of diagnosis (mm/dd/yyyy)`<-XLDateToPOSIXct(Clinical$`Date of diagnosis (mm/dd/yyyy)`)

Clinical[, 17:42] <- lapply(17:42, function(x) as.numeric(Clinical[[x]]))
str(Clinical)
Clinical$`total # of Nivo mono doses (maintenance)`<-as.numeric(Clinical$`total # of Nivo mono doses (maintenance)`)
Clinical[, 53:86] <- lapply(53:86, function(x) as.numeric(Clinical[[x]]))

Clinical$`Date of first scan`<-as.numeric(Clinical$`Date of first scan`)
Clinical$`Date of first scan`<-XLDateToPOSIXct(Clinical$`Date of first scan`)

Clinical$`Date of scan at best response`<-as.numeric(Clinical$`Date of scan at best response`)
Clinical$`Date of scan at best response`<-XLDateToPOSIXct(Clinical$`Date of scan at best response`)

Clinical$`Date of first scan showing disease progression (mm/dd/yyyy)`<-as.numeric(Clinical$`Date of first scan showing disease progression (mm/dd/yyyy)`)
Clinical$`Date of first scan showing disease progression (mm/dd/yyyy)`<-XLDateToPOSIXct(Clinical$`Date of first scan showing disease progression (mm/dd/yyyy)`)
Clinical$`Days btwn Tx t0 and Scan at progression`<-as.numeric(Clinical$`Days btwn Tx t0 and Scan at progression`)

Clinical$`Date of Scan Before Next Tx/Most Recent`<-as.numeric(Clinical$`Date of Scan Before Next Tx/Most Recent`)
Clinical$`Date of Scan Before Next Tx/Most Recent`<-XLDateToPOSIXct(Clinical$`Date of Scan Before Next Tx/Most Recent`)

Clinical$`Target lesion % change from baseline`<-as.numeric(Clinical$`Target lesion % change from baseline`)

Clinical$`Follow-up 1 lesion measurement/Baseline lesion measurement - 1`<-as.numeric(Clinical$`Follow-up 1 lesion measurement/Baseline lesion measurement - 1`)
Clinical$`Follow-up 2 lesion measurement/Baseline lesion measurement - 1`<-as.numeric(Clinical$`Follow-up 2 lesion measurement/Baseline lesion measurement - 1`)
Clinical$`Follow-up 3 lesion measurement/Baseline lesion measurement - 1`<-as.numeric(Clinical$`Follow-up 3 lesion measurement/Baseline lesion measurement - 1`)

saveRDS(Clinical, "Clinical.rds")
write.csv(Clinical,"Clinical.csv")

#Filter samples to remove NE and N/A RECIST--------------------------------------
table(Clinical$`Best Response by RECIST`)
# CR N/A  NE  PD  PR  SD 
##6   3   1  54  18  13 

sum(is.na(Clinical$`Best Response by RECIST`))
#none that are NAs


#remove samples with NE
rm(x)
x = which(!is.na(Clinical$`Best Response by RECIST`) &
            Clinical$`Best Response by RECIST`=='NE')
print(paste0('index of NE row is:', x))

Clinical$`Best Response by RECIST`[x] = NA

table(Clinical$`Best Response by RECIST`)
# CR N/A  PD  PR  SD 
###6   3  54  18  13

#need to remove the 5 samples that have N/A as their RECIST score

Clinical<-Clinical %>% filter(Clinical$`Best Response by RECIST`!="N/A")

table(Clinical$`Best Response by RECIST`)
#CR PD PR SD 
#6 54 18 13 
##Collapse Responses into NR or R----------------------------------------
#according to Zane, keep SD as responder
Clinical$BRR = NA
Clinical$BRR[which(Clinical$`Best Response by RECIST` %in% c('PD'))] = 'NR'
Clinical$BRR[which(Clinical$`Best Response by RECIST` %in% c('PR','CR','SD'))] = 'R'

table(Clinical$BRR)
#NR  R 
#54 37 

table(Clinical$BRR)[2] / nrow(Clinical[!is.na(Clinical$BRR),])
##0.4065934 

sum(is.na(Clinical$BRR))
## 0 (removed?)

saveRDS(Clinical, "Clinical.rds")

#add a new column with P001 as 1st entry, etc.
Clinical$PID = paste0('P',stringr::str_pad(1:nrow(Clinical), 3, pad = "0"))

data.plot = NULL
data.plot = Clinical 

my.levels<-c("PD","SD","PR","CR")
data.plot$`Best Response by RECIST`=factor(data.plot$`Best Response by RECIST`, levels=my.levels)
data.plot<-data.plot %>% filter(!is.na(`Target lesion % change from baseline`))
table(data.plot$`Target lesion % change from baseline`, data.plot$BRR)

##----------------------------------------------------------------
#find shrinkage of tumors at indivudal target lesion. Already the Target lesion % change from baseline has been calculated

plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
#'R' = '#CC0000', 'NR' = '#0000CC'

p1 = ggplot(data.plot, aes(x=reorder(PID, -`Target lesion % change from baseline`), y=`Target lesion % change from baseline`)) +
  geom_bar(aes(fill = `Best Response by RECIST`), color = '#000000', width = 0.6, stat = 'identity') +
  scale_fill_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90)+
  xlab("Patients")+
  geom_hline(yintercept = -30, linetype = 'dashed')
#scale_fill_discrete(drop=TRUE)

p1 

saveRDS(data.plot,"data.plot.rds")

tiff("Zane_and_Alex_Waterfall_Change_in_Target_Lesion_RECIST.tiff", res=300, width=12, height=6, units='in')
p1
dev.off()

my.levels<-c("PD","SD","PR","CR")
data.plot$`Best Response by RECIST`=factor(data.plot$`Best Response by RECIST`, levels=my.levels)

#Do the same but this time by NR vs R

p2 = ggplot(data.plot, aes(x=BRR, y=`Target lesion % change from baseline`)) +
  geom_boxplot(outlier.shape = NA) +
  theme_pubr(x.text.angle = 90)+
  xlab("Patients")+
  geom_hline(yintercept = -30, linetype = 'dashed')+
  geom_point(aes(fill = `Best Response by RECIST`), size = 5,  color = '#000000', shape = 21, position = "jitter")+
  scale_fill_manual(values = plot.colors) +
  stat_compare_means(comparisons = list(c('NR','R')),
                     method = 'wilcox.test')



p2

tiff("Zane_and_Alex_Boxplot_Change_in_Target_Lesion_RECIST_NR_v_R.tiff", res=300, width=12, height=6, units='in')
p2
dev.off()


p3 = ggplot(data.plot, aes(x=`Best Response by RECIST`, y=`Target lesion % change from baseline`)) +
  geom_boxplot(outlier.shape = NA) +
  theme_pubr(x.text.angle = 90)+
  xlab("Patients")+
  geom_hline(yintercept = -30, linetype = 'dashed')+
  geom_point(aes(fill = BRR), size = 5, shape = 21, position = position_jitterdodge())


p3


tiff("Zane_and_Alex_Boxplot_Change_in_Target_Lesion_RECIST_all_RECIST.tiff", res=300, width=8, height=6, units='in')
p3
dev.off()


saveRDS(data.plot,"data.plot.rds")

my.site<-c("adrenal", "adrenal 2.adrenal.pct", "visceral other","visceral other 2", "peritoneum", "soft tissue",
           "soft tissue 2", "nodes", "nodes 2", "brain", "brain 2", "lung", "lung 2", "liver", "liver 2")


#for some reason, Sunny's loop wasn't working so I'm doing it by hand
adrenal.df<-data.plot[,c("adrenal (baseline)", "adrenal (best response)")]
adrenal.df<-as.data.frame(adrenal.df)
adrenal.df$subtract= adrenal.df[,2] - adrenal.df[,1]
adrenal.df$ratio = adrenal.df[,2] / adrenal.df[,1]
adrenal.df$pct = (adrenal.df[,2] / adrenal.df[,1] - 1) * 100
colnames(adrenal.df)[-c(1:2)] = paste0("adrenal",'.',colnames(adrenal.df)[-c(1:2)])
adrenal.df<-as.data.frame(adrenal.df)

data.plot = cbind(data.plot, adrenal.df[,-c(1:2)])
#adrenal 2
adrenal.2.df<-data.plot[,c("adrenal 2 (baseline)", "adrenal 2 (best response)")]
adrenal.2.df$subtract = adrenal.2.df[,2] - adrenal.2.df[,1]
adrenal.2.df$ratio = adrenal.2.df[,2] / adrenal.2.df[,1]
adrenal.2.df$pct = (adrenal.2.df[,2] / adrenal.2.df[,1] - 1) * 100
colnames(adrenal.2.df)[-c(1:2)] = paste0("adrenal 2",'.',colnames(adrenal.2.df)[-c(1:2)])
adrenal.2.df<-as.data.frame(adrenal.2.df)
data.plot = cbind(data.plot, adrenal.2.df[,-c(1:2)])

#visceral other
visceral.other.df<-data.plot[,c("visceral other (baseline)" , "visceral other (best response)")]
visceral.other.df$subtract = visceral.other.df[,2] - visceral.other.df[,1]
visceral.other.df$ratio = visceral.other.df[,2] / visceral.other.df[,1]
visceral.other.df$pct = (visceral.other.df[,2] / visceral.other.df[,1] - 1) * 100
colnames(visceral.other.df)[-c(1:2)] = paste0("visceral other",'.',colnames(visceral.other.df)[-c(1:2)])

data.plot = cbind(data.plot, visceral.other.df[,-c(1:2)])


#visceral other
visceral.other.2.df<-data.plot[,c("visceral other 2 (baseline)" , "visceral other 2 (best response)")]
visceral.other.2.df$subtract = visceral.other.2.df[,2] - visceral.other.2.df[,1]
visceral.other.2.df$ratio = visceral.other.2.df[,2] / visceral.other.2.df[,1]
visceral.other.2.df$pct = (visceral.other.2.df[,2] / visceral.other.2.df[,1] - 1) * 100
colnames(visceral.other.2.df)[-c(1:2)] = paste0("visceral other 2",'.',colnames(visceral.other.2.df)[-c(1:2)])

data.plot = cbind(data.plot, visceral.other.2.df[,-c(1:2)])

#peritoneum
peritoneum.df<-data.plot[,c("peritoneum (baseline)" , "peritoneum (best response)")]
peritoneum.df$subtract = peritoneum.df[,2] - peritoneum.df[,1]
peritoneum.df$ratio = peritoneum.df[,2] / peritoneum.df[,1]
peritoneum.df$pct = (peritoneum.df[,2] / peritoneum.df[,1] - 1) * 100
colnames(peritoneum.df)[-c(1:2)] = paste0("peritoneum",'.',colnames(peritoneum.df)[-c(1:2)])

data.plot = cbind(data.plot, peritoneum.df[,-c(1:2)])

#softTissue1

soft.tissue.df<-data.plot[,c("soft tissue (baseline)" , "soft tissue (best response)")]
soft.tissue.df$subtract = soft.tissue.df[,2] - soft.tissue.df[,1]
soft.tissue.df$ratio = soft.tissue.df[,2] / soft.tissue.df[,1]
soft.tissue.df$pct = (soft.tissue.df[,2] / soft.tissue.df[,1] - 1) * 100
colnames(soft.tissue.df)[-c(1:2)] = paste0("soft tissue",'.',colnames(soft.tissue.df)[-c(1:2)])

data.plot = cbind(data.plot, soft.tissue.df[,-c(1:2)])


#soft tissue 2
soft.tissue.2.df<-data.plot[,c("soft tissue 2 (baseline)" , "soft tissue 2 (best response)")]
soft.tissue.2.df$subtract = soft.tissue.2.df[,2] - soft.tissue.2.df[,1]
soft.tissue.2.df$ratio = soft.tissue.2.df[,2] / soft.tissue.2.df[,1]
soft.tissue.2.df$pct = (soft.tissue.2.df[,2] / soft.tissue.2.df[,1] - 1) * 100
colnames(soft.tissue.2.df)[-c(1:2)] = paste0("soft tissue 2",'.',colnames(soft.tissue.2.df)[-c(1:2)])

data.plot = cbind(data.plot, soft.tissue.2.df[,-c(1:2)])

#nodes
nodes.df<-data.plot[,c("nodes (baseline)" , "nodes (best response)")]
nodes.df$subtract = nodes.df[,2] - nodes.df[,1]
nodes.df$ratio = nodes.df[,2] / nodes.df[,1]
nodes.df$pct = (nodes.df[,2] / nodes.df[,1] - 1) * 100
colnames(nodes.df)[-c(1:2)] = paste0("nodes",'.',colnames(nodes.df)[-c(1:2)])

data.plot = cbind(data.plot, nodes.df[,-c(1:2)])


#nodes 2
nodes.2.df<-data.plot[,c("nodes 2 (baseline)" , "nodes 2 (best response)")]
nodes.2.df$subtract = nodes.2.df[,2] - nodes.2.df[,1]
nodes.2.df$ratio = nodes.2.df[,2] / nodes.2.df[,1]
nodes.2.df$pct = (nodes.2.df[,2] / nodes.2.df[,1] - 1) * 100
colnames(nodes.2.df)[-c(1:2)] = paste0("nodes 2",'.',colnames(nodes.2.df)[-c(1:2)])

data.plot = cbind(data.plot, nodes.2.df[,-c(1:2)])


#brain 

brain.df<-data.plot[,c("brain (baseline)" , "brain (best response)")]
brain.df$subtract = brain.df[,2] - brain.df[,1]
brain.df$ratio = brain.df[,2] / brain.df[,1]
brain.df$pct = (brain.df[,2] / brain.df[,1] - 1) * 100
colnames(brain.df)[-c(1:2)] = paste0("brain",'.',colnames(brain.df)[-c(1:2)])

data.plot = cbind(data.plot, brain.df[,-c(1:2)])


#brain 2
brain.2.df<-data.plot[,c("brain 2 (baseline)" , "brain 2 (best response)")]
brain.2.df$subtract = brain.2.df[,2] - brain.2.df[,1]
brain.2.df$ratio = brain.2.df[,2] / brain.2.df[,1]
brain.2.df$pct = (brain.2.df[,2] / brain.2.df[,1] - 1) * 100
colnames(brain.2.df)[-c(1:2)] = paste0("brain 2",'.',colnames(brain.2.df)[-c(1:2)])

data.plot = cbind(data.plot, brain.2.df[,-c(1:2)])

#lung

lung.df<-data.plot[,c("lung (baseline)" , "lung (best response)")]
lung.df$subtract = lung.df[,2] - lung.df[,1]
lung.df$ratio = lung.df[,2] / lung.df[,1]
lung.df$pct = (lung.df[,2] / lung.df[,1] - 1) * 100
colnames(lung.df)[-c(1:2)] = paste0("lung",'.',colnames(lung.df)[-c(1:2)])

data.plot = cbind(data.plot, lung.df[,-c(1:2)])



#lung 2

lung.2.df<-data.plot[,c("lung 2 (baseline)" , "lung 2 (best response)")]
lung.2.df$subtract = lung.2.df[,2] - lung.2.df[,1]
lung.2.df$ratio = lung.2.df[,2] / lung.2.df[,1]
lung.2.df$pct = (lung.2.df[,2] / lung.2.df[,1] - 1) * 100
colnames(lung.2.df)[-c(1:2)] = paste0("lung 2",'.',colnames(lung.2.df)[-c(1:2)])

data.plot = cbind(data.plot, lung.2.df[,-c(1:2)])

#liver
liver.df<-data.plot[,c("liver (baseline)" , "liver (best response)")]
liver.df$subtract = liver.df[,2] - liver.df[,1]
liver.df$ratio = liver.df[,2] / liver.df[,1]
liver.df$pct = (liver.df[,2] / liver.df[,1] - 1) * 100
colnames(liver.df)[-c(1:2)] = paste0("liver",'.',colnames(liver.df)[-c(1:2)])

data.plot = cbind(data.plot, liver.df[,-c(1:2)])

#liver 2
liver.2.df<-data.plot[,c("liver 2 (baseline)" , "liver 2 (best response)")]
liver.2.df$subtract = liver.2.df[,2] - liver.2.df[,1]
liver.2.df$ratio = liver.2.df[,2] / liver.2.df[,1]
liver.2.df$pct = (liver.2.df[,2] / liver.2.df[,1] - 1) * 100
colnames(liver.2.df)[-c(1:2)] = paste0("liver 2",'.',colnames(liver.2.df)[-c(1:2)])

data.plot = cbind(data.plot, liver.2.df[,-c(1:2)])




saveRDS(data.plot,"data.plot.rds") #data.plot has all the change in sizes among tissues

percent<-c( "adrenal.pct", "adrenal 2.pct", "visceral other.pct", "visceral other 2.pct", "peritoneum.pct", "soft tissue.pct","soft tissue 2.pct","nodes.pct","nodes 2.pct","brain.pct","brain 2.pct","lung.pct", "lung 2.pct","liver.pct","liver 2.pct" )     
data.plot1 = data.plot[,c("PID","Best Response by RECIST","BRR",
                          percent)]


data.plot1 = reshape2::melt(data.plot1)
colnames(data.plot1)[(ncol(data.plot1)-1):ncol(data.plot1)] = c('Site','Change.pct')
saveRDS(data.plot1,"data.plot1")


## lesion 1 or 2
data.plot1$Lesion.number = "LS1"
data.plot1$Lesion.number[grep('2[.]pct$', data.plot1$Site)] = 'LS2'
table(data.plot1$Lesion.number)
#LS1 LS2 
#656 574

saveRDS(data.plot1,"data.plot1")
## use line to connect lesions from the same pt instead of calling it site and site 2
#data.plot$Site = gsub('[.]2[.]pct$','.pct',data.plot$Site)
#print(table(data.plot$Site))

data.plot1$Group = paste0(data.plot1$BRR,'.',data.plot1$Lesion.number)


saveRDS(data.plot1,"data.plot1.rds")

##-------------------------------------------------
## for each site, take mean of multiple lesions 

#for some reason adrenal 2 has a weird name so change it to proper
#data.plot1<- lapply(data.plot1, function(x){gsub("adrenal 2.adrenal.pct", "adrenal 2.pct", x)})
data.plot1<-as.data.frame(data.plot1)
data.plot1$Lesion.number[grep('2[.]pct$', data.plot1$Site)] = 'LS2'
#need to fix data plot 2 as it isn't taking average across multiple lesions
data.plot.2 = NULL 
data.plot1<-data.plot1 %>% group_by(PID)
saveRDS(data.plot1,"data.plot1.rds")
#i think the prolem is that theree are two sites. Change all sites to the same and rely on lesion for future separation
data.plot2<-data.frame(data.plot1)
data.plot2[data.plot2=="adrenal 2.pct"]<-"adrenal.pct"
data.plot2[data.plot2=="visceral other 2.pct"]<-"visceral other.pct" 
data.plot2[data.plot2=="soft tissue 2.pct"]<-"soft tissue.pct" 
data.plot2[data.plot2=="nodes 2.pct"]<-"nodes.pct" 
data.plot2[data.plot2=="brain 2.pct"]<-"brain.pct" 
data.plot2[data.plot2=="lung 2.pct"]<-"lung.pct" 
data.plot2[data.plot2=="liver 2.pct"]<-"liver.pct" 
saveRDS(data.plot2, "data.plot.2.rds")

data.plot2<-data.plot2 %>% group_by(Site)
data.plot2<-data.plot2[!is.na(data.plot2$Change.pct),] 
str(data.plot2)
data.plot2$Change.pct<-as.numeric(data.plot2$Change.pct)
#note two ways to perform mean expression and I verified each works (a and b below)
b<-data.plot2 %>% group_by(PID,Site) %>%
  summarise(mean = mean(Change.pct), na.rm=TRUE)
b

a<-aggregate(data.plot2$Change.pct, by=list(PID=data.plot2$PID, Site=data.plot2$Site),data=data.plot2,FUN=mean, na.rm=TRUE)


b1<-merge(data.plot2, b)
saveRDS(b1, "average_percent_change_by_patient_and_Site.rds")
b1<-readRDS("average_percent_change_by_patient_and_Site.rds")
#b1 or average_percent_change_by_patient_and_Site is average while data.plot2 is individual
dim(b1) #[1][1] 236   9
my.levels<-c("PD","SD","PR","CR")
b1$Best.Response.by.RECIST<-factor(b1$Best.Response.by.RECIST, my.levels)
b1$Best.Response.by.RECIST<-na.omit(b1$Best.Response.by.RECIST)

#editing this because Zane wants the facet wrap not to say ".pct"
b1<-b1 %>% mutate(Site.renamed=b1$Site)
b1$Site.renamed<-str_remove(b1$Site.renamed, ".pct")
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')


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


tiff("Zane_and_Alex_boxplot_Change_percent_between_lesion_and_Group.tiff", res=300, width=8, height=6, units='in')
p1
dev.off()

#edited for Zane's talk 9-26-22
b2<-b2 %>% rename(`Lesion Number`=Lesion.number)
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

tiff("Zane_and_Alex_Edited_boxplot_Change_percent_between_lesion_and_Group.tiff", res=300, width=10, height=6, units='in')
p1
dev.off()

pdf("Zane_and_Alex_Edited_boxplot_Change_percent_between_lesion_and_Group.pdf",  width=10, height=6)
p1
dev.off()


##------------------------------------------------
p2 = ggplot(b1, aes(BRR,
                    mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  stat_compare_means(comparisons = list(c('NR','R'))) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  xlab("Overall Best Response by RECIST")+
  ylab("Lesion Response (%)")

p2

tiff("Zane_and_Alex_boxplot_Change_percent_between_only_group.tiff", res=300, width=8, height=6, units='in')
p2
dev.off()


plot.colors2 = c('CR' = '#0000CC','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' ='#CC0000' ,
                 'NA' = '#C0C0C0')
b2<-b1 %>% rename(`Best Response by RECIST`=Best.Response.by.RECIST)
#edited for Zane's presentation
#•	Left and right figure, remove adrenal, visceral, peritoneum
b2<-b2 %>% filter(!Site.renamed %in% c("adrenal", "visceral other", "peritoneum"))
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


tiff("Zane_and_Alex_edited_boxplot_Change_percent_between_only_group.tiff", res=300, width=8, height=6, units='in')
p2
dev.off()



pdf("Zane_and_Alex_edited_boxplot_Change_percent_between_only_group.pdf")
p2
dev.off()



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

tiff("Zane_and_Alex_boxplot_Change_percent_site_and_Site_NR_v_R_.tiff", res=300, width=8, height=6, units='in')
p2.2
dev.off()

#Zane wanted median values for each Site so summarize that here
saveRDS(b1, "average_percent_change_by_patient_and_Site.rds")
d<-b1 %>% group_by(BRR, Site.renamed) %>% summarize(Median=median(Change.pct))
write.csv(d,"Median_for_Zane_and_Alex_boxplot_Change_percent_site_and_Site_NR_v_R_.csv")

#mean values for eosinphil count and neutorphil:leukocyte ratio
#% eos pre-RX
data.plot<-data.plot %>%mutate(Eosinophil_dec=(as.numeric(data.plot$'% eos pre-RX'))/100)
data.plot<-data.plot %>% mutate(Eosinophil_count=(data.plot$Eosinophil_dec*data.plot$`WBC count pre-Rx (x10^3)`))
#keep<-c("% eos pre-RX", "Eosinophil_dec","Eosinophil_count", "WBC count pre-Rx (x10^3)")
#d.1<-data.plot[,keep]
data.plot<-data.plot %>% rename("Eosinophil_count(x10^3)"="Eosinophil_count")
to.remove<-c(78, 24)

data.plot.3<-data.plot[-to.remove,]
e<-mean(data.plot.3$`Eosinophil_count(x10^3)`) #[1] 0.1878788

#neutrophil:leukocyte ratio 
saveRDS(data.plot.3, "data.plot.3.rds")
#%neutr pre-Rx/% lymph pre-Rx
data.plot.3$`%neutr pre-Rx`<-as.numeric(data.plot.3$`%neutr pre-Rx`)
data.plot.3$`% lymph pre-Rx`<-as.numeric(data.plot.3$`% lymph pre-Rx`)
data.plot.3<-data.plot.3 %>% mutate(NLR=(data.plot.3$`%neutr pre-Rx`/data.plot.3$`% lymph pre-Rx`))
keep<-c("PID","%neutr pre-Rx", "% lymph pre-Rx","NLR")
d.2<-data.plot.3[,keep]

write.csv(d.2, "Zane_and_Alex_Neutrophil_to_Lymphocyte.csv")

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

tiff("Zane_and_Alex_boxplot_Change_percent_site_and_Site_BRR_.tiff", res=300, width=8, height=6, units='in')
p3
dev.off()


saveRDS(b1, "average_percent_change_by_patient_and_Site.rds")
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


tiff("Zane_and_Alex_change_direction_by_patient.tiff.tiff", res=300, width=8, height=6, units='in')
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



tiff("Zane_and_Alex_change_direction_by_response.9.26.22.tiff", res=300, width=8, height=6, units='in')
p4.3
dev.off()


pdf("Zane_and_Alex_change_direction_by_response.9.26.22.pdf")
p4.3
dev.off()

write.csv(data.plot.2.3, "Change.direction.by.response.n.group.csv")
#note. data.plot has all original data

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


tiff("Zane_and_Alex_Difference_in_percent_change_between_lesions_Over_Mean.tiff", res=300, width=12, height=6, units='in')
p1
dev.off()


p2<-ggplot(B1.2, aes(x=Site, y=homog.only))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  ylab("Difference in percent change Lesion 1 and Lesion 2 per tissue")+
  stat_compare_means(method="wilcox.test", ref.group = "liver.pct", label = "p.signif") 

tiff("Zane_and_Alex_Difference_in_percent_change_between_lesions.tiff", res=300, width=12, height=6, units='in')
p2
dev.off()
#######################################################################
##Re-run all plots for split in clinical factors
#9-12-22
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
##Mutational status#################################################################
##start with mutational status (BRAF, NF1, NRAS)---------------------------------
#dataplot has all clinical, b1 has averaged change.pct
mutation<- data.plot %>% select("PID","BRAF status", "RAS (NRAS/KRAS)", "NF1", "Other (specify)")

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


WT.NA<-c("N/A_N/A_N/A", "WT_N/A_N/A","WT_WT_N/A","WT_WT_WT")
BRAF<-c("V600_N/A_N/A",  "V600_WT_N/A",   "V600_WT_WT")
NRAS<-c("WT_NRAS_N/A",   "WT_NRAS_WT")
NF1<-c("WT_WT_NF1")


Mutational.df.1$Simplified_mutational_Group=NA
Mutational.df.1$Simplified_mutational_Group[which(Mutational.df.1$Mutational_group %in% WT.NA)] = 'NAorWT' 
Mutational.df.1$Simplified_mutational_Group[which(Mutational.df.1$Mutational_group %in% BRAF)] = 'BRAF' 
Mutational.df.1$Simplified_mutational_Group[which(Mutational.df.1$Mutational_group %in% NRAS)] = 'NRAS' 
Mutational.df.1$Simplified_mutational_Group[which(Mutational.df.1$Mutational_group %in% NF1)] = 'NF1' 

saveRDS(Mutational.df.1, "Mutational.df.contains.only.mean.rds")
Mutational.df.1<-readRDS("Mutational.df.contains.only.mean.rds")

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

pdf("Zane_and_Alex_Mutational_status_vs_response.pdf")
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

respond<-a.2 %>% filter(BRR=="R" & Site.renamed %in% c("brain", "lung","nodes","soft tissue", "visceral other"))
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Simplified_mutational_Group, p.adjust.method="bonf")


b<-rbind(nr.1, respond.1) 

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_mutational_status_vs_response.csv")
######################################################################
##LDH elevated vs low---------------------------------
#Text message from Zane 9-12-22, LDH 171 and above is elevated, below is within normal limit (WNL)
LDH<- data.plot %>% select("PID","LDH_Pre-Rx")

LDH$LDH.levels[LDH$`LDH_Pre-Rx` >170] <- "elevated"
LDH$LDH.levels[LDH$`LDH_Pre-Rx` =="N/A"] <- "N/A"
#for some reason there is one that just has "n" so replace that with N/A
LDH$LDH.levels[LDH$`LDH_Pre-Rx` =="n"] <- "N/A"
LDH$LDH.levels[LDH$`LDH_Pre-Rx` <171] <- "WNL"


LDH.df<-merge(b1, LDH)
saveRDS(LDH.df, "LDH.df.contains.all.replicates.rds")
LDH.df<-readRDS("LDH.df.contains.all.replicates.rds")
LDH.df.1<-filter(LDH.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(LDH.df.1, "LDH.df.contains.only.mean.rds")
LDH.df.1<-readRDS("LDH.df.contains.only.mean.rds")

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

pdf("Zane_and_Alex_LDH_vs_response.pdf")
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

respond<-a.2 %>% filter(BRR=="R" & Site.renamed %in% c("lung","nodes","soft tissue"))
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ LDH.levels, p.adjust.method="bonf")


b<-rbind(nr.1, respond.1) 


write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_LDH_vs_response.csv")
######################################################################
##BMI-----------------------------------------------------------------
data.plot<-data.plot %>% mutate(`height (m)`=`height (cm)`*0.01)
data.plot<-data.plot %>% mutate(BMI=`weight (Kg)`/(`height (m)`)^2)
saveRDS(data.plot,"data.plot.rds")

BMI<- data.plot %>% select("PID","BMI")
BMI$BMI.25[BMI$BMI >25] <- "high"
BMI$BMI.25[BMI$BMI <=25] <- "low"


BMI$BMI.30[BMI$BMI >30] <- "high"
BMI$BMI.30[BMI$BMI <=30] <- "low"


b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
BMI.df<-merge(b1, BMI)

saveRDS(BMI.df, "BMI.df.contains.all.replicates.rds")
BMI.df<-readRDS("BMI.df.contains.all.replicates.rds")
BMI.df.1<-filter(BMI.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(BMI.df.1, "BMI.df.contains.only.mean.rds")
BMI.df.1<-readRDS("BMI.df.contains.only.mean.rds")


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

pdf("Zane_and_Alex_BMI_over_25_vs_response.pdf")
p12
dev.off()

#calculate pairwise significance levels
a<-BMI.df.1 %>% group_by(Site.renamed, BRR) %>% count(BMI.25) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(BMI.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one LDH per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ BMI.25, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_BMI_greater_than_25_vs_response.csv")




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

pdf("Zane_and_Alex_BMI_over_30_vs_response.pdf")
p13
dev.off()


a<-BMI.df.1 %>% group_by(Site.renamed, BRR) %>% count(BMI.30) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(BMI.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one LDH per group


nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ BMI.30, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R" & Site.renamed %in% c("adrenal", "liver","lung","nodes","soft tissue", "visceral other"))
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ BMI.30, p.adjust.method="bonf")


b<-rbind(nr.1, respond.1) 

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_BMI_greater_than_30_vs_response.csv")

#make density plot of BMI
d <- density(BMI.df.1$BMI) # returns the density data 
plot(d)

pdf("Zane_and_Alex_BMI_density.pdf")
plot(d)
dev.off()

length(unique(BMI.df.1$PID)) #number of unique patients(not including multiple tissues per patient)
unique.patient.df<- BMI.df.1[!duplicated(BMI.df.1$PID),]

table(unique.patient.df$BMI.25)
table(unique.patient.df$BMI.30)

#skipping Metabolic syndrome question-------------------------------

#Albumin levels Albumin <3.5
data.plot<-readRDS("data.plot.rds")
data.plot<-data.plot %>% mutate(Albumin.level=ifelse(data.plot$`Albumin pre-rx`< 3.5, "low", "WNL"))
saveRDS(data.plot, "data.plot.rds")

Albumin<- data.plot %>% select("PID","Albumin.level")


b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
Albumin.df<-merge(b1, Albumin)

saveRDS(Albumin.df, "Albumin.df.contains.all.replicates.rds")
Albumin.df.1<-filter(Albumin.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(Albumin.df.1, "Albumin.df.contains.only.mean.rds")


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

pdf("Zane_and_Alex_Albumin_vs_response.pdf")
p14
dev.off()


a<-Albumin.df.1 %>% group_by(Site.renamed, BRR) %>% count(Albumin.level) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Albumin.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one LDH per group

nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Albumin.level, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R" & Site.renamed %in% c("adrenal", "brain","lung","nodes","soft tissue", "visceral other"))
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Albumin.level, p.adjust.method="bonf")


b<-rbind(nr.1, respond.1) 



write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_Albumin_vs_response.csv")
#############################################################
#melanoma subtype
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
data.plot<-readRDS("data.plot.rds")
table(data.plot$`Melanoma subtype`)
# 
# Acral       Cutaneous         Mucosal             NOS unknown primary           Uveal 
#$$$$5              52               8              11               2               4 
#there are 2 unknown primaries that we can exclude

Subtype<- data.plot %>% select("PID","Melanoma subtype")


b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
subtype.df<-merge(b1, Subtype)

#subtype.df<-subtype.df %>% filter (`Melanoma subtype` !="unknown primary")
table(subtype.df$`Melanoma subtype`)

saveRDS(subtype.df, "subtype.df.contains.all.replicates.rds")
subtype.df.1<-filter(subtype.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(subtype.df.1, "subtype.df.contains.only.mean.rds")

#Melanoma subtype
p15 = ggplot(subtype.df.1, aes(`Melanoma subtype`,
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

pdf("Zane_and_Alex_Melanoma_subtype_vs_response.pdf")
p15
dev.off()


subtype.df.1<-subtype.df.1 %>% rename(Melanoma.subtype=`Melanoma subtype`)
a<-subtype.df.1 %>% group_by(Site.renamed, BRR) %>% count(`Melanoma.subtype`) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(subtype.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1)
a.2<- a.2 %>% filter(Melanoma.subtype !="unknown primary")


#remove out samples that only have one  per group
#pairwise wilcox is giving me so much trouble. I'm going to try to break it down to see which samples are the issue
nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Melanoma.subtype, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R", Site.renamed== "nodes")
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Melanoma.subtype, p.adjust.method="bonf")


b<-rbind(nr.1, respond.1) #note, had to do it differently wheere we only look at nodes since that is the only 
#tissue that has R with more than just one sample in comparison group.

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_melanoma_subtype_vs_response.csv")

##compare cutaneous vs all others combined (excluding NOW)--------------------------
subtype.df.1<-readRDS("subtype.df.contains.only.mean.rds")
combined<-c("Acral", "Mucosal", "Uveal")
subtype.df.1$Melanoma.subtype.combined[subtype.df.1$`Melanoma subtype` %in% combined] <- "Acral_Mucosal_Uveal"
subtype.df.1$Melanoma.subtype.combined[subtype.df.1$`Melanoma subtype`=="Cutaneous"] <- "Cutaneous"
subtype.df.1$Melanoma.subtype.combined[subtype.df.1$`Melanoma subtype`=="NOS"] <- "NOS"

saveRDS(subtype.df.1, "subtype.df.contains.only.mean.rds")

#filter out unknown primary
subtype.df.1<-subtype.df.1 %>% filter(`Melanoma subtype` !="unknown primary")
#Melanoma subtype combined
p16 = ggplot(subtype.df.1, aes(Melanoma.subtype.combined,
                               mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Melanoma Subtype Combined")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p16

pdf("Zane_and_Alex_Melanoma_subtype_combined_vs_response.pdf")
p16
dev.off()

a<-subtype.df.1 %>% group_by(Site.renamed, BRR) %>% count(Melanoma.subtype.combined) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(subtype.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one LDH per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ Melanoma.subtype.combined, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_melanoma_subtype_combined_vs_response.csv")
####################--------------------------------------------------
#GENDER
data.plot<-readRDS("data.plot.rds")
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
table(data.plot$Sex)

Gender<- data.plot %>% select("PID","Sex")


b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
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

pdf("Zane_and_Alex_Melanoma_Gender_vs_response.pdf")
p17
dev.off()

a<-gender.df.1 %>% group_by(Site.renamed, BRR) %>% count(Sex) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(gender.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ Sex, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_gender_vs_response.csv")
####################--------------------------------
#Baseline tumor size, assuming this is in cm

plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

data.plot<-readRDS("data.plot.rds")
table(data.plot$`Baseline lesion measurement`)

Tumor.size<- data.plot %>% select("PID","Baseline lesion measurement")
Tumor.size$size[Tumor.size$`Baseline lesion measurement`>5]<-"large"
Tumor.size$size[Tumor.size$`Baseline lesion measurement`<=5]<-"small"

table(Tumor.size$size, useNA="always")

b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
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

pdf("Zane_and_Alex_Melanoma_Baseline_tumor_size_vs_response.pdf")
p18
dev.off()

a<-baseline.size.df.1 %>% group_by(Site.renamed, BRR) %>% count(size) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(baseline.size.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ size, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_baseline_tumor_size_vs_response.csv")

table(baseline.size.df.1$size)
####################---------------------------------------------------

#total WBC
#Baseline tumor size, assuming this is in cm
data.plot<-readRDS("data.plot.rds")
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
table(data.plot$`WBC count pre-Rx (x10^3)`)

WBC.count<- data.plot %>% select("PID","WBC count pre-Rx (x10^3)")

#4.5-11 is normal range 
WBC.count$WBC_group[WBC.count$`WBC count pre-Rx (x10^3)`<4.5]<-"low"
WBC.count$WBC_group[WBC.count$`WBC count pre-Rx (x10^3)`>11]<-"high"
WBC.count$WBC_group[WBC.count$`WBC count pre-Rx (x10^3)`>=4.5 &WBC.count$`WBC count pre-Rx (x10^3)`<=11 ]<-"normal"


table(WBC.count$WBC_group, useNA="always") #no NAs
b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
WBC.count.df<-merge(b1, WBC.count)

table(WBC.count.df$WBC_group)

saveRDS(WBC.count.df, "WBC.count.df.contains.all.replicates.rds")
WBC.count.df.1<-filter(WBC.count.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(WBC.count.df.1, "WBC.count.df.contains.only.mean.rds")

#Melanoma subtype
p19 = ggplot(WBC.count.df.1, aes(WBC_group,
                                 mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("WBC Group")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p19

pdf("Zane_and_Alex_Melanoma_WBC_group_vs_response.pdf")
p19
dev.off()

a<-WBC.count.df.1 %>% group_by(Site.renamed, BRR) %>% count(WBC_group) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(WBC.count.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ WBC_group, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_WBC_group_vs_response.csv")

####################--------------------------------
#neutrophil lymphocyte ratio (from above)
data.plot<-readRDS("data.plot.rds")
d.2<-read.csv("Zane_and_Alex_Neutrophil_to_Lymphocyte.csv")
d.2<-d.2[,-1]
d.2<-d.2[,-c(2,3)]
d.2$NLR_group[d.2$NLR>3]<-"high"
d.2$NLR_group[d.2$NLR<=3]<-"normal"

data.plot<-merge(d.2, data.plot, by="PID")
saveRDS(data.plot, "data.plot.rds")

plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
table(data.plot$NLR)

NLR<- data.plot %>% select("PID","NLR_group")

table(NLR$NLR_group, useNA="always")
b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
NLR.df<-merge(b1, NLR)

table(NLR.df$NLR_group) #140 high 90 normal

saveRDS(NLR.df, "NLR.df.contains.all.replicates.rds")
NLR.df.1<-filter(NLR.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(NLR.df.1, "NLR.df.contains.only.mean.rds")

#Melanoma subtype
p20 = ggplot(NLR.df.1, aes(NLR_group,
                           mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Neutrophil:Lymphocyte Ratio")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p20

pdf("Zane_and_Alex_Melanoma_NLR_vs_response.pdf")
p20
dev.off()

a<-NLR.df.1 %>% group_by(Site.renamed, BRR) %>% count(NLR_group) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(NLR.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ NLR_group, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_NLR_vs_response.csv")
###--------------------------------------------------------------------
##NLR split by median
#neutrophil lymphocyte ratio (from above)
data.plot<-readRDS("data.plot.rds")

d.2<-read.csv("Zane_and_Alex_Neutrophil_to_Lymphocyte.csv")
d.2<-d.2[,-1]
d.2<-d.2[,-c(2,3)]
d.2$NLR_group_median[d.2$NLR>median(d.2$NLR)]<-"high"
d.2$NLR_group_median[d.2$NLR<=median(d.2$NLR)]<-"low"
d.2<-d.2[,-2]

data.plot<-merge(data.plot, d.2, by="PID")
NLR_median<-data.plot %>% select("PID", "NLR_group_median")

table(NLR_median$NLR_group_median, useNA="always")
b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
NLR_median.df<-merge(b1, NLR_median)

table(NLR_median.df$NLR_group) #122 high 108 low

saveRDS(NLR_median.df, "NLR_median.df.contains.all.replicates.rds")
NLR_median.df.1<-filter(NLR_median.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(NLR_median.df.1, "NLR_median.df.contains.only.mean.rds")

#Melanoma subtype
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

pdf("Zane_and_Alex_Melanoma_NLR_split_by_median_vs_response.pdf")
p21
dev.off()

a<-NLR_median.df.1 %>% group_by(Site.renamed, BRR) %>% count(NLR_group_median) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(NLR_median.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ NLR_group_median, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_NLR_split_by_median_vs_response.csv")



####################--------------------------------
#eosinophil count 
data.plot.3<-readRDS("data.plot.3.rds")
a<-data.plot.3 %>% select("PID", "Eosinophil_count(x10^3)")
a<-a %>% mutate(Eosinophil_count=(a$"Eosinophil_count(x10^3)"*1000))
mean(a$Eosinophil_count)
data.plot<-readRDS('data.plot.rds')
data.plot<-merge(data.plot, a, by="PID")
saveRDS(data.plot, "data.plot.rds")


plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')
table(data.plot$Eosinophil_count.y)
data.plot<-data.plot %>% rename(Eosinophil_count=Eosinophil_count.y)
saveRDS(data.plot, "data.plot.rds")

Eosinophil<- data.plot %>% select("PID","Eosinophil_count")
median(Eosinophil$Eosinophil_count) #117.6


Eosinophil$Eosinophil_group[Eosinophil$Eosinophil_count>median(Eosinophil$Eosinophil_count)]<-"high"
Eosinophil$Eosinophil_group[Eosinophil$Eosinophil_count<=median(Eosinophil$Eosinophil_count)]<-"low"

table(Eosinophil$Eosinophil_group, useNA="always")
b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
Eosinophil.df<-merge(b1, Eosinophil)

table(Eosinophil.df$Eosinophil_group) #110 high 120 normal

saveRDS(Eosinophil.df, "Eosinophil.df.contains.all.replicates.rds")
Eosinophil.df.1<-filter(Eosinophil.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(Eosinophil.df.1, "Eosinophil.df.contains.only.mean.rds")
Eosinophil.df.1<-readRDS("Eosinophil.df.contains.only.mean.rds")
#eosinophil
p21 = ggplot(Eosinophil.df.1, aes(Eosinophil_group,
                                  mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Eosinophil (split by median)")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p21

pdf("Zane_and_Alex_Melanoma_Eosinophil_count_vs_response.pdf")
p21
dev.off()

a<-Eosinophil.df.1 %>% group_by(Site.renamed, BRR) %>% count(Eosinophil_group) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Eosinophil.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ Eosinophil_group, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_Eosinophil_count_vs_response.csv")


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

pdf("Zane_and_Alex_Melanoma_Eosinophil_count_not_separated_by_BRR_vs_response.pdf")
p21.1
dev.off()

a<-Eosinophil.df.1 %>% group_by(Site.renamed) %>% count(Eosinophil_group) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Eosinophil.df.1, a.1, by=c("Site.renamed"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed)%>% pairwise_wilcox_test(mean ~ Eosinophil_group, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_Eosinophil_count_not_separated_by_BRR_vs_response.csv")






##sun exposure---------------------------------------
data.plot<-readRDS("data.plot.rds")


#melanoma subtype
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

# 
#Acral       Cutaneous         Mucosal        NOS         unknown primary           Uveal 
#9             153              20              32               7              15 
#there are 7 unknown primaries that we can exclude

Sun<- data.plot %>% select("PID","extent of sun exposure (If mentioned in chart)")


b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
Sun.df<-merge(b1, Sun)

table(Sun.df$`extent of sun exposure (If mentioned in chart)`)
#High  Low  N/A   NA 
#124   30   71    5 

saveRDS(Sun.df, "Sun.df.contains.all.replicates.rds")
Sun.df.1<-filter(Sun.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(Sun.df.1, "Sun.df.contains.only.mean.rds")
Sun.df.1<-Sun.df.1 %>% filter(`extent of sun exposure (If mentioned in chart)`!="NA")
Sun.df.1<-Sun.df.1 %>% filter(`extent of sun exposure (If mentioned in chart)`!="N/A")

#sun
p22 = ggplot(Sun.df.1, aes(`extent of sun exposure (If mentioned in chart)`,
                           mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Extent of Sun Exposure (if mentioned in chart)")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p22

pdf("Zane_and_Alex_Melanoma_Sun_exposure_vs_response.pdf")
p22
dev.off()


a<-Sun.df.1 %>% group_by(Site.renamed, BRR) %>% count(`extent of sun exposure (If mentioned in chart)`) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Sun.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

a.2<-a.2 %>% rename(Sun.exposure=`extent of sun exposure (If mentioned in chart)`)
#b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ Sun.exposure, p.adjust.method="bonf")
#giing trouble again try with just nr
nr<-a.2 %>% filter(BRR=="NR")
nr<-nr %>% filter(n>1)
nr.1<-nr %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Sun.exposure, p.adjust.method="bonf")

respond<-a.2 %>% filter(BRR=="R", Site.renamed %in% c("brain", "lung", "nodes", "soft tissue"))
respond<-respond %>% filter(n>1)
respond.1<-respond %>% group_by(Site.renamed, BRR) %>% pairwise_wilcox_test(mean ~ Sun.exposure, p.adjust.method="bonf")

b<-rbind(nr.1, respond.1) #note, had to do it differently wheere we only look at nodes since that is the only 
#tissue that has R with more than just one sample in comparison group.

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_sun_exposure_vs_response.csv")

###################################################################
##progression on adjuvant therapy (pd1 or braf)

data.plot<-readRDS("data.plot.rds")


#
plot.colors = c('CR' = '#CC0000','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' = '#0000CC',
                'NA' = '#C0C0C0')

# 

Prog.on.adjuvant<- data.plot %>% select("PID","Prior adjuvant immunotherapy  treatment")

b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
Prog.on.adjuvant.df<-merge(b1, Prog.on.adjuvant)

#subtype.df<-subtype.df %>% filter (`Melanoma subtype` !="unknown primary")
table(Prog.on.adjuvant.df$`Prior adjuvant immunotherapy  treatment`)
#High  Low  N/A   NA 
#124   30   71    5 

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

pdf("Zane_and_Alex_Melanoma_prog_on_adjuv_therapy_vs_response.pdf")
p23
dev.off()


a<-Prog.on.adjuvant.df.1 %>% group_by(Site.renamed, BRR) %>% count(Prog.on.therapy) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Prog.on.adjuvant.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ Prog.on.therapy, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_prog_on_adjuvant_count_vs_response.csv")

########################################################################
#line of metastatic immunotherapy (need to upload from other excel sheet.)


Line.of.tx<-readxl::read_excel('Melanoma_IpiNivo_THIS_IS_THE_MOST_UPTODATE.xlsx', 
                               sheet = "Line of Tx")

Line.of.tx<-Line.of.tx %>% select("Patient identifier", "Line of metastatic immunotherapy therapy IpiNivo")
Line.of.tx$Line.of.therapy.combiend[Line.of.tx$"Line of metastatic immunotherapy therapy IpiNivo"==1]<-1
Line.of.tx$Line.of.therapy.combiend[Line.of.tx$"Line of metastatic immunotherapy therapy IpiNivo"==2]<-"2_and_3"
Line.of.tx$Line.of.therapy.combiend[Line.of.tx$"Line of metastatic immunotherapy therapy IpiNivo"==3]<-"2_and_3"
Line.of.tx$Line.of.therapy.combiend[Line.of.tx$"Line of metastatic immunotherapy therapy IpiNivo"=="2nd"]<-"2_and_3"
Line.of.tx$Line.of.therapy.combiend<-as.character(Line.of.tx$Line.of.therapy.combiend)
str(Line.of.tx)

Line.of.tx<-Line.of.tx[,-2]
data.plot<-merge(data.plot, Line.of.tx, by="Patient identifier")
saveRDS(data.plot, "data.plot.rds")

Line.of.therapy<- data.plot %>% select("PID","Line.of.therapy.combiend")

b1<-read.csv("average_percent_change_by_patient_and_Site.csv")
b1<-b1[,-1]
Line.of.therapy.df<-merge(b1, Line.of.therapy)

#subtype.df<-subtype.df %>% filter (`Melanoma subtype` !="unknown primary")
table(Line.of.therapy.df$Line.of.therapy.combiend)
#  1      2_and_3    <NA> 
### 157      79 

saveRDS(Line.of.therapy.df, "Line.of.therapy.df.contains.all.replicates.rds")
Line.of.therapy.df.1<-filter(Line.of.therapy.df, grepl("LS1", Lesion.number)) #have to get rid of repeating lesion. NOTE that we have gotten rid of second lesion and 
#we MUST use the mean values now
saveRDS(Line.of.therapy.df.1, "Line.of.therapy.df.contains.only.mean.rds")
table(Line.of.therapy.df.1$Line.of.therapy.combiend, useNA="always")

saveRDS(Line.of.therapy.df.1, "Line.of.therapy.df.contains.only.mean.rds")
Line.of.therapy.df.1<-readRDS("Line.of.therapy.df.contains.only.mean.rds")

#line of therapy
p24 = ggplot(Line.of.therapy.df.1, aes(Line.of.therapy.combiend,
                                       mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors2) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(BRR~ Site.renamed)+
  xlab("Line of therapy (combined 2 and 3)")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p24

pdf("Zane_and_Alex_Melanoma_line_of_therapy_vs_response.pdf")
p24
dev.off()



a<-Line.of.therapy.df.1 %>% group_by(Site.renamed, BRR) %>% count(Line.of.therapy.combiend) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Line.of.therapy.df.1, a.1, by=c("Site.renamed", "BRR"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed, BRR)%>% pairwise_wilcox_test(mean ~ Line.of.therapy.combiend, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_line_on_therapy_vs_response.csv")


######edited for zane's presentaiton

plot.colors2 = c('CR' = '#0000CC','PR' = '#CC00CC', 'SD' = '#00CCCC','PD' ='#CC0000' ,
                 'NA' = '#C0C0C0')

Line.of.therapy.df.1<-readRDS("Line.of.therapy.df.contains.only.mean.rds")
Line.of.therapy.df.1<-Line.of.therapy.df.1 %>% rename(`Best Response by RECIST`=Best.Response.by.RECIST)
Line.of.therapy.df.1$Line.of.therapy.combiend[Line.of.therapy.df.1$Line.of.therapy.combiend=="2_and_3"]<-"≥2"

Line.of.therapy.df.1$Line.of.therapy.combiend<-factor(Line.of.therapy.df.1$Line.of.therapy.combiend, levels=c("1", "≥2"))
Line.of.therapy.df.1<-Line.of.therapy.df.1 %>% filter(!Site.renamed %in% c("adrenal", "visceral other", "peritoneum"))
#line of therapy
p24 = ggplot(Line.of.therapy.df.1, aes(Line.of.therapy.combiend,
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


#split into response
Response<-Line.of.therapy.df.1 %>% filter(BRR=="R")


p25 = ggplot(Response, aes(Line.of.therapy.combiend,
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


pdf("Zane_and_Alex_Melanoma_RESPONDERS.line_of_therapy_vs_response.pdf")
p25
dev.off()
b<-table(Response$Site.renamed, Response$Line.of.therapy.combiend)
b<-as.data.frame(b)
b$Var2<-as.character(b$Var2)
write.csv(b, "Responders.Line.of.response.numbers.csv")

#NR
N.R<-Line.of.therapy.df.1 %>% filter(BRR=="NR")


p26 = ggplot(N.R, aes(Line.of.therapy.combiend,
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


pdf("Zane_and_Alex_Melanoma_NR.line_of_therapy_vs_response.pdf")
p26
dev.off()

a<-table(N.R$Site.renamed, N.R$Line.of.therapy.combiend)
a<-as.data.frame(a)
a$Var2<-as.character(a$Var2)
write.csv(a, "N.R.Line.of.response.numbers.csv")

#Plot again, Zane wants without stratified by response 9-15-22--------------------------------


#line of therapy
p24.1 = ggplot(Line.of.therapy.df.1, aes(Line.of.therapy.combiend,
                                         mean)) +
  geom_boxplot(color = '#C0C0C0', width = 0.6, outlier.shape = NA,lwd = 1) +
  geom_jitter(aes(color = Best.Response.by.RECIST), width = 0.1, height = 0) +
  #   geom_line(aes(group = PID)) +
  scale_color_manual(values = plot.colors) +
  theme_pubr(x.text.angle = 90) +
  geom_hline(yintercept = -30, linetype = 'dashed') +
  facet_grid(~ Site.renamed)+
  xlab("Line of therapy (combined 2 and 3)")+
  ylab("Percent Change in Tumor Growth of averaged lesions")

p24.1

pdf("Zane_and_Alex_Melanoma_line_of_therapy_not_separated_by_BRR_vs_response.pdf")
p24.1
dev.off()


a<-Line.of.therapy.df.1 %>% group_by(Site.renamed) %>% count(Line.of.therapy.combiend) #fix this
a.1<-a %>% count(Site.renamed)
a.2<-merge(Line.of.therapy.df.1, a.1, by=c("Site.renamed"))
a.2<-a.2 %>% filter(n>1) #remove out samples that only have one per group

b<-a.2 %>% group_by(Site.renamed)%>% pairwise_wilcox_test(mean ~ Line.of.therapy.combiend, p.adjust.method="bonf")

write.csv(b, "Zane_and_Alex_Pairwise_wilcox_test_for_line_on_therapy_not_separated_by_BRR_vs_response.csv")




###########################################################################################

##STOPPED HERE TO START RADIOMICS MACHINE LEARNING. See this in other file as well.
##################################################

##Input Radiomic Features data sets #first order
getwd()
setwd("/ix/jluke/Projects/becky/Radiomics_Zane_Alex/IpiNivo Radiomic Features/CT_Features/First_Order")
#need to data wrangle so we can properly label each file with what they are
#STL
First.Order.CT.STL<-read.csv("First_Order_STL_scan_segment.nii.csv")
colnames(First.Order.CT.STL)<-paste("First.Order.CT.STL", colnames(First.Order.CT.STL), sep="_")
First.Order.CT.STL<- First.Order.CT.STL %>% rename(ID=First.Order.CT.STL_ID)
#First.Order.CT.STL<- First.Order.CT.STL[,-1]
write.csv(First.Order.CT.STL, "First_Order_CT_STL_scan_segment.nii.csv",row.names=FALSE)
First.Order.CT.STL<-First.Order.CT.STL[,3:14]
colnames(First.Order.CT.STL)<-gsub("First.Order.CT.STL_First.Order.CT.STL","First.Order.CT.STL",colnames(First.Order.CT.STL))

First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-008-1","RRC-GrayLuke-008",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-019-1","RRC-GrayLuke-019",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-028-2","RRC-GrayLuke-028",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-056-1","RRC-GrayLuke-056",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-083-1","RRC-GrayLuke-083",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-011-1","RRC-GrayLuke-011",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-011-2","RRC-GrayLuke-011",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-019-1","RRC-GrayLuke-019",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-028-1","RRC-GrayLuke-028",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("RRC-GrayLuke-053-1","RRC-GrayLuke-053",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL$ID<-gsub("-2","",as.character(First.Order.CT.STL$ID))
First.Order.CT.STL.1<-First.Order.CT.STL %>% group_by(ID) %>% summarise(across(everything(), list(mean)))
#note patient #11 and #56 has two some i'm just going to take averages across feats for that patient above

write.csv(First.Order.CT.STL, "First_Order_CT_STL_scan_segment.nii.csv",row.names=FALSE)

a<-c(122.807655,99.402530)
mean(a)

#SL

First.Order.CT.SL<-read.csv("First_Order_CT_SL_scan_segment.nii.csv")
colnames(First.Order.CT.SL)<-paste("First.Order.CT.SL", colnames(First.Order.CT.SL), sep="_")
First.Order.CT.SL<- First.Order.CT.SL %>% rename(ID=First.Order.CT.SL_First.Order.CT.SL_ID)
First.Order.CT.SL<-First.Order.CT.SL[,3:14]
colnames(First.Order.CT.SL)<-gsub("First.Order.CT.SL_First.Order.CT.SL", "First.Order.CT.SL", colnames(First.Order.CT.SL))
write.csv(First.Order.CT.SL, "First_Order_CT_SL_scan_segment.nii.csv", row.names=F)
First.Order.CT.SL<-First.Order.CT.SL[,-1]
First.Order.CT.SL$ID<-gsub("-1","",as.character(First.Order.CT.SL$ID))

#Lung

First.Order.CT.Lung<-read.csv("First_Order_CT_Lung_L_scan_segment.nii.csv")
colnames(First.Order.CT.Lung)<-paste("First.Order.CT.Lung", colnames(First.Order.CT.Lung), sep="_")
First.Order.CT.Lung<- First.Order.CT.Lung %>% rename(ID=First.Order.CT.Lung_ID)
First.Order.CT.Lung<-First.Order.CT.Lung[,3:14]
First.Order.CT.Lung<-First.Order.CT.Lung[,-1]
colnames(First.Order.CT.Lung)<-gsub("First.Order.CT.Lung_First.Order.CT.Lung","First.Order.CT.Lung",colnames(First.Order.CT.Lung))
#write.csv(First.Order.CT.Lung, "First_Order_CT_Lung_L_scan_segment.nii.csv", row.names=F)


First.Order.CT.LN<-read.csv("First_Order_CT_LN_scan_segment.nii.csv")
colnames(First.Order.CT.LN)<-paste("First.Order.CT.LN", colnames(First.Order.CT.LN), sep="_")
First.Order.CT.LN<- First.Order.CT.LN %>% rename(ID=First.Order.CT.LN_ID)
First.Order.CT.LN<-First.Order.CT.LN[,-1]
First.Order.CT.LN$ID<-gsub("-2","",as.character(First.Order.CT.LN$ID))
First.Order.CT.LN$ID<-gsub("RRC-GrayLuke-019-1","RRC-GrayLuke-019",as.character(First.Order.CT.LN$ID))
First.Order.CT.LN$ID<-gsub("RRC-GrayLuke-028-1","RRC-GrayLuke-028",as.character(First.Order.CT.LN$ID))
First.Order.CT.LN$ID<-gsub("RRC-GrayLuke-083-1","RRC-GrayLuke-083",as.character(First.Order.CT.LN$ID))
First.Order.CT.LN$ID<-gsub("RRC-GrayLuke-087-1","RRC-GrayLuke-087",as.character(First.Order.CT.LN$ID))


#write.csv(First.Order.CT.LN, "First_Order_CT_LN_scan_segment.nii.csv", row.names=F)

First.Order.CT.Liver<-read.csv("First_Order_CT_Liver_L_scan_segment.nii.csv")
colnames(First.Order.CT.Liver)<-paste("First.Order.CT.Liver", colnames(First.Order.CT.Liver), sep="_")
First.Order.CT.Liver<- First.Order.CT.Liver %>% rename(ID=First.Order.CT.Liver_ID)
First.Order.CT.Liver<-First.Order.CT.Liver[,-1]

First.Order.CT.Liver$ID<-gsub("RRC-GrayLuke-028-1","RRC-GrayLuke-028",as.character(First.Order.CT.Liver$ID))
First.Order.CT.Liver$ID<-gsub("RRC-GrayLuke-053-1","RRC-GrayLuke-053",as.character(First.Order.CT.Liver$ID))
First.Order.CT.Liver$ID<-gsub("RRC-GrayLuke-087-1","RRC-GrayLuke-087",as.character(First.Order.CT.Liver$ID))

write.csv(First.Order.CT.Liver, "First_Order_CT_Liver_L_scan_segment.nii.csv", row.names=F)


First.Order.CT.AL<-read.csv("First_Order_CT_AL_scan_segment.nii.csv")
colnames(First.Order.CT.AL)<-paste("First.Order.CT.AL", colnames(First.Order.CT.AL), sep="_")
First.Order.CT.AL<- First.Order.CT.AL %>% rename(ID=First.Order.CT.AL_ID)
First.Order.CT.AL<-First.Order.CT.AL[,-1]
First.Order.CT.AL$ID<-gsub("RRC-GrayLuke-008-1","RRC-GrayLuke-008",as.character(First.Order.CT.AL$ID))
First.Order.CT.AL$ID<-gsub("RRC-GrayLuke-019-1","RRC-GrayLuke-019",as.character(First.Order.CT.AL$ID))

#write.csv(First.Order.CT.AL,"First_Order_CT_AL_scan_segment.nii.csv", row.names=F)

getwd()
setwd("/ix/jluke/Projects/becky/Radiomics_Zane_Alex/IpiNivo Radiomic Features/CT_Features/Second_Order")
Second.Order.CT.AL<-readxl::read_excel('AL_CT_features.xlsx', sheet = 1)
Second.Order.CT.AL<-Second.Order.CT.AL %>% rename(ID="...1")
colnames(Second.Order.CT.AL)<-paste("Second.Order.CT.AL", colnames(Second.Order.CT.AL), sep="_")
Second.Order.CT.AL<-Second.Order.CT.AL %>% rename(ID="Second.Order.CT.AL_ID")
write.xlsx(Second.Order.CT.AL,"Second_Order_AL_CT_features.xlsx")
colnames(Second.Order.CT.AL)<-gsub("Second.Order.CT.AL_Second.Order.CT.Lung_","Second.Order.CT.AL_",colnames(Second.Order.CT.AL))
Second.Order.CT.AL$ID<-gsub("RRC-GrayLuke-008-1","RRC-GrayLuke-008",as.character(Second.Order.CT.AL$ID))
Second.Order.CT.AL$ID<-gsub("RRC-GrayLuke-019-1","RRC-GrayLuke-019",as.character(Second.Order.CT.AL$ID))
Second.Order.CT.AL<-Second.Order.CT.AL %>% column_to_rownames(var="ID")

Second.Order.CT.AL<-Second.Order.CT.AL[,-1] 


#read in all feats and combine to be one big file.. Note that xgboost will be able to handle missing vals
setwd("/ix/jluke/Projects/becky/Radiomics_Zane_Alex/IpiNivo Radiomic Features/CT_Features/Second_Order")

Second.Order.CT.Liver<-readxl::read_excel('Second_Order_Liver_L_CT_features.xlsx', sheet = 1)
colnames(Second.Order.CT.Liver)<-paste("Second.Order.CT.Liver", colnames(Second.Order.CT.Liver), sep="_")
Second.Order.CT.Liver<-Second.Order.CT.Liver %>% rename(ID="Second.Order.CT.Liver_...1")
Second.Order.CT.Liver<-Second.Order.CT.Liver %>% rename(ID="...1")
Second.Order.CT.Liver$ID<-gsub("RRC-GrayLuke-028-1","RRC-GrayLuke-028",as.character(Second.Order.CT.Liver$ID))
Second.Order.CT.Liver$ID<-gsub("RRC-GrayLuke-053-1","RRC-GrayLuke-053",as.character(Second.Order.CT.Liver$ID))
Second.Order.CT.Liver$ID<-gsub("RRC-GrayLuke-087-1","RRC-GrayLuke-087",as.character(Second.Order.CT.Liver$ID))


Second.Order.CT.LN<-readxl::read_excel('Second_Order_LN_CT_features.xlsx', sheet = 1)
colnames(Second.Order.CT.LN)<-paste("Second.Order.CT.LN", colnames(Second.Order.CT.LN), sep="_")
Second.Order.CT.LN<-Second.Order.CT.LN %>% rename(ID="Second.Order.CT.LN_...1")
Second.Order.CT.LN<-Second.Order.CT.LN %>% rename(ID="...1")
Second.Order.CT.LN$ID<-gsub("RRC-GrayLuke-008-2","RRC-GrayLuke-008",as.character(Second.Order.CT.LN$ID))
Second.Order.CT.LN$ID<-gsub("RRC-GrayLuke-019-1","RRC-GrayLuke-019",as.character(Second.Order.CT.LN$ID))
Second.Order.CT.LN$ID<-gsub("RRC-GrayLuke-028-1","RRC-GrayLuke-028",as.character(Second.Order.CT.LN$ID))
Second.Order.CT.LN$ID<-gsub("RRC-GrayLuke-083-1","RRC-GrayLuke-083",as.character(Second.Order.CT.LN$ID))
Second.Order.CT.LN$ID<-gsub("RRC-GrayLuke-084-2","RRC-GrayLuke-084",as.character(Second.Order.CT.LN$ID))
Second.Order.CT.LN$ID<-gsub("RRC-GrayLuke-087-1","RRC-GrayLuke-087",as.character(Second.Order.CT.LN$ID))


Second.Order.CT.Lung<-readxl::read_excel('Second_Order_Lung_L_CT_features.xlsx', sheet = 1)
colnames(Second.Order.CT.Lung)<-paste("Second.Order.CT.Lung", colnames(Second.Order.CT.Lung), sep="_")
Second.Order.CT.Lung<-Second.Order.CT.Lung %>% rename(ID="Second.Order.CT.Lung_...1")
#Second.Order.CT.Lung<-Second.Order.CT.Lung %>% rename(ID="...1")
Second.Order.CT.Lung$ID<-gsub("RRC-GrayLuke-008-1","RRC-GrayLuke-008",as.character(Second.Order.CT.Lung$ID))
Second.Order.CT.Lung$ID<-gsub("RRC-GrayLuke-019-1","RRC-GrayLuke-019",as.character(Second.Order.CT.Lung$ID))
Second.Order.CT.Lung$ID<-gsub("RRC-GrayLuke-028-2","RRC-GrayLuke-028",as.character(Second.Order.CT.Lung$ID))
Second.Order.CT.Lung$ID<-gsub("RRC-GrayLuke-056-1","RRC-GrayLuke-056",as.character(Second.Order.CT.Lung$ID))
Second.Order.CT.Lung$ID<-gsub("RRC-GrayLuke-083-1","RRC-GrayLuke-083",as.character(Second.Order.CT.Lung$ID))
Second.Order.CT.Lung$ID<-gsub("RRC-GrayLuke-084-1","RRC-GrayLuke-084",as.character(Second.Order.CT.Lung$ID))

Second.Order.CT.SL<-readxl::read_excel('Second_Order_SL_CT_features.xlsx', sheet = 1)
colnames(Second.Order.CT.SL)<-paste("Second.Order.CT.SL", colnames(Second.Order.CT.SL), sep="_")
Second.Order.CT.SL<-Second.Order.CT.SL %>% rename(ID="Second.Order.CT.SL_...1")
Second.Order.CT.SL<-Second.Order.CT.SL %>% rename(ID="...1")
Second.Order.CT.SL$ID<-gsub("RRC-GrayLuke-019-1","RRC-GrayLuke-019",as.character(Second.Order.CT.SL$ID))

Second.Order.CT.STL<-readxl::read_excel('Second_Order_STL_CT_features.xlsx', sheet = 1)
colnames(Second.Order.CT.STL)<-paste("Second.Order.CT.STL", colnames(Second.Order.CT.STL), sep="_")
Second.Order.CT.STL<-Second.Order.CT.STL %>% rename(ID="Second.Order.CT.STL_...1")
Second.Order.CT.STL$ID<-gsub("-2","",as.character(Second.Order.CT.STL$ID))
Second.Order.CT.STL$ID<-gsub("RRC-GrayLuke-008-1","RRC-GrayLuke-008",as.character(Second.Order.CT.STL$ID))
Second.Order.CT.STL$ID<-gsub("RRC-GrayLuke-011-1","RRC-GrayLuke-011",as.character(Second.Order.CT.STL$ID))
Second.Order.CT.STL$ID<-gsub("RRC-GrayLuke-028-1","RRC-GrayLuke-028",as.character(Second.Order.CT.STL$ID))
Second.Order.CT.STL$ID<-gsub("RRC-GrayLuke-053-1","RRC-GrayLuke-053",as.character(Second.Order.CT.STL$ID))
Second.Order.CT.STL$ID<-gsub("RRC-GrayLuke-056-1","RRC-GrayLuke-056",as.character(Second.Order.CT.STL$ID))
Second.Order.CT.STL.1<-Second.Order.CT.STL %>% group_by(ID) %>% summarise(across(everything(), list(mean)))

#Note. I found that the serafettin is the most updated response categories. There are two missing samples in the Melanoma 
#clinical spreadsheet so use serafettin for it
#response<-read.csv("response.csv")
#response<-response[,-1]
#Clinical<-Clinical %>% rename(Patient.identifier=`Patient identifier`)
#response<-merge(Clinical, response, by="Patient.identifier")

#colnames(Clinical)
#colnames(response)

#write all clinical data to one file
#write.csv(response, "combined_response_clinical.csv")


