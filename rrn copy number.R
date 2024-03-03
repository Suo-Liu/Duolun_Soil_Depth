setwd("C:/Users/True/OneDrive/桌面")

com.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"
clas.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\16S\\all samples\\Unoise\\taxonomy_16S.txt"
prefix.m = "Bacteria"

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"
env.data.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\03242023 Inner Mongolia TS\\env\\env.csv"
env.data = read.csv(env.data.file,row.names = 1, header = T, sep = ",")

library(ieggr)
comm=t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
name.row = rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] = "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] = "YD3L3"
rownames(comm) = name.row

clas=read.table(clas.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE)
treat = read.csv(treat.file, row.names = 1, header = T,sep = ",")

treat = subset(treat, plant.type == "TS")
sampc=match.name(rn.list = list(treat = treat,comm = comm))
dim(comm)
comm = sampc$comm
treat = sampc$treat
comm = comm[,colSums(comm)>0]
dim(comm)

spc=match.name(rn.list=list(clas=clas), cn.list = list(comm = comm))
comm=spc$comm
clas=spc$clas

# every sample get a average copy number
# weighted ####
adjusted.comm = t(comm)/(clas$copy)
adjusted.comm[1:5,1:5]

table(rownames(comm) == colnames(adjusted.comm))
w.sample.aver.number = rowSums(comm)/colSums(adjusted.comm)
# unweighted ####
adjusted.comm = t(comm>0)/(clas$copy)
adjusted.comm[1:5,1:5]

table(rownames(comm>0) == colnames(adjusted.comm))
uw.sample.aver.number = rowSums(comm>0)/colSums(adjusted.comm)

# select one
# sample.aver.number = w.sample.aver.number
sample.aver.number = uw.sample.aver.number

sample.aver.number = as.data.frame(sample.aver.number)
sampc=match.name(rn.list = list(treat = treat,sample.aver.number = sample.aver.number))
sample.aver.number = sampc$sample.aver.number
treat = sampc$treat

# LMM ####
library(ieggr)
library(lme4)
library(car)
library(tidyr)
library(dplyr)
library(magrittr)
se <- function(x){
  sd(x) / sqrt(length(x))
}

names(sample.aver.number)[1] = "copy" 

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"
treat<-read.csv(treat.file,header = T,row.names = 1)
treat = subset(treat, plant.type == "TS")

used.data <- sample.aver.number

sames = match.name(rn.list = list(used.data = used.data, treat = treat))
used.data = sames$used.data
treat = sames$treat

used.data$treat1 = treat$combined_treat1[match(rownames(used.data), rownames(treat))]
used.data$block = treat$block[match(rownames(used.data), rownames(treat))]
## All layers ####
treat$precip = ifelse(treat$precip!=0,1,0)
div.table = data.frame(stringsAsFactors = F)
mean.table = data.frame(stringsAsFactors = F)
for (i in c("RP","EP","W")){
  divindex <- used.data
  divindex = subset(divindex, treat1%in% c(i,"C"))
  # divindex = subset(divindex, treat1%in% c("RP","EP","W","C"))
  {
    if(i %in% c("RP","EP")){
      divindex$num_treat = treat$precip[match(divindex$treat1, treat$combined_treat1)]
    } else {
      divindex$num_treat = treat$warm[match(divindex$treat1, treat$combined_treat1)]
    }
  }
  divindex$layer = treat$layer[match(rownames(divindex), rownames(treat))]
  
  used.treat = divindex[,2:5]
  divindex = divindex[,-c(2:5), drop = F]
  
  used.treat$layer = (used.treat$layer-min(used.treat$layer))/(max(used.treat$layer)-min(used.treat$layer))
  # scale the alpha diversities
  # divindex<-scale(divindex)
  divs1<-sapply(1:ncol(divindex),function(j){
    message("Now j=",j," in ",ncol(divindex),". ",date())
    if (length(unique(divindex[,j]))<3){
      result<-rep(NA,38)
    } else {
      div<-data.frame(divtest=divindex[,j],used.treat)
      fm<-lmer(divtest~num_treat*layer+(1|block),data=div)
      presult<-car::Anova(fm,type=2)
      coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
      names(coefs)<-paste0(names(coefs),".mean")
      SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
      names(SEvalues)<-paste0(names(SEvalues),".se")
      tvalues<-coef(summary(fm))[ , "t value"]#t values
      names(tvalues)<-paste0(names(tvalues),".t")
      chisqP<-c(presult[,1],presult[,3])
      names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
      chisqP<-c(presult[,3])
      names(chisqP)<-c(paste0(row.names(presult),".P"))
      result<-c(coefs,tvalues,SEvalues,chisqP)
      # result<-c(coefs,SEvalues,chisqP)
    }
  })
  colnames(divs1)<-colnames(divindex)
  
  divs1 = divs1[!grepl("Intercept",rownames(divs1)),]
  divs1 = as.data.frame(divs1)
  divs1$treat = rownames(divs1)
  divs1 = separate(divs1,col = treat,sep = '\\.',into = c('treat','type'))
  divs1 %<>%
    select(treat, type, everything()) %>%
    arrange(treat, type)
  divs1$combined_treat = i
  div.table = rbind(div.table, divs1)
  
  treat.table = subset(divindex, used.treat$treat1 == i)
  treat.table$block = treat$block[match(rownames(treat.table), rownames(treat))]
  treat.table$Layer = treat$Layer[match(rownames(treat.table), rownames(treat))]
  treat.table = treat.table[order(treat.table$block),]
  treat.table = treat.table[order(treat.table$Layer),]
  
  ctrl.table = subset(divindex, used.treat$treat1 == "C")
  ctrl.table$block = treat$block[match(rownames(ctrl.table), rownames(treat))]
  ctrl.table$Layer = treat$Layer[match(rownames(ctrl.table), rownames(treat))]
  ctrl.table = ctrl.table[order(ctrl.table$block),]
  ctrl.table = ctrl.table[order(ctrl.table$Layer),]
  
  treat.table = treat.table[,-c(2:3), drop = F]
  ctrl.table = ctrl.table[,-c(2:3), drop = F]
  
  ratio.table = treat.table/ctrl.table
  ratio.table = ratio.table-1
  ratio.result = rbind(colMeans(ratio.table), sapply(ratio.table, se))
  ratio.result = as.data.frame(ratio.result)
  ratio.result$treat = i
  ratio.result$item = c("ratio","se")
  mean.table = rbind(mean.table, ratio.result)
}
save.file(div.table,prefix = prefix.m,"all layers LMM")
save.file(mean.table,prefix = prefix.m,"all layers ratio")

## every layer ####
treat$precip = ifelse(treat$precip!=0,1,0)
div.table = data.frame(stringsAsFactors = F)
mean.table = data.frame(stringsAsFactors = F)
for (i in c("RP","EP","W")){
  for (j in unique(treat$Layer)){
    divindex <- used.data
    used.treat = subset(treat, Layer == j)
    
    sames = match.name(rn.list = list(divindex = divindex, used.treat = used.treat))
    divindex = sames$divindex
    used.treat = sames$used.treat
    
    divindex = subset(divindex, treat1%in% c(i,"C"))
    {
      if(i %in% c("RP","EP")){
        divindex$num_treat = used.treat$precip[match(divindex$treat1, used.treat$combined_treat1)]
      } else {
        divindex$num_treat = used.treat$warm[match(divindex$treat1, used.treat$combined_treat1)]
      }
    }
    divindex$layer = used.treat$layer[match(rownames(divindex), rownames(used.treat))]
    
    used.treat = divindex[,2:5]
    divindex = divindex[,-c(2:5), drop = F]
    
    used.treat$layer = (used.treat$layer-min(used.treat$layer))/(max(used.treat$layer)-min(used.treat$layer))
    # scale the alpha diversities
    # divindex<-scale(divindex)
    
    divs1<-sapply(1:ncol(divindex),function(j){
      message("Now j=",j," in ",ncol(divindex),". ",date())
      if (length(unique(divindex[,j]))<3){
        result<-rep(NA,38)
      } else {
        div<-data.frame(divtest=divindex[,j],used.treat)
        fm<-lmer(divtest~num_treat+(1|block),data=div)
        presult<-car::Anova(fm,type=2)
        coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
        names(coefs)<-paste0(names(coefs),".mean")
        SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
        names(SEvalues)<-paste0(names(SEvalues),".se")
        tvalues<-coef(summary(fm))[ , "t value"]#t values
        names(tvalues)<-paste0(names(tvalues),".t")
        chisqP<-c(presult[,1],presult[,3])
        names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
        chisqP<-c(presult[,3])
        names(chisqP)<-c(paste0(row.names(presult),".P"))
        result<-c(coefs,tvalues,SEvalues,chisqP)
        # result<-c(coefs,SEvalues,chisqP)
      }
    })
    colnames(divs1)<-colnames(divindex)
    
    divs1 = divs1[!grepl("Intercept",rownames(divs1)),]
    divs1 = as.data.frame(divs1)
    divs1$treat = rownames(divs1)
    divs1 = separate(divs1,col = treat,sep = '\\.',into = c('treat','type'))
    divs1 %<>%
      select(treat, type, everything()) %>%
      arrange(treat, type)
    divs1$combined_treat = i
    divs1$Layer = j
    div.table = rbind(div.table, divs1)
    
    treat.table = subset(divindex, used.treat$treat1 == i)
    treat.table$block = treat$block[match(rownames(treat.table), rownames(treat))]
    treat.table$Layer = treat$Layer[match(rownames(treat.table), rownames(treat))]
    treat.table = treat.table[order(treat.table$block),]
    treat.table = treat.table[order(treat.table$Layer),]
    
    ctrl.table = subset(divindex, used.treat$treat1 == "C")
    ctrl.table$block = treat$block[match(rownames(ctrl.table), rownames(treat))]
    ctrl.table$Layer = treat$Layer[match(rownames(ctrl.table), rownames(treat))]
    ctrl.table = ctrl.table[order(ctrl.table$block),]
    ctrl.table = ctrl.table[order(ctrl.table$Layer),]
    
    treat.table = treat.table[,-c(2,3),drop = F]
    ctrl.table = ctrl.table[,-c(2,3),drop = F]
    
    ratio.table = treat.table/ctrl.table
    ratio.table = ratio.table-1
    ratio.result = rbind(colMeans(ratio.table), sapply(ratio.table, se))
    ratio.result = as.data.frame(ratio.result)
    ratio.result$treat = i
    ratio.result$Layer = j
    ratio.result$item = c("ratio","se")
    mean.table = rbind(mean.table, ratio.result)
  }
}
save.file(mean.table,prefix = prefix.m,"every layers ratio")
save.file(div.table,prefix = prefix.m,"every layers LMM")

## topsoil and subsoil ####
treat$precip = ifelse(treat$precip!=0,1,0)
div.table = data.frame(stringsAsFactors = F)
mean.table = data.frame(stringsAsFactors = F)
for (i in c("RP","EP","W")){
  for (j in c("L12","L34")){
    divindex <- used.data
    {
      if (j == "L12"){
        used.treat = subset(treat, depth1%in% c(0,10))
      } else {
        used.treat = subset(treat, depth1%in% c(20,30))
      }
    }
    sames = match.name(rn.list = list(divindex = divindex, used.treat = used.treat))
    divindex = sames$divindex
    used.treat = sames$used.treat
    
    divindex = subset(divindex, treat1%in% c(i,"C"))
    {
      if(i %in% c("RP","EP")){
        divindex$num_treat = used.treat$precip[match(divindex$treat1, used.treat$combined_treat1)]
      } else {
        divindex$num_treat = used.treat$warm[match(divindex$treat1, used.treat$combined_treat1)]
      }
    }
    divindex$layer = used.treat$layer[match(rownames(divindex), rownames(used.treat))]
    used.treat = divindex[,2:5]
    used.treat$layer = used.treat$layer - min(used.treat$layer)
    
    divindex = divindex[,-c(2:5), drop = F]
    
    # divindex<-scale(divindex)
    divs1<-sapply(1:ncol(divindex),function(j){
      message("Now j=",j," in ",ncol(divindex),". ",date())
      if (length(unique(divindex[,j]))<3){
        result<-rep(NA,38)
      } else {
        div<-data.frame(divtest=divindex[,j],used.treat)
        fm<-lmer(divtest~num_treat*layer+(1|block),data=div)
        # fm<-lmer(divtest~num_treat+layer+(1|block),data=div)
        presult<-car::Anova(fm,type=2)
        coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
        names(coefs)<-paste0(names(coefs),".mean")
        SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
        names(SEvalues)<-paste0(names(SEvalues),".se")
        tvalues<-coef(summary(fm))[ , "t value"]#t values
        names(tvalues)<-paste0(names(tvalues),".t")
        chisqP<-c(presult[,1],presult[,3])
        names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
        chisqP<-c(presult[,3])
        names(chisqP)<-c(paste0(row.names(presult),".P"))
        result<-c(coefs,tvalues,SEvalues,chisqP)
        # result<-c(coefs,SEvalues,chisqP)
      }
    })
    colnames(divs1)<-colnames(divindex)
    
    divs1 = divs1[!grepl("Intercept",rownames(divs1)),]
    divs1 = as.data.frame(divs1)
    divs1$treat = rownames(divs1)
    divs1 = separate(divs1,col = treat,sep = '\\.',into = c('treat','type'))
    divs1 %<>%
      select(treat, type, everything()) %>%
      arrange(treat, type)
    divs1$combined_treat = i
    divs1$Layer = j
    div.table = rbind(div.table, divs1)
    
    treat.table = subset(divindex, used.treat$treat1 == i)
    treat.table$block = treat$block[match(rownames(treat.table), rownames(treat))]
    treat.table$Layer = treat$Layer[match(rownames(treat.table), rownames(treat))]
    treat.table = treat.table[order(treat.table$block),]
    treat.table = treat.table[order(treat.table$Layer),]
    
    ctrl.table = subset(divindex, used.treat$treat1 == "C")
    ctrl.table$block = treat$block[match(rownames(ctrl.table), rownames(treat))]
    ctrl.table$Layer = treat$Layer[match(rownames(ctrl.table), rownames(treat))]
    ctrl.table = ctrl.table[order(ctrl.table$block),]
    ctrl.table = ctrl.table[order(ctrl.table$Layer),]
    
    treat.table = treat.table[,-c(2:3), drop = F]
    ctrl.table = ctrl.table[,-c(2:3), drop = F]
    
    ratio.table = treat.table/ctrl.table
    ratio.table = ratio.table-1
    ratio.result = rbind(colMeans(ratio.table), sapply(ratio.table, se))
    ratio.result = as.data.frame(ratio.result)
    ratio.result$treat = i
    ratio.result$Layer = j
    ratio.result$item = c("ratio","se")
    mean.table = rbind(mean.table, ratio.result)
  }
}
save.file(mean.table,prefix = prefix.m,"top and subsoil ratio")
save.file(div.table,prefix = prefix.m,"top and subsoil LMM")

# cor depth ####
library(dplyr)
library(Rmisc)
library(reshape2)
library(ieggr)
library(ggplot2)
library(ggpubr)
if (!("smatr" %in% installed.packages()[, "Package"])) {
  install.packages("smatr")
}
library(smatr)

source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\rsquaredglmm.r")
source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\tdcm.mean.no.compare.R")

result.data = sample.aver.number
result.data$year1 = treat$depth1
result.data$treat = treat$combined_treat1
result.data$sample = rownames(result.data)
result.data$block = treat$block

# calculate the R2 and p value
model_results <- lapply(names(result.data)[-c((ncol(result.data)-3):ncol(result.data))], function(colname) {
  result.data %>%
    do(tdcm.mean.no.compare(.[, colname, drop = F], ., rand = 1000, 
                            scale.num = F, alpha = T)) %>%
  cbind(colname) %>%
    cbind(rownames(.))
}) %>% bind_rows()
colnames(model_results)[7:8] = c("index","item")
write.csv(model_results,"copy.depth.regression.csv")

# P value
w.S1.RP = subset(model_results, item == "slope.fix")$RP
w.S1.EP = subset(model_results, item == "slope.fix")$EP
w.S1.W = subset(model_results, item == "slope.fix")$W
w.S1.ctrl = subset(model_results, item == "slope.fix")$C

Slope_test_p.table = data.frame(stringsAsFactors = F)
scale.num = F
for (i in 1:(ncol(result.data)-4)){
  column = subset(result.data, treat == "C")[,i]
  cor.depth = subset(result.data, treat == "C")$year1
  if (scale.num){
    column = scale(column)
    cor.depth = scale(cor.depth)
  }
  a = sma(column ~ cor.depth, slope.test = w.S1.RP[i])[["groupsummary"]][["Slope_test_p"]]
  b = sma(column ~ cor.depth, slope.test = w.S1.EP[i])[["groupsummary"]][["Slope_test_p"]]
  c = sma(column ~ cor.depth, slope.test = w.S1.W[i])[["groupsummary"]][["Slope_test_p"]]
  
  column1 = subset(result.data, treat == "RP")[,i]
  column2 = subset(result.data, treat == "EP")[,i]
  column3 = subset(result.data, treat == "W")[,i]
  cor.depth1 = subset(result.data, treat == "RP")$year1
  cor.depth2 = subset(result.data, treat == "EP")$year1
  cor.depth3 = subset(result.data, treat == "W")$year1
  
  if (scale.num){
    column1 = scale(column1)
    column2 = scale(column2)
    column3 = scale(column3)
    cor.depth1 = scale(cor.depth1)
    cor.depth2 = scale(cor.depth2)
    cor.depth3 = scale(cor.depth3)
  }
  d = sma(column1 ~ cor.depth1, slope.test = w.S1.ctrl[i])[["groupsummary"]][["Slope_test_p"]]
  e = sma(column2 ~ cor.depth2, slope.test = w.S1.ctrl[i])[["groupsummary"]][["Slope_test_p"]]
  f = sma(column3 ~ cor.depth3, slope.test = w.S1.ctrl[i])[["groupsummary"]][["Slope_test_p"]]
  Slope_test_p.table = rbind(Slope_test_p.table,c(a,b,c,d,e,f))
}
names(Slope_test_p.table) = c("RP1","EP1","W1","RP2","EP2","W2")
Slope_test_p.table$phylum = names(result.data)[1:(ncol(result.data)-4)]
# write.csv(Slope_test_p.table,paste0(prefix.m,".","taxa.P.richness.depth.regression.csv"))
write.csv(Slope_test_p.table,"group.P.copy.depth.regression.csv")

# linear plot ####
library(dplyr)
library(reshape2)
source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\rsquaredglmm.r")
source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\tdcm.mean.no.compare.R")
prefix = prefix.m
result.data = sample.aver.number
colnames(result.data)[1] = "copy"
treat = read.csv(treat.file, header = T, row.names = 1, sep = ",")
treat = subset(treat, plant.type == "TS")
samps = match.name(rn.list = list(result.data = result.data, treat = treat))
dim(result.data)
result.data = samps$result.data
dim(result.data)
treat = samps$treat

result.data$year1 = treat$depth1
result.data$treat = treat$combined_treat1
result.data$sample = rownames(result.data)
result.data$block = treat$block

# calculate the R2 and p value
model_results <- result.data %>%
  do(tdcm.mean.no.compare(.[, "copy", drop = F], ., rand = 1000,
                          scale.num = F, alpha = T))

model_results$item <- rownames(model_results)
# write.csv(model_results,paste0(prefix,".","total.richness.depth.regression.csv"))
p.results <- model_results[model_results$item == c("P.typeII"), ]
R.results <- model_results[model_results$item == c("R2M"), ]
w.results = model_results[model_results$item == c("slope.fix"), ]
se.results = model_results[model_results$item == c("slope.se"), ]

p.results <- melt(p.results, id = c("item"))
R.results <- melt(R.results, id = c("item"))
w.results <- melt(w.results, id = c("item"))
se.results <- melt(se.results, id = c("item"))

# if there are duplicate column names, mutate function will not work.
colnames(p.results)[2:3] <- c("treat", "P1")
colnames(R.results)[2:3] <- c("treat", "R1")
colnames(w.results)[2:3] <- c("treat", "w1")
colnames(se.results)[2:3] <- c("treat", "se1")
alpha <- 0.05
# 添加显著性列，根据 p 值判断显著性
p.results <- p.results %>%
  mutate(significant1 = ifelse(P1 < alpha, "significant", "not significant"))
# 将模型结果与原始数据合并
data_with_results <- left_join(result.data, p.results, by = c("treat"))
data_with_results <- left_join(data_with_results, R.results, by = c("treat"))
data_with_results <- left_join(data_with_results, w.results, by = c("treat"))
data_with_results <- left_join(data_with_results, se.results, by = c("treat"))

R.S1.RP <- round(unique(subset(data_with_results, treat == "RP")$R1), 3)
R.S1.EP <- round(unique(subset(data_with_results, treat == "EP")$R1), 3)
R.S1.W <- round(unique(subset(data_with_results, treat == "W")$R1), 3)
R.S1.ctrl <- round(unique(subset(data_with_results, treat == "C")$R1), 3)

P.S1.RP <- round(unique(subset(data_with_results, treat == "RP")$P1), 3)
P.S1.EP <- round(unique(subset(data_with_results, treat == "EP")$P1), 3)
P.S1.W <- round(unique(subset(data_with_results, treat == "W")$P1), 3)
P.S1.ctrl <- round(unique(subset(data_with_results, treat == "C")$P1), 3)

w.S1.RP <- round(unique(subset(data_with_results, treat == "RP")$w1), 3)
w.S1.EP <- round(unique(subset(data_with_results, treat == "EP")$w1), 3)
w.S1.W <- round(unique(subset(data_with_results, treat == "W")$w1), 3)
w.S1.ctrl <- round(unique(subset(data_with_results, treat == "C")$w1), 3)

se.S1.RP <- round(unique(subset(data_with_results, treat == "RP")$se1), 3)
se.S1.EP <- round(unique(subset(data_with_results, treat == "EP")$se1), 3)
se.S1.W <- round(unique(subset(data_with_results, treat == "W")$se1), 3)
se.S1.ctrl <- round(unique(subset(data_with_results, treat == "C")$se1), 3)
library(ggplot2)
p1 =  ggplot(subset(data_with_results, treat %in% c("EP","W","RP","C")), aes(x = year1, y = copy,group = treat, color = treat)) +
  geom_point(size = 5, alpha = 0.4,show.legend = F) +
  scale_color_manual(values = c("#1874CD", "#CD2626","#DAA520","#808080"),  limits = c("EP","W","RP","C")) +
  geom_line(
    data = data_with_results %>% filter(treat == "RP"),
    aes(linetype = significant1), show.legend = F,
    stat="smooth",method = "lm", se = F, linewidth = 2,alpha = 1
  ) +
  geom_line(
    data = data_with_results %>% filter(treat == "EP"),
    aes(linetype = significant1), show.legend = F,
    stat="smooth",method = "lm", se = F, linewidth = 2,alpha = 1
  ) +
  geom_line(
    data = data_with_results %>% filter(treat == "W"),
    aes(linetype = significant1), show.legend = F,
    stat="smooth",method = "lm", se = F, linewidth = 2,alpha = 1
  ) +
  geom_line(
    data = data_with_results %>% filter(treat == "C"),
    aes(linetype = significant1), show.legend = F,
    stat="smooth",method = "lm", se = F, linewidth = 2,alpha = 1
  ) +
  scale_linetype_manual(values = c("significant" = "solid", "not significant" = "dashed")) +
  xlab("Soil Depth")+
  ylab("rrn Copy Number")+
  annotate("text",
           x = 20,
           y = 1.78,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S1.RP, digits = 3, format = "f")) %+-% .(formatC(se.S1.RP, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S1.RP, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S1.RP, digits = 3, format = "f")))),
           colour = "#DAA520"
  ) +
  annotate("text",
           x = 20,
           y = 1.74,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S1.EP, digits = 3, format = "f")) %+-% .(formatC(se.S1.EP, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S1.EP, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S1.EP, digits = 3, format = "f")))),
           colour = "#1874CD"
  ) +
  annotate("text",
           x = 20,
           y = 1.70,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S1.W, digits = 3, format = "f")) %+-% .(formatC(se.S1.W, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S1.W, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S1.W, digits = 3, format = "f")))),
           colour = "#CD2626"
  ) +
  annotate("text",
           x = 20,
           y = 1.66,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S1.ctrl, digits = 3, format = "f")) %+-% .(formatC(se.S1.ctrl, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S1.ctrl, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S1.ctrl, digits = 3, format = "f")))),
           colour = "#808080"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(
      size = 17,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 15,
      face = "bold",
      angle = 0,
      margin = margin(t = 0)
    )
  )
p1
ggsave("rrn copy number.depth.TS.pdf", width = 6.31, height = 5.71, units = "in")

# barplot ####
library(ggbeeswarm)
library(ggplot2)
result.data = subset(result.data, treat %in% c("RP","EP","W","C"))
result.data$group = paste0(result.data$year1,".",result.data$treat)
result.data$year1 = as.character(result.data$year1)
p1 <- ggplot(data=result.data,mapping=aes(x=group,y=copy,
                                         color=treat,fill=treat))+
  stat_summary(fun = mean, geom="bar",alpha=0.0,
               # position = "dodge",
               width = 0.54,
               show.legend=F)+
  scale_color_manual(values = c("#1874CD", "#CD2626","#DAA520","#808080"),  limits = c("EP","W","RP","C"))
p1

p2 <- p1+stat_summary(fun.data=mean_sdl, fun.args = list(mult=1/2),
                      geom="errorbar", width=0.3,show.legend=F)
p2

p3 <- p2+
  geom_jitter(shape=16,alpha=0.8,size=2,width = 0.15,show.legend=F)
p3 +
  xlab("")+
  ylab("rnn Copy Number")+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(
      size = 17,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 15,
      face = "bold",
      angle = 0,
      margin = margin(t = 0)
    )
  )
ggsave("rrn copy number barplot.pdf", width = 9.98, height = 5.71, units = "in")

# cor with env ####
treat = subset(treat, plant.type == "TS"&combined_treat1 %in% c("RP","EP","W","C"))
sampc=match.name(rn.list = list(treat = treat,env.data = env.data,sample.aver.number = sample.aver.number))
env.data = sampc$env.data
treat = sampc$treat
sample.aver.number = sampc$sample.aver.number
cor1 = cor.test(env.data$BNPP,sample.aver.number$sample.aver.number,method = "pearson")
cor2 = cor.test(env.data$NO3,sample.aver.number$sample.aver.number,method = "pearson")
cor3 = cor.test(env.data$NH4,sample.aver.number$sample.aver.number,method = "pearson")
cor4 = cor.test(env.data$moisture.one.nighbor.aver,sample.aver.number$sample.aver.number,method = "pearson")
result.df = data.frame(BNPP = c(cor1[["estimate"]][["cor"]],cor1[["p.value"]]),
                          NO3 = c(cor2[["estimate"]][["cor"]],cor2[["p.value"]]),
                          NH4 = c(cor3[["estimate"]][["cor"]],cor3[["p.value"]]),
                          Moisture = c(cor4[["estimate"]][["cor"]],cor4[["p.value"]]))
write.csv(result.df, "unweighted.rrn.env.csv")

