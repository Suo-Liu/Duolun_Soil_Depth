# 不同层的物种组成
# sample的平均加权长度
# body size仅能分到50%的物种，也许不能代表群落整体；或者只在figure legend里面说一下，而且表明他们都是abundant的物种
# abundance weighted body size ####
setwd("C:/Users/True/OneDrive/桌面")

library(ieggr)
library(picante)

# protist
com.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\protist.zotu_resample_2961.txt"
tree.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_rooted_tree.nwk"
clas.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_tax.txt"
prefix.m = "Protist"

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"

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
dim(clas)
clas = clas[!is.na(clas$body.size),]
dim(clas)

treat = read.csv(treat.file, row.names = 1, header = T,sep = ",")
treat = subset(treat, plant.type == "TS"& combined_treat1 %in% c("C","RP","EP","W"))
sampc=match.name(rn.list = list(treat = treat,comm = comm))
dim(comm)
comm = sampc$comm
treat = sampc$treat
comm = comm[,colSums(comm)>0]
dim(comm)

spc=match.name(rn.list=list(clas=clas), cn.list = list(comm = comm))
comm=spc$comm
clas=spc$clas

# 是否需要重新resample
# 问天娇 rrn copy number
# 不同处理和叠加处理对这一关系的影响
comm.wt = comm
comm.wt = ifelse(comm>0,1,0)
rowSums(comm.wt)


comm.wt = comm.wt*clas$body.size
avg.size = rowSums(comm.wt)/rowSums(comm)
avg.size = as.data.frame(avg.size)
avg.size$year1 = treat$depth1[match(rownames(avg.size), rownames(treat))]
avg.size$Warm = treat$Warm[match(rownames(avg.size), rownames(treat))]
avg.size$Precip = treat$Precip[match(rownames(avg.size), rownames(treat))]
avg.size$plant.type = treat$plant.type[match(rownames(avg.size), rownames(treat))]

avg.size$treat = treat$combined_treat1[match(rownames(avg.size), rownames(treat))]

avg.size = avg.size[order(avg.size$year1),]
avg.size = avg.size[order(avg.size$treat),]
avg.size$block = treat$block[match(rownames(avg.size), rownames(treat))]

setwd("C:\\Users\\True\\OneDrive\\桌面")
library(dplyr)
library(magrittr)
library(reshape2)
library(ggplot2)
library(ggpubr)
if (!("smatr" %in% installed.packages()[, "Package"])) {
  install.packages("smatr")
}
library(smatr)
source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\rsquaredglmm.r")
source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\tdcm.mean.no.compare.R")

prefix = "Protist.unweight"
model_results <-  avg.size %>%
    tdcm.mean.no.compare(.[, 1, drop = F], ., rand = 1000, 
                            scale.num = F,alpha = T)
write.csv(model_results,paste0(prefix,".","body.size.depth.regression.csv"))

model_results$item = rownames(model_results)
# P value
w.S1.RP = subset(model_results, item == "slope.fix")$RP
w.S1.EP = subset(model_results, item == "slope.fix")$EP
w.S1.W = subset(model_results, item == "slope.fix")$W
w.S1.ctrl = subset(model_results, item == "slope.fix")$C

Slope_test_p.table = data.frame(stringsAsFactors = F)
scale.num = F
for (i in 1:(ncol(avg.size)-6)){
  column = subset(avg.size, treat == "C")[,i]
  cor.depth = subset(avg.size, treat == "C")$year1
  if (scale.num){
    column = scale(column)
    cor.depth = scale(cor.depth)
  }
  a = sma(column ~ cor.depth, slope.test = w.S1.RP[i])[["groupsummary"]][["Slope_test_p"]]
  b = sma(column ~ cor.depth, slope.test = w.S1.EP[i])[["groupsummary"]][["Slope_test_p"]]
  c = sma(column ~ cor.depth, slope.test = w.S1.W[i])[["groupsummary"]][["Slope_test_p"]]
  
  column1 = subset(avg.size, treat == "RP")[,i]
  column2 = subset(avg.size, treat == "EP")[,i]
  column3 = subset(avg.size, treat == "W")[,i]
  cor.depth1 = subset(avg.size, treat == "RP")$year1
  cor.depth2 = subset(avg.size, treat == "EP")$year1
  cor.depth3 = subset(avg.size, treat == "W")$year1
  
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
Slope_test_p.table$phylum = names(avg.size)[1:(ncol(avg.size)-6)]
write.csv(Slope_test_p.table,paste0(prefix,".","group.P.body.size.depth.regression.csv"))

# LMM ####
setwd("C:/Users/True/OneDrive/桌面")
library(ieggr)
library(lme4)
library(car)
library(tidyr)
library(dplyr)
library(magrittr)
se <- function(x){
  sd(x) / sqrt(length(x))
}
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"
treat<-read.csv(treat.file,header = T,row.names = 1)
treat = subset(treat, plant.type == "TS")

used.data<- avg.size[,1,drop = F]

sames = match.name(rn.list = list(used.data = used.data, treat = treat))
used.data = sames$used.data
treat = sames$treat

used.data$treat1 = treat$combined_treat1[match(rownames(used.data), rownames(treat))]
used.data$block = treat$block[match(rownames(used.data), rownames(treat))]
### All layers ####
treat$precip = ifelse(treat$precip!=0,1,0)
div.table = data.frame(stringsAsFactors = F)
mean.table = data.frame(stringsAsFactors = F)
for (i in c("RP","EP","W")){
  divindex <- used.data
  divindex = subset(divindex, treat1%in% c(i,"C"))
  {
    if(i %in% c("RP","EP")){
      divindex$num_treat = treat$precip[match(divindex$treat1, treat$combined_treat1)]
    } else {
      divindex$num_treat = treat$warm[match(divindex$treat1, treat$combined_treat1)]
    }
  }
  divindex$layer = treat$layer[match(rownames(divindex), rownames(treat))]
  
  used.treat = divindex[,2:5,drop=F]
  divindex = divindex[,-c(2:5),drop=F]
  
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
  
  treat.table = treat.table[,-c(2:3),drop = F]
  ctrl.table = ctrl.table[,-c(2:3),drop = F]
  
  ratio.table = treat.table/ctrl.table
  ratio.table = ratio.table-1
  ratio.result = rbind(colMeans(ratio.table), sapply(ratio.table, se))
  ratio.result = as.data.frame(ratio.result)
  ratio.result$treat = i
  ratio.result$item = c("ratio","se")
  mean.table = rbind(mean.table, ratio.result)
}
save.file(div.table,prefix = prefix.m,"all layers LMM")
save.file(mean.table,prefix = prefix.m,"all layers mean")

### every layer ####
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
    
    used.treat = divindex[,2:5,drop=F]
    divindex = divindex[,-c(2:5),drop=F]
    
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
    
    treat.table = treat.table[,-c(2:3),drop = F]
    ctrl.table = ctrl.table[,-c(2:3),drop = F]
    
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
save.file(mean.table,prefix = prefix.m,"every layers mean")
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
    used.treat = divindex[,2:5,drop=F]
    used.treat$layer = used.treat$layer - min(used.treat$layer)
    
    divindex = divindex[,-c(2:5),drop=F]
    
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
    
    treat.table = treat.table[,-c(2:3),drop = F]
    ctrl.table = ctrl.table[,-c(2:3),drop = F]
    
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
save.file(mean.table,prefix = prefix.m,"top and subsoil mean")
save.file(div.table,prefix = prefix.m,"top and subsoil LMM")

boxplot(avg.size$avg.size ~ avg.size$depth)

t.test(subset(avg.size, depth == 10)$avg.size,subset(avg.size, depth == 50)$avg.size )
t.test(subset(avg.size, depth == 20)$avg.size,subset(avg.size, depth == 50)$avg.size )
t.test(subset(avg.size, depth == 30)$avg.size,subset(avg.size, depth == 50)$avg.size )

# heatmap,表示沿不同taxa沿深度的丰度分布，并给它们标上body size
# 筛选在处理下丰度变化显著的ASV，并生成comm tree clas
setwd("C:/Users/True/OneDrive/桌面")
library(ieggr)
com.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\protist.zotu_resample_2961.txt"
tree.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_rooted_tree.nwk"
clas.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_tax.txt"
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"
# C ####
comm <- t(read.table(com.file,
                     header = TRUE, sep = "\t", row.names = 1,
                     as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                     check.names = FALSE
))
name.row = rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] = "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] = "YD3L3"
rownames(comm) = name.row

clas <- read.table(clas.file,
                   header = TRUE, sep = "\t", row.names = 1,
                   as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                   check.names = FALSE
)
tree = lazyopen(tree.file)
treat <- read.csv(treat.file, header = T, row.names = 1)
treat = subset(treat, plant.type == "TS"&combined_treat1 %in% c("C"))

library(ieggr)
samps = match.name(rn.list = list(treat = treat, comm = comm))
dim(comm)
comm = samps$comm
comm <- comm[, colSums(comm) > 0] # sometimes, removing samples will make some OTUs have no read across remained samples.
dim(comm)
treat = samps$treat
strees = match.name(cn.list = list(comm = comm), tree.list = list(tree = tree))
dim(comm)
comm = samps$comm
comm <- comm[, colSums(comm) > 0] # sometimes, removing samples will make some OTUs have no read across remained samples.
dim(comm)
tree = strees$tree

# leave the rare species
comm = comm[, colSums(comm>0) >= 6]
dim(comm)
divindex = comm
divindex <- divindex[match(row.names(treat), row.names(divindex)), ]
# scale the alpha diversities
divindex <- scale(divindex)
library(lme4)
library(car)
divs2 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat)
  fm <- lmer(divtest ~ depth + (1 | block), data = div)
  presult<-car::Anova(fm,type=2)
  coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
  names(coefs)<-paste0(names(coefs),".mean")
  SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
  names(SEvalues)<-paste0(names(SEvalues),".se")
  # tvalues<-coef(summary(fm))[ , "t value"]#t values
  # names(tvalues)<-paste0(names(tvalues),".t")
  # chisqP<-c(presult[,1],presult[,3])
  # names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
  chisqP<-c(presult[,3])
  names(chisqP)<-c(paste0(row.names(presult),".P"))
  # result<-c(coefs,tvalues,SEvalues,chisqP)
  result<-c(coefs,SEvalues,chisqP)
})
colnames(divs2) <- colnames(divindex)
divs2 = divs2[!grepl("Intercept",rownames(divs2)),]
divs2 = as.data.frame(divs2)
divs2$treat = rownames(divs2)
library(tidyr)
divs2 = separate(divs2,col = treat,sep = '\\.',into = c('treat','type'))
library(dplyr)
library(magrittr)
divs2 %<>%
  select(treat, type, everything()) %>%
  arrange(treat, type)

table(divs2[1,3:ncol(divs2)]<0.05)

sig.name = colnames(divs2)[which(divs2[1,3:ncol(divs2)]<0.05)+2]

dim(divs2)
sel.div = divs2[,colnames(divs2)%in%sig.name]
sel.div[2,] = ifelse(sel.div[2,]>0,1,0)
dim(sel.div)

dim(comm)
comm = comm[,colnames(comm)%in%sig.name]
dim(comm)

comm = as.data.frame(comm)
comm$depth = treat$depth[match(rownames(comm),rownames(treat))]

test.comm = scale(comm)

test.comm = test.comm[match(rownames(treat), rownames(test.comm)),]
test.comm = test.comm[order(treat$depth),]

library(pheatmap)
pheatmap(t(test.comm), cluster_cols = F)

# RP ####
comm <- t(read.table(com.file,
                     header = TRUE, sep = "\t", row.names = 1,
                     as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                     check.names = FALSE
))
name.row = rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] = "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] = "YD3L3"
rownames(comm) = name.row

clas <- read.table(clas.file,
                   header = TRUE, sep = "\t", row.names = 1,
                   as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                   check.names = FALSE
)
tree = lazyopen(tree.file)
treat <- read.csv(treat.file, header = T, row.names = 1)
treat = subset(treat, plant.type == "TS"&combined_treat1 %in% c("RP"))

library(ieggr)
samps = match.name(rn.list = list(treat = treat, comm = comm))
dim(comm)
comm = samps$comm
comm <- comm[, colSums(comm) > 0] # sometimes, removing samples will make some OTUs have no read across remained samples.
dim(comm)
treat = samps$treat
strees = match.name(cn.list = list(comm = comm), tree.list = list(tree = tree))
dim(comm)
comm = samps$comm
comm <- comm[, colSums(comm) > 0] # sometimes, removing samples will make some OTUs have no read across remained samples.
dim(comm)
tree = strees$tree

# leave the rare species
comm = comm[, colSums(comm>0) >= 6]
dim(comm)
divindex = comm
divindex <- divindex[match(row.names(treat), row.names(divindex)), ]
# scale the alpha diversities
divindex <- scale(divindex)
library(lme4)
library(car)
divs2 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat)
  fm <- lmer(divtest ~ depth + (1 | block), data = div)
  presult<-car::Anova(fm,type=2)
  coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
  names(coefs)<-paste0(names(coefs),".mean")
  SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
  names(SEvalues)<-paste0(names(SEvalues),".se")
  # tvalues<-coef(summary(fm))[ , "t value"]#t values
  # names(tvalues)<-paste0(names(tvalues),".t")
  # chisqP<-c(presult[,1],presult[,3])
  # names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
  chisqP<-c(presult[,3])
  names(chisqP)<-c(paste0(row.names(presult),".P"))
  # result<-c(coefs,tvalues,SEvalues,chisqP)
  result<-c(coefs,SEvalues,chisqP)
})
colnames(divs2) <- colnames(divindex)
divs2 = divs2[!grepl("Intercept",rownames(divs2)),]
divs2 = as.data.frame(divs2)
divs2$treat = rownames(divs2)
library(tidyr)
divs2 = separate(divs2,col = treat,sep = '\\.',into = c('treat','type'))
library(dplyr)
library(magrittr)
divs2 %<>%
  select(treat, type, everything()) %>%
  arrange(treat, type)

table(divs2[1,3:ncol(divs2)]<0.05)

sig.name = colnames(divs2)[which(divs2[1,3:ncol(divs2)]<0.05)+2]

dim(divs2)
sel.div = divs2[,colnames(divs2)%in%sig.name]
sel.div[2,] = ifelse(sel.div[2,]>0,1,0)
dim(sel.div)

dim(comm)
comm = comm[,colnames(comm)%in%sig.name]
dim(comm)

comm = as.data.frame(comm)
comm$depth = treat$depth[match(rownames(comm),rownames(treat))]

test.comm = scale(comm)

test.comm = test.comm[match(rownames(treat), rownames(test.comm)),]
test.comm = test.comm[order(treat$depth),]

library(pheatmap)
pheatmap(t(test.comm), cluster_cols = F)

# W ####
comm <- t(read.table(com.file,
                     header = TRUE, sep = "\t", row.names = 1,
                     as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                     check.names = FALSE
))
name.row = rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] = "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] = "YD3L3"
rownames(comm) = name.row

clas <- read.table(clas.file,
                   header = TRUE, sep = "\t", row.names = 1,
                   as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                   check.names = FALSE
)
tree = lazyopen(tree.file)
treat <- read.csv(treat.file, header = T, row.names = 1)
treat = subset(treat, plant.type == "TS"&combined_treat1 %in% c("W"))

library(ieggr)
samps = match.name(rn.list = list(treat = treat, comm = comm))
dim(comm)
comm = samps$comm
comm <- comm[, colSums(comm) > 0] # sometimes, removing samples will make some OTUs have no read across remained samples.
dim(comm)
treat = samps$treat
strees = match.name(cn.list = list(comm = comm), tree.list = list(tree = tree))
dim(comm)
comm = samps$comm
comm <- comm[, colSums(comm) > 0] # sometimes, removing samples will make some OTUs have no read across remained samples.
dim(comm)
tree = strees$tree

# leave the rare species
comm = comm[, colSums(comm>0) >= 6]
dim(comm)
divindex = comm
divindex <- divindex[match(row.names(treat), row.names(divindex)), ]
# scale the alpha diversities
divindex <- scale(divindex)
library(lme4)
library(car)
divs2 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat)
  fm <- lmer(divtest ~ depth + (1 | block), data = div)
  presult<-car::Anova(fm,type=2)
  coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
  names(coefs)<-paste0(names(coefs),".mean")
  SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
  names(SEvalues)<-paste0(names(SEvalues),".se")
  # tvalues<-coef(summary(fm))[ , "t value"]#t values
  # names(tvalues)<-paste0(names(tvalues),".t")
  # chisqP<-c(presult[,1],presult[,3])
  # names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
  chisqP<-c(presult[,3])
  names(chisqP)<-c(paste0(row.names(presult),".P"))
  # result<-c(coefs,tvalues,SEvalues,chisqP)
  result<-c(coefs,SEvalues,chisqP)
})
colnames(divs2) <- colnames(divindex)
divs2 = divs2[!grepl("Intercept",rownames(divs2)),]
divs2 = as.data.frame(divs2)
divs2$treat = rownames(divs2)
library(tidyr)
divs2 = separate(divs2,col = treat,sep = '\\.',into = c('treat','type'))
library(dplyr)
library(magrittr)
divs2 %<>%
  select(treat, type, everything()) %>%
  arrange(treat, type)

table(divs2[1,3:ncol(divs2)]<0.05)

sig.name = colnames(divs2)[which(divs2[1,3:ncol(divs2)]<0.05)+2]

dim(divs2)
sel.div = divs2[,colnames(divs2)%in%sig.name]
sel.div[2,] = ifelse(sel.div[2,]>0,1,0)
dim(sel.div)

dim(comm)
comm = comm[,colnames(comm)%in%sig.name]
dim(comm)

comm = as.data.frame(comm)

test.comm = scale(comm)

test.comm = test.comm[match(rownames(treat), rownames(test.comm)),]
test.comm = test.comm[order(treat$depth),]

library(pheatmap)
pheatmap(t(test.comm), cluster_cols = F)

# body size test ####
# 由于LMM按照effect size选择的decrease物种太少，而不显著响应的物种又太多，因此
# 需要按照热图选取响应量较大的物种；但是只有在control条件下，选择scale之后深层最大值至少为2的物种与其他有显著响应但效应量不大的物种
# 相比才有body size显著较小的结果，其他均无结果，且该结果也不好解释。
## divide group by the heatmap ####
body.data1 = colnames(test.comm)
body.data1 = as.data.frame(body.data1)
rownames(body.data1) = body.data1[,1]
body.data1$body.size = clas$body.size[match(colnames(test.comm), rownames(clas))]
body.data1 = body.data1[,-1, drop = F]

max.value = apply(test.comm[(nrow(test.comm)-2):nrow(test.comm),],2,max)
body.data1$value.in.depth = max.value[match(rownames(body.data1), names(max.value))]
body.data1$rank = ifelse(body.data1$value.in.depth>2,1,0)
table(body.data1$rank)
body.data1 = body.data1[!is.na(body.data1$body.size),]
lm(body.data1$body.size ~ body.data1$rank)
anova(lm(body.data1$body.size ~ body.data1$rank))


## divide group by LMM ####
in.sig.name = colnames(divs2)[which(divs2[1,3:ncol(divs2)]<0.05&divs2[2,3:ncol(divs2)]>0)+2]
de.sig.name = colnames(divs2)[which(divs2[1,3:ncol(divs2)]<0.05&divs2[2,3:ncol(divs2)]<0)+2]
even.sig.name = colnames(divs2)[which(divs2[1,3:ncol(divs2)]>0.05)+2]

in.sig.name = as.data.frame(in.sig.name)
de.sig.name =  as.data.frame(de.sig.name)
even.sig.name =  as.data.frame(even.sig.name)

rownames(in.sig.name) = in.sig.name[,1]
rownames(de.sig.name) = de.sig.name[,1]
rownames(even.sig.name) = even.sig.name[,1]

in.sig.name$body.size = clas$body.size[match(rownames(in.sig.name), rownames(clas))]
de.sig.name$body.size = clas$body.size[match(rownames(de.sig.name), rownames(clas))]
even.sig.name$body.size = clas$body.size[match(rownames(even.sig.name), rownames(clas))]

in.sig.name = in.sig.name[!is.na(in.sig.name$body.size),]
de.sig.name = de.sig.name[!is.na(de.sig.name$body.size),]
even.sig.name = even.sig.name[!is.na(even.sig.name$body.size),]

in.sig.name = in.sig.name[,2,drop = F]
de.sig.name = de.sig.name[,2,drop = F]
even.sig.name = even.sig.name[,2,drop = F]

in.sig.name$rank = 1
de.sig.name$rank = -1
even.sig.name$rank = 0

t.test(in.sig.name$body.size, de.sig.name$body.size)
t.test(in.sig.name$body.size, even.sig.name$body.size)
t.test(even.sig.name$body.size, de.sig.name$body.size)

mer.data = rbind(in.sig.name,de.sig.name,even.sig.name)

lm(mer.data$body.size ~ mer.data$rank)
anova(lm(mer.data$body.size ~ mer.data$rank))

