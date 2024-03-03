setwd("C:/Users/True/OneDrive/桌面")
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

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"
index.file = file.choose()
# index.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\03242023 Inner Mongolia TS\\results\\network\\Fungi\\BTS.index.spearman.Fungi.0.5.csv"
prefix.m = "Layer Bacteria 0.7"

result.data = read.csv(index.file, header = T, row.names = 1, sep = ",")
result.data = result.data[,c(1:4,9:12,17,22:26)]
treat = read.csv(treat.file, header = T, row.names = 1, sep = ",")
samps = match.name(rn.list = list(result.data = result.data, treat = treat))
dim(result.data)
result.data = samps$result.data
dim(result.data)
treat = samps$treat

result.data$year1 = treat$depth1
result.data$treat = treat$combined_treat1
result.data$sample = rownames(result.data)
result.data$block = treat$block

model_results <- lapply(names(result.data)[-c((ncol(result.data)-3):ncol(result.data))], function(colname) {
  result.data %>%
    do(tdcm.mean.no.compare(.[, colname, drop = F], ., rand = 1000,
                            scale.num = T, alpha = T)) %>%
    cbind(colname) %>%
    cbind(rownames(.))
}) %>% bind_rows()
colnames(model_results)[5:6] = c("index","item")
write.csv(model_results,paste0(prefix.m,".","network.index.depth.regression.csv"))

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
write.csv(Slope_test_p.table,paste0(prefix.m,".","group.P.network.index.depth.regression.csv"))
