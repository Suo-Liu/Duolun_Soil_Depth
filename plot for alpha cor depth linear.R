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
## richness ####
### Bacteria ####
prefix.m = "Bacteria"
divindex.file = paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03242023 Inner Mongolia TS\\results\\11142023\\alpha diversity\\",prefix.m,"\\",prefix.m,".alpha.div.csv")

prefix = prefix.m
result.data = read.csv(divindex.file, header = T, row.names = 1, sep = ",")
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
  do(tdcm.mean.no.compare(.[, "richness", drop = F], ., rand = 1000,
                          scale.num = F))

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

p1 =  ggplot(subset(data_with_results, treat %in% c("EP","W","RP","C")), aes(x = year1, y = richness,group = treat, color = treat)) +
  geom_point(size = 5, alpha = 0.4,show.legend = F) +
  # geom_errorbar(aes(ymin = richness - se, ymax = richness + se),
  #               width = 5.5, alpha = 0.3,
  #               position = position_dodge(.0), lwd = 1.2
  # ) +
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
  xlab("Soil depth")+
  ylab("Richness")+
  annotate("text",
           x = 20,
           y = 10600,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S1.RP, digits = 3, format = "f")) %+-% .(formatC(se.S1.RP, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S1.RP, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S1.RP, digits = 3, format = "f")))),
           colour = "#DAA520"
  ) +
  annotate("text",
           x = 20,
           y = 10200,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S1.EP, digits = 3, format = "f")) %+-% .(formatC(se.S1.EP, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S1.EP, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S1.EP, digits = 3, format = "f")))),
           colour = "#1874CD"
  ) +
  annotate("text",
           x = 20,
           y = 9800,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S1.W, digits = 3, format = "f")) %+-% .(formatC(se.S1.W, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S1.W, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S1.W, digits = 3, format = "f")))),
           colour = "#CD2626"
  ) +
  annotate("text",
           x = 20,
           y = 9400,
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
# 4.76 4.93
### Fungi ####
prefix.m = "Fungi"
divindex.file = paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03242023 Inner Mongolia TS\\results\\11142023\\alpha diversity\\",prefix.m,"\\",prefix.m,".alpha.div.csv")

prefix = prefix.m
result.data = read.csv(divindex.file, header = T, row.names = 1, sep = ",")
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
  do(tdcm.mean.no.compare(.[, "richness", drop = F], ., rand = 1000, 
                          scale.num = F))

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

R.S2.RP <- round(unique(subset(data_with_results, treat == "RP")$R1), 3)
R.S2.EP <- round(unique(subset(data_with_results, treat == "EP")$R1), 3)
R.S2.W <- round(unique(subset(data_with_results, treat == "W")$R1), 3)
R.S2.ctrl <- round(unique(subset(data_with_results, treat == "C")$R1), 3)

P.S2.RP <- round(unique(subset(data_with_results, treat == "RP")$P1), 3)
P.S2.EP <- round(unique(subset(data_with_results, treat == "EP")$P1), 3)
P.S2.W <- round(unique(subset(data_with_results, treat == "W")$P1), 3)
P.S2.ctrl <- round(unique(subset(data_with_results, treat == "C")$P1), 3)

w.S2.RP <- round(unique(subset(data_with_results, treat == "RP")$w1), 3)
w.S2.EP <- round(unique(subset(data_with_results, treat == "EP")$w1), 3)
w.S2.W <- round(unique(subset(data_with_results, treat == "W")$w1), 3)
w.S2.ctrl <- round(unique(subset(data_with_results, treat == "C")$w1), 3)

se.S2.RP <- round(unique(subset(data_with_results, treat == "RP")$se1), 3)
se.S2.EP <- round(unique(subset(data_with_results, treat == "EP")$se1), 3)
se.S2.W <- round(unique(subset(data_with_results, treat == "W")$se1), 3)
se.S2.ctrl <- round(unique(subset(data_with_results, treat == "C")$se1), 3)

p2 =  ggplot(subset(data_with_results, treat%in%c("EP","W","RP","C")), aes(x = year1, y = richness,group = treat, color = treat)) +
  geom_point(size = 4, alpha = 0.4,show.legend = F) +
  # geom_errorbar(aes(ymin = richness - se, ymax = richness + se),
  #               width = 5.5, alpha = 0.3,
  #               position = position_dodge(.0), lwd = 1.2
  # ) +
  scale_color_manual(values = c("#1874CD", "#CD2626","#DAA520","#808080"),  limits = c("EP","W","RP","C")) +
  geom_smooth(
    data = data_with_results %>% filter(treat == "RP"),
    aes(linetype = significant1), show.legend = F,
    method = "lm", se = F, linewidth = 2
  ) +
  geom_smooth(
    data = data_with_results %>% filter(treat == "EP"),
    aes(linetype = significant1), show.legend = F,
    method = "lm", se = F, linewidth = 2
  ) +
  geom_smooth(
    data = data_with_results %>% filter(treat == "W"),
    aes(linetype = significant1), show.legend = F,
    method = "lm", se = F, linewidth = 2
  ) +
  geom_smooth(
    data = data_with_results %>% filter(treat == "C"),
    aes(linetype = significant1), show.legend = F,
    method = "lm", se = F, linewidth = 2
  ) +
  scale_linetype_manual(values = c("significant" = "solid", "not significant" = "dashed")) +
  xlab("Soil depth")+
  ylab("Richness")+
  annotate("text",
           x = 20,
           y = 1330,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S2.RP, digits = 3, format = "f")) %+-% .(formatC(se.S2.RP, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S2.RP, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S2.RP, digits = 3, format = "f")))),
           colour = "#DAA520"
  ) +
  annotate("text",
           x = 20,
           y = 1200,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S2.EP, digits = 3, format = "f")) %+-% .(formatC(se.S2.EP, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S2.EP, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S2.EP, digits = 3, format = "f")))),
           colour = "#1874CD"
  ) +
  annotate("text",
           x = 20,
           y = 1070,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S2.W, digits = 3, format = "f")) %+-% .(formatC(se.S2.W, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S2.W, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S2.W, digits = 3, format = "f")))),
           colour = "#CD2626"
  ) +
  annotate("text",
           x = 12,
           y = 580,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S2.ctrl, digits = 3, format = "f")) %+-% .(formatC(se.S2.ctrl, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S2.ctrl, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S2.ctrl, digits = 3, format = "f")))),
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
p2

### Protist ####
prefix.m = "Protist"
divindex.file = paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03242023 Inner Mongolia TS\\results\\11142023\\alpha diversity\\",prefix.m,"\\",prefix.m,".alpha.div.csv")

prefix = prefix.m
result.data = read.csv(divindex.file, header = T, row.names = 1, sep = ",")
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
  do(tdcm.mean.no.compare(.[, "richness", drop = F], ., rand = 1000,
                          scale.num = F))

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

R.S3.RP <- round(unique(subset(data_with_results, treat == "RP")$R1), 3)
R.S3.EP <- round(unique(subset(data_with_results, treat == "EP")$R1), 3)
R.S3.W <- round(unique(subset(data_with_results, treat == "W")$R1), 3)
R.S3.ctrl <- round(unique(subset(data_with_results, treat == "C")$R1), 3)

P.S3.RP <- round(unique(subset(data_with_results, treat == "RP")$P1), 3)
P.S3.EP <- round(unique(subset(data_with_results, treat == "EP")$P1), 3)
P.S3.W <- round(unique(subset(data_with_results, treat == "W")$P1), 3)
P.S3.ctrl <- round(unique(subset(data_with_results, treat == "C")$P1), 3)

w.S3.RP <- round(unique(subset(data_with_results, treat == "RP")$w1), 3)
w.S3.EP <- round(unique(subset(data_with_results, treat == "EP")$w1), 3)
w.S3.W <- round(unique(subset(data_with_results, treat == "W")$w1), 3)
w.S3.ctrl <- round(unique(subset(data_with_results, treat == "C")$w1), 3)

se.S3.RP <- round(unique(subset(data_with_results, treat == "RP")$se1), 3)
se.S3.EP <- round(unique(subset(data_with_results, treat == "EP")$se1), 3)
se.S3.W <- round(unique(subset(data_with_results, treat == "W")$se1), 3)
se.S3.ctrl <- round(unique(subset(data_with_results, treat == "C")$se1), 3)

p3 =  ggplot(subset(data_with_results, treat%in%c("EP","W","RP","C")), aes(x = year1, y = richness,group = treat, color = treat)) +
  geom_point(size = 4, alpha = 0.4,show.legend = F) +
  # geom_errorbar(aes(ymin = richness - se, ymax = richness + se),
  #               width = 5.5, alpha = 0.3,
  #               position = position_dodge(.0), lwd = 1.2
  # ) +
  scale_color_manual(values = c("#1874CD", "#CD2626","#DAA520","#808080"),  limits = c("EP","W","RP","C")) +
  geom_smooth(
    data = data_with_results %>% filter(treat == "RP"),
    aes(linetype = significant1), show.legend = F,
    method = "lm", se = F, linewidth = 2
  ) +
  geom_smooth(
    data = data_with_results %>% filter(treat == "EP"),
    aes(linetype = significant1), show.legend = F,
    method = "lm", se = F, linewidth = 2
  ) +
  geom_smooth(
    data = data_with_results %>% filter(treat == "W"),
    aes(linetype = significant1), show.legend = F,
    method = "lm", se = F, linewidth = 2
  ) +
  geom_smooth(
    data = data_with_results %>% filter(treat == "C"),
    aes(linetype = significant1), show.legend = F,
    method = "lm", se = F, linewidth = 2
  ) +
  scale_linetype_manual(values = c("significant" = "solid", "not significant" = "dashed")) +
  xlab("Soil depth")+
  ylab("Richness")+
  annotate("text",
           x = 20,
           y = 420,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S3.RP, digits = 3, format = "f")) %+-% .(formatC(se.S3.RP, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S3.RP, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S3.RP, digits = 3, format = "f")))),
           colour = "#DAA520"
  ) +
  annotate("text",
           x = 20,
           y = 350,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S3.EP, digits = 3, format = "f")) %+-% .(formatC(se.S3.EP, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S3.EP, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S3.EP, digits = 3, format = "f")))),
           colour = "#1874CD"
  ) +
  annotate("text",
           x = 20,
           y = 280,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S3.W, digits = 3, format = "f")) %+-% .(formatC(se.S3.W, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S3.W, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S3.W, digits = 3, format = "f")))),
           colour = "#CD2626"
  ) +
  annotate("text",
           x = 12,
           y = 650,
           size = 5.5,
           label = bquote(atop("w" ~ "=" ~ .(formatC(w.S3.ctrl, digits = 3, format = "f")) %+-% .(formatC(se.S3.ctrl, digits = 3, format = "f"))~","~ italic(r)^2 ~ "=" ~ .(formatC(R.S3.ctrl, digits = 3, format = "f")) ~ ",",italic(P) ~ "=" ~ .(formatC(P.S3.ctrl, digits = 3, format = "f")))),
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
p3

### plot ####
ggarrange(p1,p2,p3,
          ncol = 2, nrow = 2,
          labels = c("(a)", "(b)","(c)"))
ggsave("richness.depth.TS.pdf", width = 12, height = 12, units = "in")
