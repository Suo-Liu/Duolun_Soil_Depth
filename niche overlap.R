library(vegan)
# protist
com.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\protist.zotu_resample_2961.txt"
tree.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_rooted_tree.nwk"
clas.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_tax.txt"
prefix.m = "Protist"

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"
result = list()
z = 1
# for (i in c("C","W","RP","EP")){
  for (j in c("L1","L2","L3","L4")){
treat<-read.csv(treat.file,header = T,row.names = 1)
# treat = subset(treat, plant.type == "TS" & Layer == "L4" & combined_treat1 %in% c("C","W","RP","EP"))
treat = subset(treat, plant.type == "TS" & Layer %in% j & combined_treat1 %in% c("C","W","RP","EP"))

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


sampc=match.name(rn.list = list(treat = treat,comm = comm))
dim(comm)
comm = sampc$comm
treat = sampc$treat
comm = comm[,colSums(comm)>0]
dim(comm)

spc=match.name(rn.list=list(clas=clas), cn.list = list(comm = comm))
comm=spc$comm
clas=spc$clas

library(spaa)

niche_overlap <- niche.overlap(comm, method = 'levins')
niche_overlap

#上述结果默认以 dist 类型存储，可转换为 matrix 类型输出到本地
niche_overlap <- as.matrix(niche_overlap)
library(ieggr)
niche_overlap.col = dist.3col(niche_overlap)
# mean(niche_overlap.col$dis)
# write.table(niche_overlap, 'niche_overlap.txt', sep = '\t', col.names = NA, quote = FALSE)
result[[z]] = niche_overlap.col$dis
z = z+1
  }

set.seed(123)
niche_overlap_boot <- niche.overlap.boot(dune, method = 'levins', times = 10, quant = c(0.025, 0.975))
head(niche_overlap_boot)



mean(result[[1]])
mean(result[[2]])
mean(result[[3]])
mean(result[[4]])

mean(result[[1]][which(result[[1]]!=0)])
mean(result[[2]][which(result[[2]]!=0)])
mean(result[[3]][which(result[[3]]!=0)])
mean(result[[4]][which(result[[4]]!=0)])


wilcox.test(result[[1]], result[[2]])
wilcox.test(result[[1]], result[[3]])
wilcox.test(result[[1]], result[[4]])
wilcox.test(result[[2]], result[[3]])
wilcox.test(result[[2]], result[[4]])
wilcox.test(result[[3]], result[[4]])


