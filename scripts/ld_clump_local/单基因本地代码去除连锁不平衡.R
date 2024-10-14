# 本地代码去除连锁不平衡----

library(data.table)
source("ld_clump.R") 
source("ld_matrix.R") 
source("afl2.r")
source("api.R")
source("backwards.R")
source("query.R")
source("utils-pipe.R")
source("variants.R")
source("zzz.R")
library(TwoSampleMR)
library(ieugwasr)

# 读取单基因eQTL数据
expo_rt=fread("ATAD3B.txt", header = T)
head(expo_rt)

# 进行显著性检测
expo_rt=expo_rt[expo_rt$p<1e-5,]#5e-8

# 选取id列和p值列并重命名
expo_rt2=expo_rt[,c("rsids2","p")]
colnames(expo_rt2)=c("rsid", "pval")

# 本地代码去除连锁不平衡，需要完整路径
clumdf <- ld_clump_local(dat = expo_rt2, clump_kb = 100, clump_r2 = 0.1,clump_p=1,
                     bfile ="C:/Users/84241/Desktop/网课/孟德尔随机化/生信狂人/生信狂人-109_基因孟德尔随机化课程/codes_生信狂人_基因孟德尔随机化课程/单基因MR分析和敏感性分析/基因去除连锁不平衡/data_maf0/data_maf0.01_rs_ref", 
                     plink_bin = "C:/Users/84241/Desktop/网课/孟德尔随机化/生信狂人/生信狂人-109_基因孟德尔随机化课程/codes_生信狂人_基因孟德尔随机化课程/单基因MR分析和敏感性分析/基因去除连锁不平衡/plink_win64_20231018/plink.exe")

# 提取符合条件的基因的全部数据并写出
expo_rt3=expo_rt[which(expo_rt$rsids2%in%clumdf$rsid),]
write.table(expo_rt3,"expo_rt.txt",row.names = F,sep = "\t",quote = F)



