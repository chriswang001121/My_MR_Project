# ----install package----

# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")

# install.packages("data.table")

# install.packages("tidyverse")

# install.packages("devtools")
# devtools::install_github("mrcieu/ieugwasr", force = TRUE)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("VariantAnnotation")

# devtools::install_github("mrcieu/gwasglue",force = TRUE)


# ----library package----

library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(ieugwasr)
library(VariantAnnotation)
library(gwasglue)
library(MRPRESSO)

# ----input exposure gwas summary data----

## 1.deCODE_protein

ExposureGwas <- fread("6895_1_TFRC_TR.txt.gz",sep="\t",data.table = F) # input
head(ExposureGwas)

Exposure_P_MAF <- filter(ExposureGwas,Pval<5e-08 & ImpMAF>0.01) # filtered by pvalue and maf
write.csv(Exposure_P_MAF, file=ExposurePvalue, row.names=F)

Exposure_Data<-read_exposure_data(filename="Exposure_P_MAF.csv",
                                 sep = ",",
                                 snp_col = "rsids",
                                 beta_col = "Beta",
                                 se_col = "SE",
                                 effect_allele_col = "effectAllele",
                                 other_allele_col = "otherAllele",
                                 eaf_col = "ImpMAF",
                                 pval_col = "Pval",
                                 samplesize_col = "N",
                                 chr_col="Chrom",
                                 pos_col = "Pos",
                                 clump = F)

## 2.UK Biobank

ExposureGwas <- fread(file = '100015_raw.gwas.imputed_v3.both_sexes.tsv',sep="\t",data.table = F) # input
head(ExposureGwas)

temp1 <- strsplit(ExposureGwas$variant,split = ":",fixed = T)

ExposureGwas <- ExposureGwas %>%
  dplyr::mutate(
    chr = sapply(temp1,function(x){x[1]}),
    pos = sapply(temp1,function(x){x[2]}),
    ref = sapply(temp1,function(x){x[3]}),
    alt = sapply(temp1,function(x){x[4]}),
    eaf = ifelse(minor_AF < 0.5, minor_AF, 1 - minor_AF))

# position to rsid
# posTorsid <- fread(file = 'variants.tsv',sep="\t",data.table = F)
# head(posTorsid)
# tail(posTorsid)
# 
# posTorsid <- posTorsid %>%
#   dplyr::select(variant,rsid)
# 
# write.csv(posTorsid,file = 'uk_biobank_posTosid.csv',row.names = F)

posTorsid <- fread(file = 'uk_biobank_posTosid.csv',sep = ',',data.table = F)

ExposureGwas <- ExposureGwas %>%
  left_join(posTorsid,by = 'variant') 

Exposure_P_MAF <- filter(ExposureGwas,pval<5e-6 & eaf>1e-100) # filtered by pvalue and maf

write.csv(Exposure_P_MAF, file="Exposure_P_MAF.csv", row.names=F)

Exposure_Data <- read_exposure_data(filename="Exposure_P_MAF.csv",
                                  sep = ",",
                                  snp_col = "rsid",
                                  beta_col = "beta",
                                  se_col = "se",
                                  effect_allele_col = "alt",
                                  other_allele_col = "ref",
                                  eaf_col = "eaf",
                                  pval_col = "pval",
                                  samplesize_col = "n_complete_samples",
                                  chr_col="chr",
                                  pos_col = "pos",
                                  clump = F)
Exposure_Data$exposure <- 'Vitamin C'

Exposure_Data$id.exposure <- 'Vitamin C'


### 1.1 ld_clump_online
# Exposure_Filtered <- clump_data(Exposure_Data, clump_kb=10000, clump_r2=0.001)

### 1.2 ld_clump_local
source("ld_clump_local/ld_clump.R") 
source("ld_clump_local/ld_matrix.R") 
source("ld_clump_local/afl2.r")
source("ld_clump_local/api.R")
source("ld_clump_local/backwards.R")
source("ld_clump_local/query.R")
source("ld_clump_local/utils-pipe.R")
source("ld_clump_local/variants.R")
source("ld_clump_local/zzz.R")

Exposure_Data_Id_Pvalue <- Exposure_Data %>%
  dplyr::select(SNP,pval.exposure) %>% # select rsids and pvalue
  dplyr::rename(rsid='SNP',pval="pval.exposure")

clumdf <- ld_clump_local(dat = Exposure_Data_Id_Pvalue, clump_kb = 10000, clump_r2 = 0.001,clump_p=1,
                         bfile ="C:/Users/84241/Desktop/MR_Omics_Database/MR_Own/ld_clump_local/data_maf0/data_maf0.01_rs_ref", 
                         plink_bin = "C:/Users/84241/Desktop/MR_Omics_Database/MR_Own/ld_clump_local/plink_win64_20231018/plink.exe")

Exposure_Filtered=Exposure_Data[which(Exposure_Data$SNP%in%clumdf$rsid),]

write.csv(Exposure_Filtered, file="ExposureFiltered.csv", row.names=F)

Exposure_Filtered_F <- Exposure_Filtered %>%
  mutate(R2 = 2*beta.exposure*beta.exposure*eaf.exposure*(1-eaf.exposure)/(2*beta.exposure*beta.exposure*eaf.exposure*(1-eaf.exposure)+2*se.exposure*se.exposure*samplesize.exposure*eaf.exposure*(1-eaf.exposure)), # calculate R2
         F = R2*(samplesize.exposure-2)/(1-R2)) %>% # calculate F
  filter(F>10)
 
write.csv(Exposure_Filtered_F, file="Exposure_Filtered_F.csv", row.names=F)

## 3.ieu

Exposure_Data <- extract_instruments(outcomes = "ukb-b-20261", p1 = 5e-08,clump = F,r2=0.001,kb=10000, p2 = 5e-08,access_token = NULL)

Exposure_Data$exposure <- "Ever smoked"

source("ld_clump_local/ld_clump.R") 
source("ld_clump_local/ld_matrix.R") 
source("ld_clump_local/afl2.r")
source("ld_clump_local/api.R")
source("ld_clump_local/backwards.R")
source("ld_clump_local/query.R")
source("ld_clump_local/utils-pipe.R")
source("ld_clump_local/variants.R")
source("ld_clump_local/zzz.R")

Exposure_Data_Id_Pvalue <- Exposure_Data %>%
  dplyr::select(SNP,pval.exposure) %>% # select rsids and pvalue
  dplyr::rename(rsid='SNP',pval="pval.exposure")

clumdf <- ld_clump_local(dat = Exposure_Data_Id_Pvalue, clump_kb = 10000, clump_r2 = 0.001,clump_p=1,
                         bfile ="C:/Users/84241/Desktop/MR_Omics_Database/MR_Own/ld_clump_local/data_maf0/data_maf0.01_rs_ref", 
                         plink_bin = "C:/Users/84241/Desktop/MR_Omics_Database/MR_Own/ld_clump_local/plink_win64_20231018/plink.exe")

Exposure_Filtered <- Exposure_Data[which(Exposure_Data$SNP%in%clumdf$rsid),]

write.csv(Exposure_Filtered, file="ExposureFiltered.csv", row.names=F)

Exposure_Filtered_F <- Exposure_Filtered %>%
  mutate(R2 = 2*beta.exposure*beta.exposure*eaf.exposure*(1-eaf.exposure)/(2*beta.exposure*beta.exposure*eaf.exposure*(1-eaf.exposure)+2*se.exposure*se.exposure*samplesize.exposure*eaf.exposure*(1-eaf.exposure)), # calculate R2
         F = R2*(samplesize.exposure-2)/(1-R2)) %>% # calculate F
  filter(F>10)

write.csv(Exposure_Filtered_F, file="Exposure_Filtered_F.csv", row.names=F)

EXP <- Exposure_Filtered_F




# ----input outcome gwas summary data----

## 1.1 local txt file or txt.gz
OutcomeFile <- "finngen_R10_C3_SQUOMOUS_CELL_CARCINOMA_SKIN_EXALLC.gz" # file name
ExposurePvalue <- "Exposure_P_MAF.csv"
ExposureFiltered <- "ExposureFiltered.csv" # output file name
OutcomeGwas <- fread(OutcomeFile,data.table = F)
head(OutcomeGwas)
OutcomeGwas <- OutcomeGwas %>%
  dplyr::rename(
    effect_allele.outcome = 'alt',
    other_allele.outcome = 'ref', 
    SNP = 'rsids',
    pval.outcome = 'pval',
    beta.outcome = 'beta',
    se.outcome = 'sebeta',
    eaf.outcome = 'af_alt') %>%
  mutate(id.outcome = 'SCC',
         outcome = "SCC",
         samplesize.outcome = 317724)

## 1.2 local txt file or txt.gz
OutcomeGwas <- read_outcome_data(snps=exposure_dat$SNP,
                               filename=OutcomeFile, 
                               sep = "\t",
                               snp_col = "rsids",
                               beta_col = "beta",
                               se_col = "sebeta",
                               effect_allele_col = "alt",
                               other_allele_col = "ref",
                               pval_col = "pval",
                               eaf_col = "af_alt")

## 1.3 local vcf file

OutcomeGwas <- readVcf("ieu-a-7.vcf.gz" )

OutcomeGwas <- gwasvcf_to_TwoSampleMR(OutcomeGwas,type="outcome")

head(OutcomeGwas)

## 1.4 ieu

OUT <- extract_outcome_data(snps=EXP$SNP,outcomes="ieu-a-987",proxies=T,maf_threshold = 0.01,access_token = NULL)
OUT$outcome<-"Lung cancer"
OUT <- OUT[!duplicated(OUT$SNP),]


# ----harmonise exposure and outcome----

total <- merge(OutcomeGwas,Exposure_Filtered_F,by.x="SNP",by.y="SNP",all = F)
total <- subset(total,pval.outcome>5e-08)
total <- total[!duplicated(total$SNP),]

EXP <- total[,c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure", "beta.exposure","se.exposure", "pval.exposure","id.exposure","exposure","samplesize.exposure")]
OUT <- total[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome", "beta.outcome","se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

MR_dat <- harmonise_data(exposure_dat=EXP,outcome_dat=OUT,action=2)
write.csv(OUT, file="MR_Outcome.csv", row.names=F)

MR_dat_keep <- MR_dat[MR_dat$mr_keep == "TRUE",]
write.csv(MR_dat_keep, file="MR_SNP.csv", row.names=F)


# ----Confound----
# library(phenoscanner)
# dim(MR_dat)[1] # number of SNP
# PhenoScan <- phenoscanner(snpquery=MR_dat$SNP[1:82],pvalue = 5e-08)
# # 100 once
# write.csv(PhenoScan$result,file="PhenoScan.csv")

# MR_dat <- MR_dat %>%
#   dplyr::filter(SNP %in% SNP_Confound)

# ----Steiger filtering----

MR_dat_steiger <- steiger_filtering(MR_dat) |>
  subset(steiger_dir==TRUE)

MR_dat_steiger <- MR_dat_steiger[!duplicated(MR_dat_steiger$SNP),]

MR_dat <- MR_dat_steiger

# ----MR & heterogeneity & pleiotropy----

MR_Result <- mr(MR_dat)

# MR PRESSO for Outlier value
presso=run_mr_presso(MR_dat)
write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))

# MR PRESSO
library(MRPRESSO)
MRPRESSO <- mr_presso(BetaOutcome = "beta.outcome",BetaExposure="beta.exposure", 
          SdOutcome = "se.outcome", SdExposure = "se.exposure",
          OUTLIERtest = T, 
          DISTORTIONtest = T,
          data=MR_dat)

# MRPRESSO$`MR-PRESSO results`$`Global Test`,若p < 0.05,则具有水平多效性
# MRPRESSO$`MR-PRESSO results`$`Outlier Test`，若p < 0.05,则具有水平多效性
# 水平多效性指SNP可以不经过exposure影响结局

# 对结果进行OR值的计算
MR_Result_OR <- generate_odds_ratios(MR_Result)
write.csv(MR_Result_OR, file="MR_Result_OR.csv", row.names=F)

# 异质性分析
MR_heterogeneity <- mr_heterogeneity(MR_dat)
write.csv(MR_heterogeneity , file="MR_heterogeneity .csv", row.names=F)

# 多效性检验
MR_pleiotropy=mr_pleiotropy_test(MR_dat)
write.csv(MR_pleiotropy, file="MR_pleiotropy.csv", row.names=F)

# 绘制散点图
pdf(file="MR_Scatter_Plot.pdf", width=7.5, height=7)
mr_scatter_plot(MR_Result, MR_dat)
dev.off()

# 森林图
pdf(file="MR_Forest.pdf", width=7, height=5.5)
mr_forest_plot(mr_singlesnp(MR_dat))
dev.off()

# 漏斗图
pdf(file="MR_Funnel_Plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = mr_singlesnp(MR_dat))
dev.off()

# 留一法敏感性分析
pdf(file="MR_Leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(MR_dat))
dev.off()
