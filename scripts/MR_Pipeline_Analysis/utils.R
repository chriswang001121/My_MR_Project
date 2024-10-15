# utils.R


# 加载所需R包
load_libraries <- function() {
  library(TwoSampleMR)
  library(data.table)
  library(tidyverse)
  library(ieugwasr)
  library(VariantAnnotation)
  library(gwasglue)
  library(MRPRESSO)
  library(stringr)
}


# 处理暴露数据
process_exposure_data <- function(exposure_source, 
                                            file_path, 
                                            exposure_name = "exposure",
                                            exposure_id = "exposure",
                                            pval_threshold = 5e-08, 
                                            maf_threshold = 0.01,  
                                            pos_to_rsid = "data/processed/uk_biobank_posTosid.csv", 
                                            output_dir = "data/processed/") {
  
  # 确保目标文件夹存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 根据数据来源处理暴露数据
  if (exposure_source == "decode_protein") {
    # 读取 deCODE_protein 数据
    ExposureGwas <- fread(file_path, sep = "\t", data.table = FALSE)
    Exposure_P_MAF <- dplyr::filter(ExposureGwas, 
                             Pval < pval_threshold & ImpMAF > maf_threshold)
    
    # 构建输出文件路径并保存到 data/processed 目录
    output_file <- file.path(output_dir, "Exposure_P_MAF.csv")
    write.csv(Exposure_P_MAF, file = output_file, row.names = FALSE)
    
    # 读取并格式化为标准MR暴露数据格式
    Exposure_Data <- read_exposure_data(
      filename = output_file,
      sep = ",",
      snp_col = "rsids",
      beta_col = "Beta",
      se_col = "SE",
      effect_allele_col = "effectAllele",
      other_allele_col = "otherAllele",
      eaf_col = "ImpMAF",
      pval_col = "Pval",
      samplesize_col = "N",
      chr_col = "Chrom",
      pos_col = "Pos",
      clump = FALSE
    )
    
  } else if (exposure_source == "UK_Biobank") {
    # 读取 UK Biobank 数据
    ExposureGwas <- fread(file = file_path, sep = "\t", data.table = FALSE)
    
    ExposureGwas <- ExposureGwas %>%
      dplyr::mutate(
        # 使用 str_split_fixed 一次性分割 variant 字段
        parsed = str_split_fixed(variant, ":", 4),
        chr = parsed[, 1],
        pos = as.numeric(parsed[, 2]),
        ref = parsed[, 3],
        alt = parsed[, 4],
        
        # 使用 pmin 简化 eaf 计算逻辑
        eaf = pmin(minor_AF, 1 - minor_AF)
      ) %>%
      dplyr::select(-parsed)  # 最后删除临时列

    # pos到rsid 的转换
    
    posTorsid <- fread(file = pos_to_rsid, sep = ',', data.table = FALSE)
    ExposureGwas <- ExposureGwas %>%
        left_join(posTorsid, by = 'variant')

    
    # 过滤数据
    Exposure_P_MAF <- dplyr::filter(ExposureGwas, 
                             pval < pval_threshold & eaf > maf_threshold)
    
    # 构建输出文件路径并保存到 data/processed 目录
    output_file <- file.path(output_dir, "Exposure_P_MAF.csv")
    write.csv(Exposure_P_MAF, file = output_file, row.names = FALSE)
    
    # 读取并格式化为标准MR暴露数据格式
    Exposure_Data <- read_exposure_data(
      filename = output_file,
      sep = ",",
      snp_col = "rsid",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      eaf_col = "eaf",
      pval_col = "pval",
      samplesize_col = "n_complete_samples",
      chr_col = "chr",
      pos_col = "pos",
      clump = FALSE
    )
    
  } else if (exposure_source == "ieu") {
    # 直接从 IEU 提取数据
    Exposure_Data <- extract_instruments(outcomes = file_path, 
                                         p1 = pval_threshold, 
                                         clump = FALSE, 
                                         r2 = 0.001, 
                                         kb = 10000)
    
  } else {
    stop("Unknown data source")
  }
  
  # 补充暴露名字和ID
  Exposure_Data$exposure <- exposure_name
  Exposure_Data$id.exposure <- exposure_id
  
  return(Exposure_Data)
}


# 去除连锁不平衡
perform_ld_clumping <- function(Exposure_Data, 
                                method = c("online", "local"),
                                clump_kb = 10000, 
                                clump_r2 = 0.001, 
                                plink_bfile = NULL, 
                                plink_bin = NULL) {
  
  method <- match.arg(method)
  
  if (method == "online") {
    # 在线LD Clumping
    cat("Performing online LD clumping...\n")
    
    Exposure_Filtered <- clump_data(Exposure_Data, 
                                    clump_kb = clump_kb, 
                                    clump_r2 = clump_r2)
    
  } else if (method == "local") {
    # 本地LD Clumping
    if (is.null(plink_bfile) || is.null(plink_bin)) {
      stop("Please provide valid PLINK file path and binary for local LD clumping.")
    }
    
    cat("Performing local LD clumping...\n")
    
    Exposure_Data_Id_Pvalue <- Exposure_Data %>%
      dplyr::select(SNP, pval.exposure) %>%
      dplyr::rename(rsid = 'SNP', pval = "pval.exposure")
    
    # 调用本地的LD Clumping函数
    clumdf <- ld_clump_local(dat = Exposure_Data_Id_Pvalue, 
                             clump_kb = clump_kb, 
                             clump_r2 = clump_r2, 
                             clump_p = 1, 
                             bfile = plink_bfile, 
                             plink_bin = plink_bin)
    
    Exposure_Filtered <- Exposure_Data %>%
      dplyr::filter(SNP %in% clumdf$rsid)
    
  }
  
  return(Exposure_Filtered)
}


# 封装 R² 和 F 值的计算函数，直接处理 Exposure_Filtered 数据集
calculate_R2_F <- function(Exposure_Filtered) {
  Exposure_Filtered_F <- Exposure_Filtered %>%
    dplyr::mutate(
      R2 = 2 * beta.exposure^2 * eaf.exposure * (1 - eaf.exposure) / 
        (2 * beta.exposure^2 * eaf.exposure * (1 - eaf.exposure) + 
           2 * se.exposure^2 * samplesize.exposure * eaf.exposure * (1 - eaf.exposure)),  # 计算 R²
      F = R2 * (samplesize.exposure - 2) / (1 - R2)  # 计算 F 值
    ) %>%
    dplyr::filter(F > 10)  # 过滤 F 值大于 10 的行
  
  return(Exposure_Filtered_F)
}


# 封装结局数据处理函数
process_outcome_data <- function(outcome_source, 
                                 file_path,
                                 outcome_name = "outcome", 
                                 outcome_id = "outcome",
                                 samplesize.outcome = 12345) {
  
  # 处理不同来源的结局数据
  if (outcome_source == "finngen") {
    # 读取本地 txt 文件或 txt.gz 文件
    OutcomeGwas <- fread(file_path, data.table = FALSE)
    
    # 重命名列并补充必要的信息
    OutcomeGwas <- OutcomeGwas %>%
      dplyr::rename(
        effect_allele.outcome = 'alt',
        other_allele.outcome = 'ref', 
        SNP = 'rsids',
        pval.outcome = 'pval',
        beta.outcome = 'beta',
        se.outcome = 'sebeta',
        eaf.outcome = 'af_alt'
      )
    
  } else if (outcome_source == "ieu") {
    # 读取本地 VCF 文件
    OutcomeGwas <- readVcf(file_path)
    
    # 转换为 TwoSampleMR 格式
    OutcomeGwas <- gwasvcf_to_TwoSampleMR(OutcomeGwas, type = "outcome")
    
  } else {
    stop("未知的数据来源类型。")
  }
  
  # 补充 id.outcome 和 outcome 列
  OutcomeGwas <- OutcomeGwas %>%
    dplyr::mutate(
      id.outcome = outcome_id,
      outcome = outcome_name,
      samplesize.outcome = samplesize.outcome
    )
  
  return(OutcomeGwas)
}


# harmonise函数封装
harmonise_exposure_outcome <- function(OutcomeGwas, 
                                       ExposureFiltered, 
                                       pval_threshold = 5e-08, 
                                       output_dir = "data/processed") {
  
  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 合并暴露和结局数据
  total <- merge(
    OutcomeGwas, ExposureFiltered, 
    by.x = "SNP", by.y = "SNP", 
    all = FALSE
  )
  
  # 过滤结局数据的p值
  total <- subset(total, pval.outcome > pval_threshold)
  
  # 去除重复SNP
  total <- total[!duplicated(total$SNP), ]
  
  # 选择暴露数据字段
  EXP <- total %>%
    dplyr::select(
      SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure,
      beta.exposure, se.exposure, pval.exposure, id.exposure,
      exposure, samplesize.exposure
    )
  
  # 选择结局数据字段
  OUT <- total %>%
    dplyr::select(
      SNP, effect_allele.outcome, other_allele.outcome, eaf.outcome,
      beta.outcome, se.outcome, pval.outcome, id.outcome,
      outcome, samplesize.outcome
    )
  
  # 调和暴露和结局数据
  MR_dat <- harmonise_data(
    exposure_dat = EXP, 
    outcome_dat = OUT, 
    action = 2
  )
  
  # 保存结局数据为CSV文件
  write.csv(
    OUT, 
    file = file.path(output_dir, "MR_Outcome.csv"), 
    row.names = FALSE
  )
  
  # 过滤掉不符合MR分析标准的数据
  MR_dat_keep <- MR_dat[MR_dat$mr_keep == "TRUE", ]
  
  # 保存处理过的MR SNP数据为CSV文件
  write.csv(
    MR_dat_keep, 
    file = file.path(output_dir, "MR_SNP.csv"), 
    row.names = FALSE
  )
  
  # 返回调和后的MR数据
  return(MR_dat)
}


# 混杂因子分析函数
perform_confounding_analysis <- function(MR_dat, 
                                         output_dir = "data/processed", 
                                         snp_limit = 100) {
  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 使用 PhenoScanner 查询SNP信息，每次最多查询 snp_limit 个 SNP
  library(phenoscanner)
  
  num_snps <- nrow(MR_dat)
  snp_batches <- split(MR_dat$SNP, 
                       ceiling(seq_along(MR_dat$SNP) / snp_limit))
  
  pheno_scan_results <- list()
  
  for (batch in snp_batches) {
    pheno_scan <- phenoscanner(snpquery = batch, pvalue = 5e-08)
    pheno_scan_results[[length(pheno_scan_results) + 1]] <- pheno_scan$result
  }
  
  # 合并所有查询结果
  PhenoScan <- do.call(rbind, pheno_scan_results)
  
  # 保存结果至输出目录
  output_file <- file.path(output_dir, "PhenoScan.csv")
  write.csv(PhenoScan, output_file, row.names = FALSE)
  
  return(PhenoScan)
}


# Steiger 方向性过滤函数
perform_steiger_filtering <- function(MR_dat) {
  # 执行 Steiger 过滤，过滤掉方向性不一致的 SNP
  MR_dat_steiger <- steiger_filtering(MR_dat) %>%
    dplyr::filter(steiger_dir == TRUE) %>%
    dplyr::distinct(SNP, .keep_all = TRUE)  # 去重
  
  return(MR_dat_steiger)
}


run_mr_analysis <- function(MR_dat, output_dir = "results/summary_stats/Pipeline/") {
  MR_Result <- mr(MR_dat)
  
  # 保存 OR 值的计算结果
  MR_Result_OR <- generate_odds_ratios(MR_Result)
  output_file <- file.path(output_dir, "MR_Result_OR.csv")
  write.csv(MR_Result_OR, file = output_file, row.names = FALSE)  
  return(MR_Result)
}


run_heterogeneity_analysis <- function(MR_dat, output_dir = "results/summary_stats/Pipeline/") {
  MR_heterogeneity <- mr_heterogeneity(MR_dat)
  output_file <- file.path(output_dir, "MR_heterogeneity.csv")
  write.csv(MR_heterogeneity, file = output_file, row.names = FALSE)
  return(MR_heterogeneity)
}


run_pleiotropy_test <- function(MR_dat, output_dir = "results/summary_stats/Pipeline/") {
  MR_pleiotropy <- mr_pleiotropy_test(MR_dat)
  output_file <- file.path(output_dir, "MR_pleiotropy.csv")
  write.csv(MR_pleiotropy, file = output_file, row.names = FALSE) 
  return(MR_pleiotropy)
}


plot_mr_results <- function(MR_Result, 
                            MR_dat, 
                            plot_type = c("scatter", 
                                          "forest", 
                                          "funnel", 
                                          "leaveoneout"), 
                            output_dir = "results/figures/Pipeline/") {
  
  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 根据图类型生成不同的图表
  for (type in plot_type) {
    if (type == "scatter") {
      file_name <- file.path(output_dir, "MR_Scatter_Plot.pdf")
      pdf(file = file_name, width = 7.5, height = 7)
      mr_scatter_plot(MR_Result, MR_dat)
      dev.off()
    }
    
    if (type == "forest") {
      file_name <- file.path(output_dir, "MR_Forest.pdf")
      pdf(file = file_name, width = 7, height = 5.5)
      mr_forest_plot(mr_singlesnp(MR_dat))
      dev.off()
    }
    
    if (type == "funnel") {
      file_name <- file.path(output_dir, "MR_Funnel_Plot.pdf")
      pdf(file = file_name, width = 7, height = 6.5)
      mr_funnel_plot(singlesnp_results = mr_singlesnp(MR_dat))
      dev.off()
    }
    
    if (type == "leaveoneout") {
      file_name <- file.path(output_dir, "MR_Leaveoneout.pdf")
      pdf(file = file_name, width = 7, height = 5.5)
      mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(MR_dat))
      dev.off()
    }
  }
}


run_mr_presso_analysis <- function(MR_dat, 
                                   beta_outcome_col = "beta.outcome", 
                                   beta_exposure_col = "beta.exposure", 
                                   se_outcome_col = "se.outcome", 
                                   se_exposure_col = "se.exposure", 
                                   outlier_test = TRUE, 
                                   distortion_test = TRUE,
                                   output_dir = "results/summary_stats/Pipeline/") {
  
  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 运行 MR-PRESSO 分析
  presso_result <- mr_presso(BetaOutcome = beta_outcome_col,
                             BetaExposure = beta_exposure_col,
                             SdOutcome = se_outcome_col,
                             SdExposure = se_exposure_col,
                             OUTLIERtest = outlier_test,
                             DISTORTIONtest = distortion_test,
                             data = MR_dat, 
                             NbDistribution = 1000)
  
  # 保存 MR-PRESSO 结果到文件
  global_test_path <- file.path(output_dir, "MR-PRESSO_Global_Test.csv")
  outlier_test_path <- file.path(output_dir, "MR-PRESSO_Outlier_Test.csv")
  
  if (!is.null(presso_result[[2]]$`Global Test`)) {
    write.csv(presso_result[[2]]$`Global Test`, 
              file = global_test_path, row.names = FALSE)
  }
  
  if (!is.null(presso_result[[2]]$`Outlier Test`)) {
    write.csv(presso_result[[2]]$`Outlier Test`, 
              file = outlier_test_path, row.names = FALSE)
  }
  
  # 返回结果
  return(presso_result)
}
