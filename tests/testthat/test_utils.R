library(testthat)
source("scripts/MR_Pipeline_Analysis/utils.R")

# ---- Test load_libraries() ----
test_that("load_libraries works correctly", {
  # 测试包是否能正确加载，而不会出现错误
  expect_silent(load_libraries())
})

# ---- Test process_exposure_data() with decode_protein data ----
test_that("process_exposure_data works with decode_protein data", {
  
  # 设置文件路径，使用本地的 decode_protein 文件
  file_path <- "data/raw/6895_1_TFRC_TR.txt.gz"  # 这是你提供的本地文件路径
  
  # 处理暴露数据
  processed_data <- process_exposure_data(
    exposure_source = "decode_protein",
    file_path = file_path,
    exposure_name = "TFRC Protein",
    exposure_id = "6895_TFRC",
    pval_threshold = 5e-08,
    maf_threshold = 0.01,
    output_dir = "data/processed"
  )
  
  # 检查暴露名称是否一致（整列比较）
  expect_true(all(processed_data$exposure == "TFRC Protein"))
  
  # 检查暴露ID是否一致（整列比较）
  expect_true(all(processed_data$id.exposure == "6895_TFRC"))
  
  # 确保返回数据不为空
  expect_true(nrow(processed_data) > 0)
  
  # 检查是否生成了正确的输出文件
  expect_true(file.exists("data/processed/Exposure_P_MAF.csv"))
})

# ---- Test process_exposure_data() with UK_Biobank data ----
test_that("process_exposure_data works with UK_Biobank data", {
  
  # 设置 UK Biobank 数据的本地路径
  file_path <- "data/raw/C43.gwas.imputed_v3.both_sexes.tsv"  # 这是你提供的本地文件路径
  
  # 处理暴露数据
  processed_data <- process_exposure_data(
    exposure_source = "UK_Biobank",
    file_path = file_path,
    exposure_name = "UK Biobank GWAS",
    exposure_id = "C43",
    pval_threshold = 5e-08,
    maf_threshold = 0.01,
    output_dir = "data/processed"
  )
  
  # 检查暴露名称是否一致（整列比较）
  expect_true(all(processed_data$exposure == "UK Biobank GWAS"))
  
  # 检查暴露ID是否一致（整列比较）
  expect_true(all(processed_data$id.exposure == "C43"))
  
  # 确保返回数据不为空
  expect_true(nrow(processed_data) > 0)
  
  # 检查是否生成了正确的输出文件
  expect_true(file.exists("data/processed/Exposure_P_MAF.csv"))
})

# ---- Test process_exposure_data() with IEU data ----
test_that("process_exposure_data works with IEU data", {
  skip("由于 extract_instruments 依赖外部数据库，该测试跳过或模拟")
  
  # 测试 extract_instruments 的代码片段可以直接根据需要补充
  # file_path <- "ukb-b-20261"  # 通过 IEU API 提取的数据集 ID
  
  # processed_data <- process_exposure_data(
  #   exposure_source = "ieu",
  #   file_path = file_path,
  #   exposure_name = "Ever Smoked",
  #   exposure_id = "ieu_20261",
  #   pval_threshold = 5e-08,
  #   maf_threshold = 0.01,
  #   output_dir = "data/processed"
  # )
  
  # 检查数据是否正确生成
  # expect_equal(processed_data$exposure, "Ever Smoked")
  # expect_equal(processed_data$id.exposure, "ieu_20261")
  # expect_true(nrow(processed_data) > 0)
})


# test_ld_clumping.R

test_that("Online LD clumping works correctly", {
  # 加载依赖库和测试数据

  load_libraries()
  
  # 模拟的输入数据
  test_data <- data.frame(
    SNP = c("rs12345", "rs67890", "rs24680"),
    pval.exposure = c(5e-07, 1e-06, 4e-07),
    beta.exposure = c(0.02, 0.03, -0.01),
    se.exposure = c(0.01, 0.01, 0.02),
    eaf.exposure = c(0.3, 0.4, 0.25),
    stringsAsFactors = FALSE
  )
  
  # 运行在线 LD Clumping
  Exposure_Filtered <- perform_ld_clumping(
    Exposure_Data = test_data, 
    method = "online", 
    clump_kb = 10000, 
    clump_r2 = 0.001
  )
  
  # 测试返回结果是否有数据
  expect_true(nrow(Exposure_Filtered) > 0)
  # 测试返回结果是否包含必要的列
  expect_true("SNP" %in% colnames(Exposure_Filtered))
})

test_that("Local LD clumping works correctly", {
  # 加载依赖库和测试数据
  
  load_libraries()
  
  # 模拟的输入数据
  test_data <- data.frame(
    SNP = c("rs12345", "rs67890", "rs24680"),
    pval.exposure = c(5e-07, 1e-06, 4e-07),
    beta.exposure = c(0.02, 0.03, -0.01),
    se.exposure = c(0.01, 0.01, 0.02),
    eaf.exposure = c(0.3, 0.4, 0.25),
    stringsAsFactors = FALSE
  )
  
  # 定义本地 PLINK 文件和可执行文件路径
  plink_bfile_path <- "F:/Projects/My_MR_Project/data/raw/plink_ref/data_maf0.01_rs_ref"
  plink_bin_path <- "F:/Projects/My_MR_Project/tools/plink_win64_20231018/plink.exe"
  
  # 运行本地 LD Clumping
  Exposure_Filtered <- perform_ld_clumping(
    Exposure_Data = test_data, 
    method = "local", 
    clump_kb = 10000, 
    clump_r2 = 0.001, 
    plink_bfile = plink_bfile_path, 
    plink_bin = plink_bin_path
  )
  
  # 测试返回结果是否有数据
  expect_true(nrow(Exposure_Filtered) > 0)
  # 测试返回结果是否包含必要的列
  expect_true("SNP" %in% colnames(Exposure_Filtered))
})


library(testthat)
library(data.table)
library(dplyr)
library(VariantAnnotation)
library(TwoSampleMR)


# 测试 process_outcome_data 函数
test_that("process_outcome_data works correctly for finngen source", {
  
  # 设置测试的输入文件路径
  file_path <- "data/raw/finngen_R10_C3_SQUOMOUS_CELL_CARCINOMA_SKIN_EXALLC.gz"
  
  # 调用处理函数
  outcome_data <- process_outcome_data(
    outcome_source = "finngen",
    file_path = file_path,
    outcome_name = "SCC",
    outcome_id = "SCC_ID"
  )
  
  # 断言数据加载成功并且有内容
  expect_s3_class(outcome_data, "data.frame")
  expect_gt(nrow(outcome_data), 0)
  
  # 验证列是否存在
  expect_true(all(c("SNP", "effect_allele.outcome", "other_allele.outcome",
                    "pval.outcome", "beta.outcome", "se.outcome", "eaf.outcome") %in% colnames(outcome_data)))
  
  # 验证 outcome 名字和 ID 是否正确
  expect_equal(outcome_data$outcome[1], "SCC")
  expect_equal(outcome_data$id.outcome[1], "SCC_ID")
})

# 测试 process_outcome_data 函数 for ieu source
test_that("process_outcome_data works correctly for ieu source", {
  
  # 设置测试的输入文件路径
  file_path <- "data/raw/ieu-b-4969.vcf"
  
  # 调用处理函数
  outcome_data <- process_outcome_data(
    outcome_source = "ieu",
    file_path = file_path,
    outcome_name = "Lung Cancer",
    outcome_id = "LC_ID"
  )
  
  # 断言数据加载成功并且有内容
  expect_s3_class(outcome_data, "data.frame")
  expect_gt(nrow(outcome_data), 0)
  
  # 验证列是否存在
  expect_true(all(c("SNP", "effect_allele.outcome", "other_allele.outcome",
                    "pval.outcome", "beta.outcome", "se.outcome", "eaf.outcome") %in% colnames(outcome_data)))
  
  # 验证 outcome 名字和 ID 是否正确
  expect_equal(outcome_data$outcome[1], "Lung Cancer")
  expect_equal(outcome_data$id.outcome[1], "LC_ID")
})

