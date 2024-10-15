
source("scripts/MR_Pipeline_Analysis/utils.R")

# ----加载所需的 R 包----
load_libraries()

# ----处理暴露数据----
Exposure_Data <- process_exposure_data(
  exposure_source = "decode_protein",  # 可选： "UK_Biobank", "ieu"
  file_path = "data/raw/6895_1_TFRC_TR.txt.gz" ,
  exposure_name = "TFRC",
  exposure_id = "6895_1_TFRC_TR",
  pval_threshold = 5e-08,
  maf_threshold = 0.01
)

# ----去除连锁不平衡 (LD Clumping)----
# 如果你使用在线 LD clumping
Exposure_Filtered <- perform_ld_clumping(
  Exposure_Data = Exposure_Data, 
  method = "online",  # 或 "local"
  clump_kb = 10000, 
  clump_r2 = 0.001
)

# 如果使用本地 LD clumping，需要提供 PLINK 文件和可执行文件的路径
Exposure_Filtered <- perform_ld_clumping(
  Exposure_Data = Exposure_Data, 
  method = "local", 
  clump_kb = 10000, 
  clump_r2 = 0.001, 
  plink_bfile = "F:/Projects/My_MR_Project/data/raw/plink_ref/data_maf0.01_rs_ref",
  plink_bin = "F:/Projects/My_MR_Project/tools/plink_win64_20231018/plink.exe"
)

# ----计算 R² 和 F 值----
Exposure_Filtered_F <- calculate_R2_F(Exposure_Filtered)

# ----处理结局数据 (Outcome Data)----
Outcome_Data <- process_outcome_data(
  outcome_source = "finngen",  # 或 "ieu"
  file_path <- "data/raw/finngen_R10_C3_SQUOMOUS_CELL_CARCINOMA_SKIN_EXALLC.gz",
  outcome_name = "SCC",
  outcome_id = "SCC",
  samplesize.outcome = 22215
)

# ----调和暴露和结局数据 (Harmonise Exposure and Outcome)----
MR_dat <- harmonise_exposure_outcome(
  OutcomeGwas = Outcome_Data, 
  ExposureFiltered = Exposure_Filtered_F,
  pval_threshold = 5e-08,
  output_dir = "data/processed"
)

# ----Steiger 过滤 (Steiger Filtering)----
MR_dat_steiger <- perform_steiger_filtering(MR_dat)

# ----运行 MR 分析 (Mendelian Randomization Analysis)----
MR_Result <- run_mr_analysis(MR_dat_steiger, output_dir = "results/summary_stats/Pipeline")

# ----异质性分析 (Heterogeneity Analysis)----
MR_heterogeneity <- run_heterogeneity_analysis(MR_dat_steiger, output_dir = "results/summary_stats/Pipeline")

# ----多效性检验 (Pleiotropy Test)----
MR_pleiotropy <- run_pleiotropy_test(MR_dat_steiger, output_dir = "results/summary_stats/Pipeline")

# ----绘制 MR 结果图表 (Plotting MR Results)----
plot_mr_results(
  MR_Result = MR_Result, 
  MR_dat = MR_dat_steiger, 
  plot_type = c("scatter", "forest", "funnel", "leaveoneout"), 
  output_dir = "results/figures/Pipeline/"
)

# ----运行 MR-PRESSO 分析 (MR-PRESSO Analysis)----
presso_result <- run_mr_presso_analysis(
  MR_dat = MR_dat_steiger, 
  output_dir = "results/summary_stats/Pipeline/"
)
