# 安装 CRAN 包
install_cran_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing CRAN package:", pkg))
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# 安装 GitHub 包
install_github_packages <- function(packages) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  for (pkg in packages) {
    message(paste("Installing GitHub package:", pkg))
    remotes::install_github(pkg)
  }
}

# 安装 Bioconductor 包
install_bioc_packages <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing Bioconductor package:", pkg))
      BiocManager::install(pkg)
    }
  }
}

# 定义要安装的包
cran_packages <- c("data.table", "tidyverse")
github_packages <- c("MRCIEU/TwoSampleMR", "mrcieu/ieugwasr", "mrcieu/gwasglue")
bioc_packages <- c("VariantAnnotation")

# 执行安装函数
install_cran_packages(cran_packages)
install_github_packages(github_packages)
install_bioc_packages(bioc_packages)

message("All packages installed successfully!")
