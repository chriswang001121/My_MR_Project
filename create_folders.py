import os

# 当前目录为 My_MR_Project，不再创建 My_MR_Project 文件夹
project_root = os.getcwd()

# 定义要创建的文件夹列表
folders = [
    "config",                      # 配置文件文件夹
    "data/raw",                    # 原始数据文件夹
    "data/processed",              # 处理后的数据文件夹
    "scripts",                     # 存放脚本
    "results/figures",             # 图表结果
    "results/summary_stats",       # 统计结果
    "docs",                        # 项目文档
    "tests"                        # 测试文件夹
]

# 创建文件夹
for folder in folders:
    path = os.path.join(project_root, folder)
    os.makedirs(path, exist_ok=True)
    print(f"Created folder: {path}")

# 创建文件（如config.yaml、README.md等）
files = {
    "config/config.yaml": "# 通用配置文件\n",
    "config/data_paths.yaml": "# 数据路径的配置文件\n",
    "README.md": "# 项目说明文档\n",
    ".gitignore": "# 忽略规则\n",
    "requirements.txt": "# 依赖列表\n"
}

# 创建并写入内容到文件
for file, content in files.items():
    file_path = os.path.join(project_root, file)
    with open(file_path, 'w') as f:
        f.write(content)
        print(f"Created file: {file_path}")
