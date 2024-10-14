import os

def generate_file_tree(root_dir):
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # 跳过.git和.venv文件夹
        dirnames[:] = [d for d in dirnames if d not in ['.git', '.venv']]
        
        # 打印文件夹路径
        level = dirpath.replace(root_dir, '').count(os.sep)
        indent = ' ' * 4 * level
        print(f'{indent}{os.path.basename(dirpath)}/')
        
        # 打印文件路径
        subindent = ' ' * 4 * (level + 1)
        for f in filenames:
            print(f'{subindent}{f}')

if __name__ == "__main__":
    project_root = os.path.dirname(os.path.abspath(__file__))  # 获取当前项目的根目录
    generate_file_tree(project_root)
