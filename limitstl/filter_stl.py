import trimesh

def filter_largest_component(input_file, output_file):
    """
    读取STL文件，保留最大的连通组件，移除其他部分，并保存为新STL文件。
    
    :param input_file: 输入STL文件路径
    :param output_file: 输出STL文件路径
    """
    # 加载STL文件
    mesh = trimesh.load(input_file)
    
    # 分离连通组件（不要求水密）
    components = mesh.split(only_watertight=False)
    
    if not components:
        print("No components found in the STL file.")
        return
    
    # 找到最大的组件（基于面数）
    largest_component = max(components, key=lambda c: len(c.faces))
    
    # 保存最大的组件
    largest_component.export(output_file)
    
    print(f"Filtered STL saved to {output_file}. Original components: {len(components)}, Largest has {len(largest_component.faces)} faces.")

# 示例使用（替换为您的STL文件路径）
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python filter_stl.py <input.stl> <output.stl>")
        sys.exit(1)
    input_stl = sys.argv[1]
    output_stl = sys.argv[2]
    filter_largest_component(input_stl, output_stl)
