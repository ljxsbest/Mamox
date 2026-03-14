import re

mechname = 'completeDetailedMech.inp'
PLOGs_dict = {}
TROE_dict = {}

def format_A(s):
    try:
        return "{0:>12.2E}".format(float(s))
    except ValueError:
        return s

def format_n(s):
    try:
        return "{0:>8.2f}".format(float(s))
    except ValueError:
        return s

def format_E(s):
    try:
        return "{0:>8}".format(int(float(s)))
    except ValueError:
        return s

# 读取 PLOG 数据并存入字典
with open('PLOGs.txt', 'r') as f:
    lines = f.readlines()

i = 0
while i < len(lines):
    line = lines[i].strip()
    if line.startswith('Label'):
        label = line.split()[-1]  # 提取 Label
        i += 1
        PLOG_entries = []
        TROE_entries = []
        
        # 如果 Label 是 12，需要存储后面 3 行作为 TROE 数据
        if label == "12":
            TROE_entries = [lines[i].strip(), lines[i+1].strip(), lines[i+2].strip()]
            TROE_dict[label] = TROE_entries
            i += 3  # 跳过 3 行
        else:
            while i < len(lines) and lines[i].strip().startswith('PLOG'):
                plog_line = lines[i].strip()
                parts = plog_line.split('/')
                if len(parts) > 1:
                    values = re.split(r'\s+', parts[1].strip())  # 以空格或缩进分割数值部分
                    if len(values) >= 3:
                        P_value = values[0]  # 压力值 P
                        A_value = float(values[1])  # 频率因子 A
                        n_value = values[2]  # 指数因子 n
                        E_value = values[3] if len(values) > 3 else '0.0'  # 激活能 E，可能为空
                        PLOG_entries.append((plog_line, P_value, A_value, n_value, E_value))  # 存储原始数据
                i += 1
            PLOGs_dict[label] = PLOG_entries
    else:
        i += 1

# 读取机理文件
with open(mechname, 'r') as f:
    mechs = f.readlines()

# 处理并插入 PLOG 数据或 TROE 数据
new_mechs = []
for j, mechline in enumerate(mechs):
    new_mechs.append(mechline)
    mechline = mechline.strip()

    if not mechline or mechline.startswith('!'):
        continue

    if 'RR' in mechline:
        parts = mechline.split('!RR')
        if len(parts) > 1:
            label_parts = parts[1].split()
            if len(label_parts) >= 3:
                Label1 = label_parts[0]  # 反应物标签
                Multiplier = float(label_parts[2])  # 读取 Multiplier 并转换为浮点数

                if Label1 == "12":  # 如果是 Label 12，插入 TROE 数据
                    if Label1 in TROE_dict:
                        updated_TROE = []
                        for line in TROE_dict[Label1]:
                            if line.startswith('LOW'):
                                # 将 LOW 行以 '/' 分隔
                                parts = line.split('/')
                                if len(parts) > 1:
                                    values = re.split(r'\s+', parts[1].strip())  # 以空格或缩进分割
                                    if len(values) > 1:
                                        A_value = float(values[0])  # 原始 A 值
                                        new_A_value = A_value * Multiplier  # 新 A 值
                                        # 重建该行，将 A 值替换为新的 A
                                        updated_line = f"LOW / {new_A_value:.6e} " + " ".join(values[1:]) + " /"
                                        updated_TROE.append(updated_line)
                                    else:
                                        updated_TROE.append(line)  # 如果没有 A 值，保持原样
                            else:
                                updated_TROE.append(line)  # 不是以 LOW 开头的行保持原样
                        new_mechs.extend([line + '\n' for line in updated_TROE])

                elif Label1 in PLOGs_dict:  # 处理其他 Label 的 PLOG 数据
                    updated_PLOGs = []
                    for plog_line, P_value, A_value, n_value, E_value in PLOGs_dict[Label1]:
                        new_A_value = A_value * Multiplier  # 计算新的 A
                        
                        # 重新构造 PLOG 语句
                        updated_plog = f"PLOG / {float(P_value):<8.2E} {new_A_value:>12.2E} {float(n_value):>8.2f} {float(E_value):>8.0f} /"
                        updated_PLOGs.append(updated_plog)

                    # **确保 PLOG 直接跟在 `!RR` 之后**
                    new_mechs.extend([plog + '\n' for plog in updated_PLOGs])

for i, line in enumerate(new_mechs):
    # 检查该行不以 '!' 开头且包含 '! RR M'
    if not line.startswith('!') and '!RR 12 M' in line:
        # 在 <=> 前后加上 (+M)
        line = re.sub(r'(<=>)', r'(+M) \1', line)
        
        # 以 !RR 为分隔符分割行
        parts = line.split('!RR')

        if len(parts) > 1:
            # 分割第一个部分为新列表
            first_part_elements = re.split(r'\s+', parts[0].strip())
            if len(first_part_elements) >= 3:
                # 将 (+M) 放在新列表的倒数第三个元素前面
                first_part_elements.insert(-3, '(+M)')
            
            # 重新合并修改后的部分
            modified_first_part = ' '.join(first_part_elements)
            # 更新 line 的第一部分
            line = modified_first_part + '   !RR' + parts[1]
        
        # 更新列表中的当前行
        new_mechs[i] = line




with open('completeDetailedMech.MECH', 'w', encoding='utf-8') as file:
    file.writelines(new_mechs)

print("PLOGs done!")


