import re

mechname = 'test1.txt'
RR_label, A_list, n_list, E_list, PLOGs, extra_lines = [], [], [], [], [], []

with open('PLOGs.txt', 'r') as f:
    lines = f.readlines()
with open(mechname, 'r') as f1:
    mechs = f1.readlines()

i = 0
for i, line in enumerate(lines):
    line = lines[i].strip()
    if line.startswith('Label'):
        label = line.split()[-1]
        PLOGs, extra_lines = [], []
        is_label_12 = label == '12'
        
        while i + 1 < len(lines) and lines[i + 1].strip().startswith('PLOG'):
            PLOGs.append(lines[i + 1].strip() + '\n')
            i += 1
        
        if is_label_12 and i + 3 < len(lines):
            extra_lines = [lines[i + j].strip() + '\n' for j in range(1, 4)]
            i += 3
        
        inserted = False  # 确保内容只插入一次
        for j, mechline in enumerate(mechs):
            mechline = mechline.strip()
            if not mechline or mechline.startswith('!'):
                continue
            
            if '!RR' in mechline:
                parts = mechline.split('!RR')[1].split()
                if len(parts) >= 3:
                    Label1, Multiplier = parts[0], float(parts[2])
                    if Label1 == label and not inserted:
                        modified_PLOGs = []
                        for plog in PLOGs:
                            segments = plog.split('/')
                            if len(segments) > 1:
                                values = re.split(r'[\s\t]+', segments[1].strip())
                                if len(values) >= 4:  # A 是原位置的下一个数据
                                    A_new = float(values[1]) * Multiplier
                                    values[1] = f'{A_new:.6e}'
                                    new_plog = f'{segments[0]}/ ' + ' '.join(values) + ' /\n'
                                    modified_PLOGs.append(new_plog)
                        
                        insert_content = modified_PLOGs + extra_lines  # 先插入 PLOGs，再插入额外内容
                        mechs[j + 1:j + 1] = insert_content  # 确保内容被正确插入到 !RR 所在行的下一行
                        inserted = True  # 标记已插入，防止重复
                        break  # 退出循环，避免二次插入
        
        PLOGs, extra_lines = [], []

with open('test2.txt', 'w', encoding='utf-8') as file:
    file.writelines(mechs)