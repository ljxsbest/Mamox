# Replace the species names of RMG thermo outputs and sort the information
# Jiaxin Liu, j.liu8@universityofgalway.ie
# 1. Get the real names and temporary names in SPECIES block;
# 2. Replace the temporary names in the THERM block and keep the name length same;
#   2.1. Change the location of species;
#   2.2. Append the current and the following 3 lines;

filename = "chem_annotated.inp"

TEPO = []
NAME = []
inside_block = False
with open(filename, 'r') as f:
    for line in f:
        line = line.strip()  # 去除行首尾的空白字符
        if line.startswith("SPECIES"): 
            inside_block = True
            continue  
        elif line.startswith("END") and inside_block:  
            inside_block = False
            break
        elif inside_block:  # 只处理区块内的内容
            parts = line.split('!')
            if len(parts) > 1:
                TEPO.append(parts[0].split()[0])
                #if len(first) > 1:
                NAME.append(parts[1].split()[0])
#print(NAME)

Marker = False
count = 0
result = []
with open(filename, 'r') as f:
    lines = f.readlines()  
with open('Thermo_clean_RMG.out', 'w') as fw:
	for i in range(len(lines)):
	    line = lines[i].strip()  
	    if line.startswith("THERM"): 
	        Marker = True
	        continue  
	    elif Marker:  # 只处理区块内的内容
	        if len(line) > 0:
	            if line[-1] == '1':
	            	count += 1
	            	species = line.split()[0]
	            	#print(species)
	            	idx = TEPO.index(species)
	            	n = len(NAME[idx]) - len(species)
	            	#print(n)
	            	n1 = len(NAME[idx].split('(')[1])+1
	            	line = NAME[idx].split('(')[0] + n1*' ' +line[len(species)+n:]
	            	print(line[len(species)+n-n1:])
	            	fw.write(line + '\n')
	            	fw.write(lines[i+1])
	            	fw.write(lines[i+2])
	            	fw.write(lines[i+3])
	        		#result.append(block)  # 添加到结果列表中
print('Number of cleaned species: ' + str(count))
            #if len(parts) > 0:
#
#
            #    first = parts[0].split('!')
            #    if len(first) > 0:
            #        TEPO.append(first[0])
            #    if len(first) > 1:
            #        NAME.append(first[1])