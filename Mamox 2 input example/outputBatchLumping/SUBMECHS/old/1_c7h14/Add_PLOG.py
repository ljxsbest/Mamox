'''
Algorithm:
1. Find all the labels in PLOGs.txt;
2. Get the PLOGs behind the label;
3. Find the label from the mechs;
4. Write the PLOGs to the mech
'''

mechname = 'kinetic_detailed.inp'
RR_label, A_list, n_list, E_list, PLOGs = [], [], [], [], []
with open('PLOGs.txt', 'r') as f:
	lines = f.readlines()
with open(mechname, 'r') as f1:
	mechs = f1.readlines()

i = 0
#while i < len(lines):
for i, line in enumerate(lines):
	line = lines[i].strip()
	#print(line)
	if line.startswith('Label'):
		label = line.split()[-1]
		#print('label ' +label)
		#RR_label.append(label)
		print(lines[i+1].strip())
		while lines[i+1].strip().startswith('PLOG') and (i+1 < len(lines)-1):
			PLOGs.append(lines[i+1].strip())  # print PLOG
			#print(lines[i+1].strip())  # print PLOG
			A_list.append(float(lines[i+1].strip().split('/')[1].split()[1]))
			n_list.append(float(lines[i+1].strip().split('/')[1].split()[2]))
			E_list.append(float(lines[i+1].strip().split('/')[1].split()[3]))
			i += 1

		#print(PLOGs)
		for j, mechline in enumerate(mechs):
			#print('1')
			mechline = mechline.strip()
			if not mechline or mechline.startswith('!'):
				continue
			if 'RR' in mechline:
				print(mechline.strip().split('!RR')[1].split())
				Label1 = mechline.strip().split('!RR')[1].split()[0]
				Multiplier = mechline.strip().split('!RR')[1].split()[2]
				if Label1 == label:
					mechs[j+1:j+1] = [item + '\n' for item in PLOGs]
		PLOGs = []

with open('kinetic_detailed_PLOG.inp', 'w', encoding='utf-8') as file:
	file.writelines(mechs)