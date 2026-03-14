import os
from tqdm import tqdm
import sys
import re

#input info
Mech_name = 'kinetic'

#return the number of atoms
def atom_find(species, mech):
    atom_counts = {'C': 0.0, 'H': 0.0, 'O': 0.0, 'N': 0.0}
    atom_fields = []

    # Open the .THERM file safely
    with open(f"{mech}.THERM", 'r') as f:
        for line in f:
            if not line.strip():
                continue
            if line.split()[0] == species:
                # Extract the fixed-width fields (atoms are located at specific positions)
                atom_fields = [line[24:29], line[29:34], line[34:39], line[39:44]]
                break

    # Parse the atom information
    for field in atom_fields:
        parts = field.strip().split()
        if len(parts) != 2:
            continue
        atom, count = parts
        if atom in atom_counts:
            atom_counts[atom] = float(count)

    return atom_counts['C'], atom_counts['H'], atom_counts['O'], atom_counts['N']

#########################################################################################
#                READING SPECIES FROM REACTION BLOCK AND SORT IT                        #
#########################################################################################
SPECIES_LIST = []
# Open the mechanism file safely
with open(f"{Mech_name}.MECH", 'r') as f:
    for line in f:
        # Only consider lines that contain reactions (before '!' comment)
        reaction_part = line.split('!')[0]
        if '=' not in reaction_part:
            continue
            
        # Extract Arrhenius A-coefficient (assumed to be third-to-last item)
        try:
            arrhenius_A = reaction_part.strip().split()[-3]
        except IndexError:
            continue  # skip if line is too short

        # Clean the reaction line before the A-coefficient
        cleaned_line = line.split(arrhenius_A)[0]
        cleaned_line = cleaned_line.replace(' ', '').replace('>', '').replace('<', '').replace('\t', '')
        cleaned_line = re.sub(r'\(\+[^)]*\)', '', cleaned_line)  # remove (+M) or similar

        # Split into reactants and products
        reactants, products = cleaned_line.split('=')

        for sp in reactants.split('+') + products.split('+'):
            if sp and (sp[0].isdigit() or sp == 'M'):
                continue  # skip numeric prefixes and 'M'
            SPECIES_LIST.append(sp)

# Add inert species
SPECIES_LIST.extend(['AR', 'HE'])

#remove the duplicate elements
SPECIES_LIST = list(set(SPECIES_LIST))

#obtain the C/H/O/N atom number for each species
C_list = []
H_list = []
O_list = []
N_list = []

print('Reading the atom numbers...')
for i,species in enumerate(tqdm(SPECIES_LIST)):
    C_list.append(atom_find(species,Mech_name)[0])
    H_list.append(atom_find(species,Mech_name)[1])
    O_list.append(atom_find(species,Mech_name)[2])
    N_list.append(atom_find(species,Mech_name)[3])
    
# Combine all data into a list of tuples, and sort it
combined = list(zip(SPECIES_LIST, C_list, H_list, O_list, N_list))
sorted_combined = sorted(combined, key=lambda x: (x[1], x[2], x[3], x[4], x[0]))

# Extract the sorted species
sorted_species = [item[0] for item in sorted_combined]


#########################################################################################
#                                 WRITING .MECH FILE                                    #
#########################################################################################
print('Writing MECH file...',end='')
# ELEMENTS BLOCK
f_clean = open(f"{Mech_name}_clean.MECH", 'w')
f_clean.write('ELEMENTS\n')
f_clean.write('C H O N AR HE\n')
f_clean.write('END\n')
f_clean.write('\n')

# SPECIES BLOCK
f_clean.write('SPECIES\n')
for i,SPECIES in enumerate(sorted_species):
    f_clean.write('{0:30s}'.format(SPECIES))
    if (i+1)%3 == 0: 
        f_clean.write('\n')
f_clean.write('END\n')
f_clean.write('\n')

# REACTION BLOCK
with open(f"{Mech_name}.MECH", 'r') as f:
    reaction_found = False
    for line in f:
        if line.strip().upper().startswith('REACTIONS'):
            reaction_found = True
        
        if reaction_found == True:
            if '=' in line.split('!')[0]:    
            
                if '!' in line:
                    comment = '!'+line.split('!', 1)[1].replace('\n','').strip()
                else:
                    comment = ''                
                
                A = line.split('!')[0].strip().split()[-3]
                n = float(line.split('!')[0].strip().split()[-2])
                E = float(line.split('!')[0].strip().split()[-1])
                Reaction = line.replace('\t',' ').split(A+' ')[0].replace(' ','')
                A = float(A)
                
                Reaction_clean = Reaction.replace(' ','')
                if '<=>' in Reaction_clean or '=>' in Reaction_clean:
                    Reaction_clean = Reaction_clean
                else:
                    Reaction_clean = Reaction_clean.replace('=', '<=>')
                    
                line = '{0:50s} {1:.8E} {2:.8E} {3:.8E} {4}\n'.format(Reaction_clean,A,n,E,comment)
                
            if line.upper().strip().startswith('PLOG'):
                pressure = float(line.split('/')[1].strip().split()[0])
                A        = float(line.split('/')[1].strip().split()[1])
                n        = float(line.split('/')[1].strip().split()[2])
                E        = float(line.split('/')[1].strip().split()[3])
                line = 'PLOG /                              {0:.8E} {1:.8E} {2:.8E} {3:.8E}/\n'.format(pressure,A,n,E)
            f_clean.write(line)
            
f_clean.close()
print('Done')

#########################################################################################
#                                WRITING .THERM FILE                                    #
#########################################################################################

f_clean = open(f"{Mech_name}_clean.THERM", 'w')
f_clean.write('THERMO ALL\n')
f_clean.write('300.   1000.   5000.\n')

print('Writing THERMO file...')
for i,SPECIES in enumerate(tqdm(sorted_species)):

    f = open(f"{Mech_name}.THERM", 'r')
    while True:
        line = f.readline()
        if line == '':
            break
        if len(line.split('!')[0].split())>5:
            if SPECIES == line.split()[0]:
                f_clean.write(line)
                line = f.readline()
                f_clean.write(line)
                line = f.readline()
                f_clean.write(line)
                line = f.readline()
                f_clean.write(line)
                break
    f.close()

f_clean.write('END')
f_clean.close()
    
