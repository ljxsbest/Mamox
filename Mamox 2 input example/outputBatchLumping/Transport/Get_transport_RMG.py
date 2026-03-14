import pandas as pd
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.transport import TransportDatabase

# 读取 CSV
species_df = pd.read_csv('speciesList_SMILES.csv')

# 存储输出行
output_lines = []

database = TransportDatabase()
# species     = Species() Do not use this , use 'Species().' directly when needed.
database.load(path='/home/chemlab/Desktop/rmg_env/RMG-database/input/transport', libraries=['GRI-Mech'])
for index, row in species_df.iterrows():
    label = row['MECH_NAME']
    smiles = row['SMILES']
    print(smiles)
    # 构建分子对象
    try:
        sp = Species().from_smiles(smiles)
        print(sp)
    except Exception as e:
        print(f"Error parsing SMILES {smiles}: {e}")
        continue
    
    # 估算输运参数
    try:
        transport_data = database.get_transport_properties(sp)
        print(transport_data)
        
    except Exception as e:
        print(f"Error estimating transport for {label}: {e}")
        continue
    

