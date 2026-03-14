import pandas as pd
from rdkit import Chem
import os

# 处理文件路径
#PATH = "C:/Users/22230357/source/repos/Mamox-2/Mamox 2 input example/outputBatchLumping/"
input_file = "speciesList_clean.csv"
output_file = "speciesList_SMILES.csv"
jump_file = "jump.txt"

df = pd.read_csv(input_file)  # 不再 skiprows 和 header=None


required_columns = ["Mech. name", "InChI"]
for col in required_columns:
    if col not in df.columns:
        raise ValueError(f"CSV 文件缺少必要的列: {col}")

name_column = "Mech. name"
inchi_column = "InChI"

for _, row in df.iterrows():
    mech_name = row[name_column]
    inchi = row[inchi_column]
    print(mech_name, inchi)




processed_inchis = set()
processed_data = []
jumped_entries = []

def inchi_to_smiles(inchi):
    try:
        if pd.notna(inchi):  # 确保不是 NaN
            mol = Chem.MolFromInchi(inchi)
            if mol:  # 确保 RDKit 解析成功
                return Chem.MolToSmiles(mol)
            else:
                print(f"警告: 无法解析 InChI -> {inchi}")
    except Exception as e:
        print(f"错误: 解析失败 {inchi} -> {e}")  # 记录错误
    return None  # 


for index, row in df.iterrows():
    mech_name = row[name_column]
    inchi = row[inchi_column]

    
    if inchi in processed_inchis:
        jumped_entries.append(mech_name)
        continue  


    smiles = inchi_to_smiles(inchi)
    if smiles:  # 只有成功转换的才加入
        processed_inchis.add(inchi)  # 记录已处理的 InChI
        processed_data.append([mech_name, smiles])


output_df = pd.DataFrame(processed_data, columns=["MECH_NAME", "SMILES"])


output_df.to_csv(output_file, index=False)


if jumped_entries:
    with open(jump_file, "w", encoding="utf-8") as f:
        f.write("\n".join(jumped_entries))
    
