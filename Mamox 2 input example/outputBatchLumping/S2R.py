# S2R.py: Smiles to rmg input files
# Jiaxin Liu: j.liu8@universityofgalway.ie
# 2024-11-28
#
# 1. Get Mech_names and Smiles first;
# 2. Convert into the RMG format.
import os
import pandas as pd

pwd = os.path.dirname(os.path.abspath(__file__))
os.chdir(pwd)
print(os.getcwd())

path = pwd+'\\'+"speciesList_SMILES.csv"
df = pd.read_csv(path, usecols=["MECH_NAME", "SMILES"])

Name   = df["MECH_NAME"].tolist()  
Smiles = df["SMILES"].tolist()   

#print(Name)
#print(Smiles)
with open('input_RMG.py', 'w') as f:
	f.write("database(thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0'])")
	for i, name in enumerate(Name):
		f.write('\n')
		f.write("species(label='"+str(name)+"',structure=SMILES('"+str(Smiles[i])+"'))")

#database(thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0'])
#species(label='ethane',structure=SMILES('CC'),)