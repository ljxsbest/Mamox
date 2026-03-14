#####################################################
'''
1. Convert transport data of RMG outputs to CHEMKIN format

RMG format: 
(TransportData(shapeIndex=X1, epsilon=(X2,'J/mol'), sigma=(X3,'angstroms'), 
dipoleMoment=(X4,'C*m'), polarizability=(X5,'angstroms^3'), 
rotrelaxcollnum=X6, comment="""Epsilon & sigma estimated with Tc=796.41 K, Pc=48.23 bar (from Joback method)"""), None, None)

CHEMKIN format:
Species Y1 Y2 Y3 Y4 Y5 Y6
Y1:        : Shape Index (0: single-atom; 1: linear; 2: non-linear molecule)
Y2 (K)     : Lennard-Jones Potential Well Depth (epsilon)
Y3 (angs)  : Lennard-Jones Collision Diameter (sigma)
Y4 (debye) : Dipole moment
Y5 (angs^3): Polarizability
Y6         : Rotational Relaxation Collision Number at 298 K
'''
#####################################################

import re

# 输入输出文件名
rmg = "out"
outfile = "trans.dat"

rows = []
with open(rmg, "r", encoding="utf-8") as f:
    lines = f.readlines()

species = None
for line in lines:
    line = line.strip()
    if line.startswith("Species:"):
        species = line.split()[1]   # 取第二个元素作为名字
    if "TransportData" in line:
        shapeIndex = line.split("shapeIndex=")[1].split(",")[0]
        epsilon = float(line.split("epsilon=(")[1].split(",")[0])
        sigma = line.split("sigma=(")[1].split(",")[0]
        dipoleMoment = line.split("dipoleMoment=(")[1].split(",")[0]
        polarizability = line.split("polarizability=(")[1].split(",")[0]
        rotrelaxcollnum = line.split("rotrelaxcollnum=")[1].split(",")[0].replace(")","")
        
        rows.append([species, shapeIndex, epsilon, sigma, dipoleMoment, polarizability, rotrelaxcollnum])

with open(outfile, "w", encoding="utf-8") as f:
    f.write("!Species    shapeID   eps(K)   sigma(angs)  dipoleMoment(debye)  polarizability(angs^3)  rotrelaxcollnum\n")
    for row in rows:
        species, shapeIndex, epsilon, sigma, dipoleMoment, polarizability, rotrelaxcollnum = row
        f.write(f"{species:<22}{shapeIndex:>1}{(float(epsilon)/8.314):12.2f}{float(sigma):10.2f}{float(dipoleMoment):10.2f}{float(polarizability):10.2f}{float(rotrelaxcollnum):10.2f}\n")

print(f"Complete {outfile}")
