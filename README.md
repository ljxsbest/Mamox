# MAMOX 2

Authors:
- Jiaxin Liu (J.liu8@universityofgalway.ie / liujiaxin9702@gmail.com)
- Sirio Brunalti (sirio.brunialti@kaust.edu.sa)
- Enrico Garavaglia
- Tiziano Faravelli (tiziano.faravelli@polimi.it)
- Eliseo Ranzi
  
For any inquire please contact Sirio Brunialti - sirio.brunialti@kaust.edu.sa

HOW TO CITE  
If a publication is produced using this software the following references must be reported
- J. Liu, et al., An ignition delay time and automated chemical kinetic modeling study of three heptene isomers, Under review.
- S. Brunialti, X. Zhang, T. Faravelli, A. Frassoldati, S.M. Sarathy, Automatically generated detailed and lumped reaction mechanisms for low- and high-temperature oxidation of alkanes, Proceedings of the Combustion Institute 39 (2023) 335-344.
- E. Ranzi, A. Frassoldati, S. Granata, T. Faravelli, Wide-Range Kinetic Modeling Study of the Pyrolysis, Partial Oxidation, and Combustion of Heavy n-Alkanes, Ind. Eng. Chem. Res. 44 (2005) 5170-5183
- E. Ranzi, T. Faravelli, P. Gaffuri, A. Sogaro, Low-Temperature Combustion - Automatic-Generation of Primary Oxidation Reactions and Lumping Procedures, Combustion and Flame 102 (1995) 179-192.

WARNING  
This version of MAMOX2 and its source code are intended for educational purposes. Use for commercial purposes is not permitted. For any commercial issue please contact	Sirio Brunialti (sirio.brunialti@kaust.edu.sa) and Tiziano Faravelli (tiziano.faravelli@polimi.it).

LIMITED WARRANTY  
The Software and related documentation are provided “AS IS” and without	any warranty of any kind and Seller EXPRESSLY DISCLAIMS ALL WARRANTIES,	EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES	OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

## HOW TO USE
From source: Compile the source file An example is provided in the "Mamox 2 input example.zip" file. Unzip the file in the Documents folder of your computer. To test the example compile Mamox2, provide "libnlopt.dll" from NLopt library (https://nlopt.readthedocs.io/en/latest/) and "libinchi.dll" from InChI library (https://www.inchi-trust.org/downloads/), and run the following command:

"Mamox 2.exe" "C:\users%username%\Documents\Mamox 2 input example\config.txt"


## INPUT FILE
The input file contains all the information to run the kinetic mechanism generation job.  
The comment character is "!".
The information is provided through keywords. A keyword has to be followed by "=" and then by the information to provide, the line has to be closed with the ";" character.
The keywords are the following:  
`MOL_FOLDER` absolute path to the folder containing the molecules file. The molecules files need to have the .hyd extension.  
`BASE_MECH` absolute path to the file containing the base mechanism.  
`BASE_THERMO` absolute path to the file containing the reference thermochemical properties.  
`BASE_NAMES` absolute path to the reference glossary.  
`THERMO_GROUP_FILE` absolute path to the file containing the Benson group additivy values.  
`THERMO_HBI_FILE` absolute path to the file containing the HBI values.  
`RATE_RULES_FILE` absolute path to the file containing the rate rules.  
`FUELS` the list of fuels the full kinetic submechanism is generated for. The names have to be the same of the file in the MOL_FOLDER without the .hyd extension. More molecules can be passed and they have to be separated by the "," character.  
`OUT_FOLDER` absolute path to the output folder.  
`GENERATE_LUMPED` set value to `TRUE` to generate the lumped mechanisms, otherwise set value to `FALSE`.  
`SIM_TYPE` set to `CSTR` to perform CSTR simulations, set `BATCH` to perform batch simulations.  
`LUMP_WITH+CORE_MECH` set to `TRUE` to merge the base mechanism to the fuel submechanism when running the simulations for the lumping procedure, otherwise set to `FALSE`.  
`LUMPING_T` the set of temperatures in K at which to perform the simulations to get the reference concentration for the lumping methodology. Values are separated by the ',' character.  
`LUMPING_P` the set of pressures in atm at which to perform the simulations to get the reference concentration for the lumping methodology. Values are separated by the ',' character.  
`EQ_RATIO` the equivalence ratio of the fuel/air mixture for the simulations.  
`TA` the residence time in s for the CSTR simulations.  
`NUM_CORES` the number of cores to use in parallel to run the simulations for the lumping methodology.  


