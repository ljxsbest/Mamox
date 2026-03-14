@echo off
call "C:\Program Files\ANSYS Inc\v193\reaction\chemkinpro.win64\bin\run_chemkinpro_env_setup.bat"
SET CHEMKIN_MODE=Pro
chem -i kinetic.kin -d thermo.ter -o chem.out
tran < tran.dat > tran.out
exit