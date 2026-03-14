@echo off
call "C:\Program Files\ANSYS Inc\v193\reaction\chemkinpro.win64\bin\run_chemkinpro_env_setup.bat"
SET CHEMKIN_MODE=Pro
chem -i completeDetailedMech.MECH -d completeDetailedMech.THERM -o chem.out
tran < tran.dat > tran.out
exit