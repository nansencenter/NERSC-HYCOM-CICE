FES2004 package adapded to for the hycom model with tides at the boundary. 
Routines created by Inti - some setup changes done by Knut 

1--- go to  your model directory
     need regional.grid and regional.depth files
     run fes2mod
     get the FES_obs_elev.dat
     (KAL: in HYCOM 2.2 version the routines are setup under main model dir)

2--- copy FES_obc_elev.dat to /work/yourlogin/yourwmodeldir/Data/ 
     KAL: handled transparently in HYCOM 2.2 (no need to copy)

3--- in infile.in: put the tide on with the flag FES : T FES F F

4---If need more info about FES:
             read ./FES_Readme.txt 
             ./nersc_src/Inti_notes.doc describes the list of excutables available in
             src directory

###############################################################################
#Knut - these were the original instructions by Inti

1--- go to ./src
     need depthsYYYxYYY.uf and newpos.uf
     run fes2mod
     get the FES_obs_elev.dat

2--- copy FES_obc_elev.dat to /work/yourlogin/yourwmodeldir/Data/

3--- in infile.in: put the tide on with the flag FES : T FES F F


4---If need more info about FES:
             read ./FES_Readme.txt 
             ./src/Inti_notes.doc describes the list of excutables available in
             src directory
