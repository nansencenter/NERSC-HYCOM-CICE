Convert restart files from TOPAZ5 to TOPAZ6 grid 
#
echo "Convert TP5 restart  file to TP6  restart file"
|source file     | purpose|
|-------- | -------------|
|main_convert_TP5restart_2_TP6restart.sh : main script  convert TP5 restart  file to TP6  restart file
| " "
| "Example:"
| "../bin/Initialize/main_convert_TP5restart_2_TP6restart.sh TP5restart.2024_183_00_0000_mem000.a"
| " "
# --- input  TP5  restart  file:
# --- output TP6  restart  file
|---------------------
| "../bin/Initialize/main_convert_TP5restart_2_TP6restart.sh TP5restart.2024_183_00_0000_mem000.a"
|---------------------
echo "Convert TP5 restart  file to TP6  restart file"
|-------------------------------------------------------------------------
# It calls:
|source file     | purpose|
|-------- | -------------|
| 1) topaz_restart2archv_TP5_CML.sh |    to extract TP5 archv file from TP5 restart file
| 2) topaz_TP5archv_2_TP6archv_CML.sh |  to interpolate TP5 archv file to TP6 archv file
| 3) topaz_archv2restart_TP6_CML.sh   |  to create TP6 restart file from TP6 archv file 
| 4) convert_cice_restart.sh          | to interpolate CICE TP5 restart file to TP6 grid
|    - this a wraper which uses: cice_convert_restart.py
| New restart files are created under |TP6a0.03/initialize/New_Restart directory
| "-------------------------------------------------------"
