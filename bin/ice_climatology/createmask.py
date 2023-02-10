import sys 
from os.path import exists
import abfile.abfile as abf
import numpy as np # import necessary modules

Na=len(sys.argv)
if Na<3:
   print ("createmask.py <inputpath> <outfile>")
   sys.exit()
#print("Argument list:",str(sys.argv))
inputpath=sys.argv[1]
outfile=sys.argv[2]
if len(inputpath)<2 or len(outfile)<2:
   print ("Are you ensure for inputs: ",inputpath," and ", outfile,"?")
   sys.exit()

print ("Input model grid file path is ", inputpath)
print ("Output mask file is ", outfile)

# Topography directory and check its accessibility
grdname = inputpath + '/regional.grid.a'
file_exists=exists(grdname)
if file_exists:
   print(f'The model directory pass a check')
else:
   print(f'The input model directory without reginal.grid.a')
   sys.exit()

grdfile = abf.ABFileGrid(grdname,"r")

# read in the necessary grid variables
plon=np.transpose(grdfile.read_field("plon"))
plat=np.transpose(grdfile.read_field("plat"))
qlon=np.transpose(grdfile.read_field("qlon"))
qlat=np.transpose(grdfile.read_field("qlat"))
ulon=np.transpose(grdfile.read_field("ulon"))
ulat=np.transpose(grdfile.read_field("ulat"))
vlon=np.transpose(grdfile.read_field("vlon"))
vlat=np.transpose(grdfile.read_field("vlat"))
scpx=np.transpose(grdfile.read_field("scpx")) 
scpy=np.transpose(grdfile.read_field("scpy"))
scux=np.transpose(grdfile.read_field("scux"))
scuy=np.transpose(grdfile.read_field("scuy"))
scvx=np.transpose(grdfile.read_field("scvx"))
scvy=np.transpose(grdfile.read_field("scvy"))

#write the grid information to the gridinfofile

fileout = open(outfile,"w")

idm,jdm=np.shape(plon)
print(idm,jdm,idm*jdm,(idm-1)*(jdm-1))
fileout.write('gridtype = curvilinear\n');
fileout.write('gridsize = ' + str((idm)*(jdm)) + '\n');
fileout.write('xsize    = ' + str(np.size(plon,0)) + '\n');
fileout.write('ysize    = ' + str(np.size(plon,1)) + '\n');

fileout.write('xvals    = ');
for j in range(np.size(plon,1)):
    for i in range(np.size(plon,0)):
        fval = "{:.3f}".format(plon[i,j])
        fileout.write(fval + " ") 
    fileout.write('\n');   
fileout.write('\n'); 

fileout.write('yvals    = ');
for j in range(np.size(plon,1)):
    for i in range(np.size(plon,0)):
        fval = "{:.3f}".format(plat[i,j])
        fileout.write(fval + " ") 
    fileout.write('\n'); 
fileout.write('\n');

print("The mask file is created and dumped into ", outfile)
