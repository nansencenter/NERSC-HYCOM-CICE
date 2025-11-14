import numpy as np
import abfile
import warnings
import pickle
import argparse
warnings.filterwarnings('ignore')

# CAGLAR, 11.Nov.2025
# You only need to run this once when the in situ data files are updated, or region definitions have changed.
# The script ignores depths below 10.0 meters.
# usage: python convert_txt_pckl_insitu_preobs_OMmask.py ~/NERSC-HYCOM-CICE/ 1997 2023

def main(folder,year1,year2):

   f=open('/cluster/projects/nn9481k/BGC.Validation/OMmaskTP2.pckl','rb')
   REGmaskTP2 = pickle.load(f,encoding="latin1")
   f.close()
   f=open('/cluster/projects/nn9481k/BGC.Validation/OMmaskTP5.pckl','rb')
   REGmaskTP5 = pickle.load(f,encoding="latin1")
   f.close()
   f=open('/cluster/projects/nn9481k/BGC.Validation/OMmaskTP0.pckl','rb')
   REGmaskTP0 = pickle.load(f,encoding="latin1")
   f.close()

   # get domain dimensions
   abgrid = abfile.ABFileGrid(folder+"TP2a0.10/topo/regional.grid","r")
   plonTP2=abgrid.read_field("plon")
   platTP2=abgrid.read_field("plat")

   abgrid = abfile.ABFileGrid(folder+"TP5a0.06/topo/regional.grid","r")
   plonTP5=abgrid.read_field("plon")
   platTP5=abgrid.read_field("plat")

   abgrid = abfile.ABFileGrid(folder+"TP0a1.00/topo/regional.grid","r")
   plonTP0=abgrid.read_field("plon")
   platTP0=abgrid.read_field("plat")

   # get domain depth

   jdm,idm=plonTP5.shape
   abdepth = abfile.ABFileBathy(folder+"/TP5a0.06/topo/regional.depth.b", \
           "r",idm=idm,jdm=jdm)
   depthTP5=abdepth.read_field("depth")

   jdm,idm=plonTP2.shape
   abdepth = abfile.ABFileBathy(folder+"/TP2a0.10/topo/regional.depth.b", \
           "r",idm=idm,jdm=jdm)
   depthTP2=abdepth.read_field("depth")

   jdm,idm=plonTP0.shape
   abdepth = abfile.ABFileBathy(folder+"/TP0a1.00/topo/regional.depth.b", \
           "r",idm=idm,jdm=jdm)
   depthTP0=abdepth.read_field("depth")

   # get data
   INSITU = {\
   "provider" : [],
   "expocode" : [],
   "year" : [],
   "month" : [],
   "day" : [],
   "hour" : [],
   "longitude" : [],
   "latitude" : [],
   "depth" : [],
   "temp" : [],
   "psal" : [],
   "doxy" : [],
   "phos" : [],
   "ntra" : [],
   "slca" : [],
   "cphl" : [],
   "subregionTP2" : [],
   "subregionTP5" : [],
   "subregionTP0" : [],
   "IITP2" : [],
   "JJTP2" : [],
   "IITP0" : [],
   "JJTP0" : [],
   "IITP5" : [],
   "JJTP5" : [] \
   }


   count = 0
   not_the_same = 0
   for fileyears in range(year1,year2+1):
    print(fileyears)
    with open('/cluster/projects/nn9481k/BGCDATA/prepobs_bgc/bgc_in_situ_'+str(fileyears)+'0101.txt', 'r') as f:
     # Skip the first two header lines
     next(f)
     next(f)
 
      # Loop through each line in the file
     for line in f:
        count = count + 1
        if np.mod(count,5000) == 0 : print(fileyears,count)
        # Split the line into a list
        row = line.split()

       # cooINDEX = abs( platTP2-float(row[7]) ) + abs( plonTP2-float(row[6]) )
       # JJ,II = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
       # if depthm[JJ,II] >= 10.0:
 
           # Add each column to the appropriate list
        INSITU["provider"].append(row[0])
        INSITU["expocode"].append(row[1])
        INSITU["year"].append(int(row[2]))
        INSITU["month"].append(int(row[3]))
        INSITU["day"].append(int(row[4]))
        INSITU["hour"].append(float(row[5]))
        INSITU["longitude"].append(float(row[6]))
        INSITU["latitude"].append(float(row[7]))
        INSITU["depth"].append(-float(row[8]))
        INSITU["temp"].append(float(row[9]))
        INSITU["psal"].append(float(row[10]))
        INSITU["doxy"].append(float(row[11]))
        INSITU["phos"].append(float(row[12]))
        INSITU["ntra"].append(float(row[13]))
        INSITU["slca"].append(float(row[14]))
        INSITU["cphl"].append(float(row[15]))

        # TP0 domain
        cooINDEX = abs( platTP0-float(row[7]) ) + abs( plonTP0-float(row[6]) )
        JJ,II = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
        append_TP0 = "OUTSIDE"
        INSITU["IITP0"].append(int(II))
        INSITU["JJTP0"].append(int(JJ))
        if depthTP0[JJ,II] < 10.0:
           INSITU["subregionTP0"].append(append_TP0)
        else:
           not_inside = True
           found = None
           for key in REGmaskTP0:
             if found is None:
              if REGmaskTP0[key][JJ,II] == 1:
                 #subregionTP0.append(key)
                 found = key
                 append_TP0 = key
#                 not_inside = False
#           if not_inside:
#              subregionTP0.append("OUTSIDE")
#           else:
#              subregionTP0.append(found)
           INSITU["subregionTP0"].append(append_TP0)

        # TP2 domain
        cooINDEX = abs( platTP2-float(row[7]) ) + abs( plonTP2-float(row[6]) )
        JJ,II = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
        append_TP2 = "OUTSIDE"
        INSITU["IITP2"].append(int(II))
        INSITU["JJTP2"].append(int(JJ))
        if depthTP2[JJ,II] < 10.0:
           INSITU["subregionTP2"].append(append_TP2)
        else:
           not_inside = True
           found = None
           for key in REGmaskTP2:
             if found is None:
              if REGmaskTP2[key][JJ,II] == 1: 
                 #subregionTP2.append(key)
                 found = key
                 append_TP2 = key
#                 not_inside = False
#           if not_inside:
#              subregionTP2.append("OUTSIDE")
#           else:
#              subregionTP2.append(found)
           INSITU["subregionTP2"].append(append_TP2)

        # TP5 domain
        cooINDEX = abs( platTP5-float(row[7]) ) + abs( plonTP5-float(row[6]) )
        JJ,II = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)                  
        append_TP5 = "OUTSIDE"
        INSITU["IITP5"].append(int(II))
        INSITU["JJTP5"].append(int(JJ))
        if depthTP5[JJ,II] < 10.0:
           INSITU["subregionTP5"].append(append_TP5)
        else:
           not_inside = True
           found=None
           for key in REGmaskTP5:
             if found is None:               
               if REGmaskTP5[key][JJ,II] == 1: 
                  #subregionTP5.append(key)
                  found = key
                  append_TP5 = key
                  #not_inside = False
#           if not_inside: 
#              subregionTP5.append("OUTSIDE")
#           else:
#              subregionTP5.append(found)
           INSITU["subregionTP5"].append(append_TP5)

        if append_TP5 != append_TP2 : 
           not_the_same = not_the_same + 1
           print("not same region", "TP5:" ,append_TP5,"TP2:" ,append_TP2,"lon:" ,row[6],"lat:", row[7],"total mismatch of", str(not_the_same)+" out of "+str(count) )


   f = open("./INSITU.OMmask.pckl","wb")
   pickle.dump(INSITU,f)
   f.close()

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('folder',  type=str, help='a folder that likely has regions e.g. ~/NERSC-HYCOM-CICE/')
    parser.add_argument('year1',  type=int)
    parser.add_argument('year2',  type=int)
    args = parser.parse_args()

    main(args.folder,args.year1,args.year2)
