import numpy as np
import abfile
import warnings
import pickle
warnings.filterwarnings('ignore')

f=open('../VALIDATION/OMmaskTP2.pckl','rb')
REGmaskTP2 = pickle.load(f,encoding="latin1")
f.close()
f=open('../VALIDATION/OMmaskTP5.pckl','rb')
REGmaskTP5 = pickle.load(f,encoding="latin1")
f.close()
# f=open('../VALIDATION/maskTP2extended.pckl','rb')
# REGmaskTP2ext = pickle.load(f,encoding="latin1")
# f.close()
f=open('../VALIDATION/OMmaskTP0.pckl','rb')
REGmaskTP0 = pickle.load(f,encoding="latin1")
f.close()

# get domain dimensions
abgrid = abfile.ABFileGrid("/cluster/work/users/cagyum/TP2a0.10/topo/regional.grid","r")
plonTP2=abgrid.read_field("plon")
platTP2=abgrid.read_field("plat")

abgrid = abfile.ABFileGrid("/cluster/work/users/cagyum/TP5a0.06/topo/regional.grid","r")
plonTP5=abgrid.read_field("plon")
platTP5=abgrid.read_field("plat")

# abgrid = abfile.ABFileGrid("/cluster/work/users/cagyum/TP2extended/regional.grid","r")
# plonTP2ext=abgrid.read_field("plon")
# platTP2ext=abgrid.read_field("plat")

abgrid = abfile.ABFileGrid("/cluster/work/users/cagyum/TP0a1.00/topo/regional.grid","r")
plonTP0=abgrid.read_field("plon")
platTP0=abgrid.read_field("plat")

# get domain depth
jdm,idm=plonTP5.shape
abdepth = abfile.ABFileBathy("/cluster/work/users/cagyum/TP5a0.06/topo/depth_TP5a0.06_05.b", \
           "r",idm=idm,jdm=jdm)
depthTP5=abdepth.read_field("depth")

jdm,idm=plonTP2.shape
abdepth = abfile.ABFileBathy("/cluster/work/users/cagyum/TP2a0.10/topo/regional.depth.b", \
           "r",idm=idm,jdm=jdm)
depthTP2=abdepth.read_field("depth")

# jdm,idm=plonTP2ext.shape
# abdepth = abfile.ABFileBathy("/cluster/work/users/cagyum/TP2extended/regional.depth.b", \
#            "r",idm=idm,jdm=jdm)
# depthTP2ext=abdepth.read_field("depth")

jdm,idm=plonTP0.shape
abdepth = abfile.ABFileBathy("/cluster/work/users/cagyum/TP0a1.00/topo/regional.depth.b", \
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
"subregionTP2ext" : [],
"subregionTP5" : [],
"subregionTP0" : [],
"IITP2" : [],
"JJTP2" : [],
"IITP0" : [],
"JJTP0" : [],
"IITP2ext" : [],
"JJTP2ext" : [],
"IITP5" : [],
"JJTP5" : [] \
}


count = 0
not_the_same = 0
for fileyears in range(1997,2023):
 print(fileyears)
 with open('/nird/projects/NS9481K/BGCDATA/prepobs_bgc/bgc_in_situ_'+str(fileyears)+'0101.txt', 'r') as f:
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

# # TP2 extended domain
#         cooINDEX = abs( platTP2ext-float(row[7]) ) + abs( plonTP2ext-float(row[6]) )
#         JJ,II = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
#         append_TP2ext = "OUTSIDE"
#         INSITU["IITP2ext"].append(int(II))
#         INSITU["JJTP2ext"].append(int(JJ))
#         if depthTP2ext[JJ,II] < 10.0:
#            INSITU["subregionTP2ext"].append(append_TP2ext)
#         else:
#            not_inside = True
#            found = None
#            for key in REGmaskTP2ext:
#              if found is None:
#               if REGmaskTP2ext[key][JJ,II] == 1: 
#                  found = key
#                  append_TP2ext = key
#            INSITU["subregionTP2ext"].append(append_TP2ext)

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
           print("not same region", append_TP5, append_TP2, row[6], row[7], str(not_the_same)+" of "+str(count) )


f = open("INSITU.2025.OMmask.pckl","wb")
pickle.dump(INSITU,f)
f.close()
'''
# example climatolgy 0 - 10 meters
fig=plt.figure(figsize=(12,9))
pcount=1
for t in range(5): 
 if t == 0: tracer = "ntra"; titletrc="nitrate (mmol m-3)"; ymin=0; ymax=15
 if t == 1: tracer = "slca"; titletrc="silicate (mmol m-3)"; ymin=0; ymax=7.5
 if t == 2: tracer = "phos"; titletrc="phosphate (mmol m-3)"; ymin=0; ymax=1
 if t == 3: tracer = "cphl"; titletrc="chl-a (mgChl m-3)"; ymin=0; ymax=4
 if t == 4: tracer = "doxy"; titletrc="oxygen (mmol m-3)"; ymin=240; ymax=380
 for drange in range(4):
#   plt.subplot(5, 4, drange + t*5 +1)
   plt.subplot(5, 4, pcount)
   if drange == 0: depmin=0.; depmax=10.
   if drange == 1: depmin=10.; depmax=30.
   if drange == 2: depmin=30.; depmax=60.
   if drange == 3: depmin=60.; depmax=110.

   for reg in ["NorwegianN","NorwegianS","Barents","Greenland"]:
    clim = np.zeros((12))
    count = 0
    for m in range(1,13):

       index = [index for index, (element1, element2, element3, element4) \
               in enumerate(zip(INSITU["year"],INSITU["subregionTP5"], INSITU["depth"], INSITU["month"])) \
                    if 1990 <= element1 <= 2023 and element2 == reg and depmax > element3 >= depmin  and element4 == m ]

       print(tracer,reg,m,len(index),"depth range :"+str(depmin)+" - "+str(depmax) )
       suffix = " "+str(int(depmin))+"m - "+str(int(depmax))+"m"
       if len(index)>=50:
          clim[count] = np.nanmean(np.array(INSITU[tracer])[index])
       else:
          clim[count] = -999.
       count = count + 1
    clim = np.ma.masked_where(clim<-900,clim)
    plt.plot(np.arange(12)+1,clim,label=reg)
    plt.title(titletrc+suffix,fontsize=10)

   plt.ylim(ymin,ymax)
   plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12])
   if pcount == 16: plt.legend()
   plt.grid(color='grey', linestyle=':')
   pcount = pcount + 1

plt.tight_layout()
'''
