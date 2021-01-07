import numpy
import modeltools.hycom
import abfile
import os,fnmatch
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import pickle

# this code prepares the regression (slope, intercept) map. It's output is already copied to
# hycom input folder. This code is here for future modifications on the domain.


def getvarib(abfile,fieldname,fieldlevel) :
    fld = abfile.read_field(fieldname,fieldlevel)
    return fld

def plthycom(varib,title,figx,figy,cmin,cmax,types) :

  fig=plt.figure(title,figsize=(figx,figy),facecolor='w')
  ax = fig.add_subplot(1,1,1)
  if (types=='diff'):
    cmap = plt.get_cmap('RdYlBu_r')
  else:
    cmap = plt.get_cmap('Spectral_r')
  ax.set_facecolor('xkcd:gray')

  levels = MaxNLocator(nbins=40).tick_values(cmin, cmax)
  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

  pmesh = plt.pcolormesh(varib,cmap=cmap,norm=norm)
  pmesh.set_clim(cmin,cmax)
  cb=plt.colorbar(pmesh)

  return plt.show()

# get free run ssh and montg1 for regression #################
# calculate mean dt and ssh anomaly
# regress anomaly and montg1
day=0
#directory = '/cluster/work/users/annettes/TP5a0.06/expt_01.1/data/'
#ndays = numpy.size(sorted(fnmatch.filter(os.listdir(directory), 'archm.%s*.b' %('201[01123456]') ) ))

directory = '/cluster/work/users/achoth/TP2a0.10/expt_01.0/data/'
ndays = numpy.size(sorted(fnmatch.filter(os.listdir(directory), 'archm.%s*.b' %('201[01234]') ) ))

NX = 380; NY = 400
montg = numpy.zeros((ndays,NX,NY))
ssh   = numpy.zeros((ndays,NX,NY))
a     = numpy.zeros((NX,NY))
b     = numpy.zeros((NX,NY))


for f1 in sorted(fnmatch.filter(os.listdir(directory), 'archm.%s*.b' %('201[01234]') ) ) :
 print(day)
 f = abfile.ABFileArchv('%s%s' %(directory,f1),"r")

 mon = getvarib(f,'montg1',0)
 srf = getvarib(f,'srfhgt',0)

 montg[day,:,:] = mon
 ssh[day,:,:]   = srf


 day = day + 1
 f.close

t=0
ssh1 = numpy.ma.masked_where(ssh>1E10,ssh)
meanssh = numpy.mean(ssh1)
anomaly = ssh1 - meanssh
for II in range(NX):
  for JJ in range(NY):

     print('steps left: %s' %(str(NX*NY-t)) )
     slope, intercept, r_value, p_value, std_err = stats.linregress(anomaly[:,II,JJ],montg[:,II,JJ])
#     meandt[II,JJ] = numpy.mean(ssh1[:,II,JJ])
#     slope, intercept, r_value, p_value, std_err = stats.linregress(ssh1[:,II,JJ]-meandt[II,JJ],montg[:,II,JJ])   
     a[II,JJ] = slope
     b[II,JJ] = intercept 

     t = t + 1

del montg,ssh,mon,srf,ssh1
##############################################################

# get nemo archive files interpolated to hycom grid ##########
# calculate mean dt
nestdir = '/cluster/work/users/achoth/TP2a0.10/nest/010/'
nestdays = numpy.size(sorted(fnmatch.filter(os.listdir(nestdir), 'archv*.a') ))
nest = 0
nestmean = 0.
for f1 in sorted(fnmatch.filter(os.listdir(nestdir), 'archv*.a') ) :
   print(nestdays-nest)
   f = abfile.ABFileArchv('%s%s' %(nestdir,f1),"r")

   nestsrf = getvarib(f,'srfhgt',0)
   nestmean = nestmean + numpy.mean(nestsrf) / nestdays

   nest = nest + 1
   f.close
##############################################################
#f = open("/cluster/work/users/cagyum/model_output/montg_regress.pckl","wb")
f = open("/cluster/work/users/achoth/TP2a0.10/TP2_montg_regress.pckl","wb")
pickle.dump([a,b,nestmean],f)
f.close()

# test
#file = '%s%s' %(directory,'archm.2013_250_12.b')
#file = abfile.ABFileArchv(file,"r")

#mont_test = getvarib(file,'montg1',0)
#ssh_test  = getvarib(file,'srfhgt',0)
#anom_test = ssh_test - meanssh

#mont_regr = anom_test * a + b
#cmin = numpy.min([mont_test.min(),mont_regr.min()])
#cmax = numpy.max([mont_test.max(),mont_regr.max()])

#plthycom(mont_test,'from model',7,6,cmin,cmax,'')
#plthycom(mont_regr,'calculated',7,6,cmin,cmax,'')

#cmin = numpy.min([mont_test[0:20,300:720].min(),mont_regr[0:20,300:720].min()])
#cmax = numpy.max([mont_test[0:20,300:720].max(),mont_regr[0:20,300:720].max()])

#plthycom(mont_test[0:20,300:720],'south from model',7,6,cmin,cmax,'')
#plthycom(mont_regr[0:20,300:720],'south calculated',7,6,cmin,cmax,'')

#cmin = numpy.min([mont_test[740:759,0:210].min(),mont_regr[740:759,0:210].min()])
#cmax = numpy.max([mont_test[740:759,0:210].max(),mont_regr[740:759,0:210].max()])

#plthycom(mont_test[740:759,0:210],'east from model',7,6,cmin,cmax,'')
#plthycom(mont_regr[740:759,0:210],'east calculated',7,6,cmin,cmax,'')
