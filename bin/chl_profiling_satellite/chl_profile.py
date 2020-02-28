import numpy as np
from netCDF4 import Dataset as NetCDFFile
import pickle
import abfile
import modeltools.hycom
# needed only plotting #
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import ticker
#
import datetime
import warnings
########################
warnings.filterwarnings('ignore')
# TP5 domain related stuff ##
I1 = 760; I2 = 800; dm = 50
bfile=abfile.ABFileBathy("/cluster/work/users/cagyum/TP5a0.06/topo/depth_TP5a0.06_08.a","r",idm=I2,jdm=I1)
depthm=bfile.read_field("depth")
gfile=abfile.ABFileGrid("/cluster/work/users/cagyum/TP5a0.06/topo/regional.grid",'r')
plon = gfile.read_field("plon")
plat = gfile.read_field("plat")
gfile.close()
#############################
#glob = '/cluster/work/users/cagyum/program_files/dataset-oc-glo-bio-multi-l4-chl_interpolated_4km_daily-rt_1580117362903.nc'
glob = '/cluster/work/users/cagyum/OCCCI/daily/chl/2013/ESACCI-OC-L3S-CHLOR_A-MERGED-1D_DAILY_4km_GEO_PML_OCx-20130515-fv4.0.nc'
nc = NetCDFFile(glob)
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
time= nc.variables['time'][:]
#satchl = nc.variables['CHL'][0,:,:]
satchl = nc.variables['chlor_a'][0,:,:]
#glob = '/cluster/work/users/cagyum/program_files/dataset-oc-glo-opt-multi-l4-kd490_interpolated_4km_daily-rt_1580123609065.nc'
#nc = NetCDFFile(glob)
#satkd  = nc.variables['KD490'][0,:,:]

# load restart file #
print 'reading restart file'
resfile = '/cluster/work/users/cagyum/TP5a0.06/expt_02.0/data/restart.2013_135_00_0000.a'
i_abfile = abfile.ABFileRestart(resfile,"r",idm=I2,jdm=I1)
dia3D = np.zeros((dm,I1,I2))
fla3D = np.zeros((dm,I1,I2))
diachl3D = np.zeros((dm,I1,I2))
flachl3D = np.zeros((dm,I1,I2))
dp3D     = np.zeros((dm,I1,I2))
depth3D  = np.zeros((dm,I1,I2))

dianew    = np.zeros((dm,I1,I2))
diachlnew = np.zeros((dm,I1,I2))
flanew    = np.zeros((dm,I1,I2))
flachlnew = np.zeros((dm,I1,I2))
for k in range(dm):
    dia3D[k,:,:]    = i_abfile.read_field('ECO_dia',k+1)
    fla3D[k,:,:]    = i_abfile.read_field('ECO_fla',k+1)
    diachl3D[k,:,:] = i_abfile.read_field('ECO_diac',k+1)
    flachl3D[k,:,:] = i_abfile.read_field('ECO_flac',k+1)
    dp3D[k,:,:]     = i_abfile.read_field('dp',k+1)/modeltools.hycom.onem
    if k == 0:
       depth3D[k,:,:] = dp3D[k,:,:]/2.
    else:
       depth3D[k,:,:] = depth3D[k-1,:,:] + dp3D[k,:,:]/2.

    dianew[k,:,:] = dia3D[k,:,:]
    flanew[k,:,:] = fla3D[k,:,:]
    diachlnew[k,:,:] = diachl3D[k,:,:]
    flachlnew[k,:,:] = flachl3D[k,:,:]
#####################


print 'mapping satellite chl to TP5 domain'
mapped = np.zeros((I1,I2))
for j in range(I1):
  for i in range(I2):
     #print(i,j)
     II = np.abs(plon[j,i]-lon).argmin()
     JJ = np.abs(plat[j,i]-lat).argmin()

     if np.ma.is_masked(depthm[j,i]):
        mapped[j,i] = -999.
     elif np.ma.is_masked(satchl[JJ,II]) :
        mapped[j,i] = -999.
     else:
        mapped[j,i] = satchl[JJ,II]

mapped = np.ma.masked_where(mapped == -999.,mapped)
mapped = np.ma.masked_where(depthm<100.,mapped)

'''
m = Basemap(width=8500000,height=8500000,
resolution='c',projection='stere',\
lat_ts=80,lat_0=80,lon_0=-40.,round='True')
x,y=m(plon,plat)
satlon,satlat = np.meshgrid(lon,lat)
x1,y1 = m(satlon,satlat)

cmin = 0.; cmax = 5.
warnings.filterwarnings('ignore')
cmap = plt.get_cmap('Spectral_r')
figure=plt.figure(figsize=(12,6))
ax=figure.add_subplot(211)
ax.set_position([0.0,0.05,0.5,0.95])

levels = MaxNLocator(nbins=15).tick_values(cmin, cmax)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
pmesh = m.pcolormesh(x,y,mapped,cmap=cmap,norm=norm)

m.drawcoastlines(linewidth=0.25)
m.fillcontinents(color='lightgrey')
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='whitesmoke')

ax2=figure.add_subplot(212)
ax2.set_position([0.5,0.05,0.5,0.95])

pmesh2 = m.pcolormesh(x1,y1,satchl,cmap=cmap,norm=norm)
m.drawcoastlines(linewidth=0.25)
m.fillcontinents(color='lightgrey')
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='whitesmoke')

cbaxes = figure.add_axes([0.95, 0.2, 0.015, 0.175])
cb = ax.figure.colorbar(pmesh, cax = cbaxes)
'''
##################################################################
# read mixed water profile curves
Mchlsigma = np.zeros((20,5)) ; Msigma = np.zeros((20,5))
k=0
with open("./M_profiles", "r") as filestream:
     for line in filestream:
         currentline = line.split(",")

         Msigma[k,0] = float(currentline[1]); Msigma[k,1] = float(currentline[3]);
         Msigma[k,2] = float(currentline[5]); Msigma[k,3] = float(currentline[7]); Msigma[k,4] = float(currentline[9])
         Mchlsigma[k,0] = float(currentline[0]); Mchlsigma[k,1] = float(currentline[2]);
         Mchlsigma[k,2] = float(currentline[4]); Mchlsigma[k,3] = float(currentline[6]); Mchlsigma[k,4] = float(currentline[8]);
         k=k+1

Mchlave = np.array([0.244,0.592,0.885,1.881,6.32])
##################################################################

##################################################################
# read mld climatology
ncmld = NetCDFFile('./mld_DR003_c1m_reg2.0.nc')
ncmld.set_auto_mask(False)
mldlat = ncmld.variables['lat'][:]
mldlon = ncmld.variables['lon'][:]
mld3D = ncmld.variables['mld'][:,:,:]
ncmld.close()
for l in range(mldlon.shape[0]):
    if mldlon[l]>180. :
       mldlon[l] = -( 360. - mldlon[l] )
date = datetime.date(1970, 1, 1) + datetime.timedelta(int(time))
dd = int(date.strftime('%-j'))
mldt = np.array([-15,15,46,74,105,135,166,196,227,258,288,319,349,380])
##################################################################

print 'adjusting profiles'
for j in range(I1):
   for i in range(I2):
#      print j,i
##################################################################
# get mld from climatology
      JJ = np.abs(plon[j,i]-mldlon).argmin()
      II = np.abs(plat[j,i]-mldlat).argmin()
      mld = np.interp(dd,mldt,np.concatenate((mld3D[11,II,JJ],mld3D[:,II,JJ],mld3D[0,II,JJ]),axis=None))
      # masked mld have very high values, I keep it that way, so later water column is assigned "mixed" to avoid code malfunction
###################################################################

##################################################################
# set whether mixed or stratified
      surf = mapped[j,i] # copy the local satellite chl
      if np.ma.is_masked(surf):
        pass
      else:
        eup  = 4.61 / (0.041 + 0.04*surf) # taken from ECOSMO - kw=0.041 kb=0.04 ln(1%)=4.61

        stratified = 'false'
#                       # initial estimation of euphotic depth to determine whether mixed or not
#                       # later on Morel and Maritorena,2001 estimation will be used
        if eup > mld :
           stratified = 'true'
##################################################################

        if stratified :

          if surf <= 1.0 :
             A10 = 36.1 ; B10 = 0.357
             A15 = 42.0 ; B15 = 0.248
          else :
             A10 = 37.7 ; B10 = 0.615
             A15 = 43.5 ; B15 = 0.847

          Cbc = np.zeros((9)); sc = np.zeros((9)); Cmaxc = np.zeros((9)); sigma_maxc = np.zeros((9)); deltac = np.zeros((9)); sigma_avec = np.zeros((9))
          Cbc[0] = 0.471; sc[0] = 0.135; Cmaxc[0] = 1.572; sigma_maxc[0] = 0.969; deltac[0] = 0.393; sigma_avec[0] = 0.032
          Cbc[1] = 0.533; sc[1] = 0.172; Cmaxc[1] = 1.194; sigma_maxc[1] = 0.921; deltac[1] = 0.435; sigma_avec[1] = 0.062
          Cbc[2] = 0.428; sc[2] = 0.138; Cmaxc[2] = 1.015; sigma_maxc[2] = 0.905; deltac[2] = 0.630; sigma_avec[2] = 0.098
          Cbc[3] = 0.570; sc[3] = 0.173; Cmaxc[3] = 0.766; sigma_maxc[3] = 0.814; deltac[3] = 0.586; sigma_avec[3] = 0.158
          Cbc[4] = 0.611; sc[4] = 0.214; Cmaxc[4] = 0.676; sigma_maxc[4] = 0.663; deltac[4] = 0.539; sigma_avec[4] = 0.244
          Cbc[5] = 0.390; sc[5] = 0.109; Cmaxc[5] = 0.788; sigma_maxc[5] = 0.521; deltac[5] = 0.681; sigma_avec[5] = 0.347
          Cbc[6] = 0.569; sc[6] = 0.183; Cmaxc[6] = 0.608; sigma_maxc[6] = 0.452; deltac[6] = 0.744; sigma_avec[6] = 0.540
          Cbc[7] = 0.835; sc[7] = 0.298; Cmaxc[7] = 0.382; sigma_maxc[7] = 0.512; deltac[7] = 0.625; sigma_avec[7] = 1.235
          Cbc[8] = 0.188; sc[8] = 0.   ; Cmaxc[8] = 0.885; sigma_maxc[8] = 0.378; deltac[8] = 1.081; sigma_avec[8] = 2.953

          if surf <= sigma_avec[0] :
             Cb = Cbc[0]; s = sc[0]; Cmax = Cmaxc[0]; sigma_max = sigma_maxc[0]; delta = deltac[0]; sigma_ave = sigma_avec[0]
          elif surf >= sigma_avec[8] :
             Cb = Cbc[8]; s = sc[8]; Cmax = Cmaxc[8]; sigma_max = sigma_maxc[8]; delta = deltac[8]; sigma_ave = sigma_avec[8]
          else :
             Cb        = np.interp( surf,sigma_avec,Cbc )
             s         = np.interp( surf,sigma_avec,sc )
             sigma_max = np.interp( surf,sigma_avec,sigma_maxc) 
             Cmax      = np.interp( surf,sigma_avec,Cmaxc )
             delta     = np.interp( surf,sigma_avec,deltac )

          sigma   = np.arange(21)/10.
          c_sigma = Cb - s*sigma + Cmax*np.exp( -( (sigma - sigma_max)/delta )**2 ) # the shape of the dimensionless profile for stratified waters based on surface chl
          '''
          plt.plot(c_sigma,-sigma)
          plt.scatter(c_sigma,-sigma)
          '''
          chl_int_eup       = A10 * surf**B10 # integrated chl within the euphotic depth for stratified waters
          chl_ave_eup       = chl_int_eup / eup # mean chl within the euphotic depth for stratified waters
          chl_profile       = chl_ave_eup * c_sigma # chl profile estimated from satellite chl for stratified waters
          chl_profile_depth = sigma * eup # chl profile sample depth for stratified waters

          '''
          plt.plot(chl_profile,-chl_profile_depth)
          plt.scatter(chl_profile,-chl_profile_depth)
          '''
        else:

          A10 = 42.1 ; B10 = 0.538
          A15 = 58.5 ; B15 = 0.546

          sigma = np.zeros((20)) ; c_sigma = zeros((20))
          if surf <= Mchlave[0] :
            sigma   = Msigma[:,0]
            c_sigma = Mchlsigma[:,0]
          elif surf >= Mchlave[4] :
            sigma   = Msigma[:,4]
            c_sigma = Mchlsigma[:,4]
          else:
            for k in range(20):
                sigma[k]   = np.interp( surf,Mchlave,Msigma[k,:] )
                c_sigma[k] = np.interp( surf,Mchlave,Mchlsigma[k,:] ) # the shape of the dimensionless profile for mixed waters based on surface chl

          '''
          plt.plot(c_sigma,-sigma)
          plt.scatter(c_sigma,-sigma)
          '''
          chl_int_eup = A10 * surf**B10 # integrated chl within the euphotic depth for mixed waters
          chl_ave_eup = chl_int_eup / eup # mean chl within the euphotic depth for stratified waters
          chl_profile       = chl_ave_eup * c_sigma # chl profile estimated from satellite chl for mixed waters
          chl_profile_depth = sigma * eup # chl profile sample depth for mixed waters
          '''
          plt.plot(chl_profile,-chl_profile_depth)
          plt.scatter(chl_profile,-chl_profile_depth)
          '''



################################################
# take model data in 1D
        dia    = dia3D[:,j,i]
        fla    = fla3D[:,j,i]
        diachl = diachl3D[:,j,i]
        flachl = flachl3D[:,j,i]
        chl    = diachl + flachl
        d1D    = depth3D[:,j,i]

# interpolate estimated profile to model depth
        chl_estimated = np.interp(d1D,chl_profile_depth,chl_profile)

        chl_adjusted = np.zeros((dm))
        flachl_adjsuted = np.zeros((dm))
        diachl_adjsuted = np.zeros((dm))
        fla_adjsuted = np.zeros((dm))
        dia_adjsuted = np.zeros((dm))
        for k in range(dm):
            if d1D[k] <= np.max(chl_profile_depth) :
               chl_adjusted[k] = chl_estimated[k] * 0.3 + chl[k] * 0.7
            else :
               chl_adjusted[k] = chl[k]

            flachl_adjsuted[k] = ( flachl[k] / (flachl[k] + diachl[k]) ) * chl_adjusted[k]
            diachl_adjsuted[k] = ( diachl[k] / (flachl[k] + diachl[k]) ) * chl_adjusted[k]
            fla_adjsuted[k] = flachl_adjsuted[k] * ( fla[k] / flachl[k] )
            dia_adjsuted[k] = diachl_adjsuted[k] * ( dia[k] / diachl[k] )

        dianew[:,j,i] = dia_adjsuted
        flanew[:,j,i] = fla_adjsuted
        diachlnew[:,j,i] = diachl_adjsuted
        flachlnew[:,j,i] = flachl_adjsuted





newfile = '/cluster/work/users/cagyum/program_files/restart.new'
new_abfile = abfile.ABFileRestart(newfile,"w",idm=I2,jdm=I1)
new_abfile.write_header(i_abfile._iexpt,i_abfile._iversn,i_abfile._yrflag,i_abfile._sigver,i_abfile._nstep,i_abfile._dtime,i_abfile._thbase)

for keys in sorted( i_abfile.fields.keys() ) :
    fieldname = i_abfile.fields[keys]["field"]
    k         = i_abfile.fields[keys]["k"]
    t         = i_abfile.fields[keys]["tlevel"]
    field     = i_abfile.read_field(fieldname,k,t)
    if fieldname == "ECO_dia" :
       print("MODIFYING  %10s at level %3d at time=%d"%(fieldname,k,t))
       new_abfile.write_field(dianew[k-1,:,:],None,fieldname,k,t)
    elif fieldname == "ECO_fla" :
       print("MODIFYING  %10s at level %3d at time=%d"%(fieldname,k,t))
       new_abfile.write_field(flanew[k-1,:,:],None,fieldname,k,t)
    elif fieldname == "ECO_flac" :
       print("MODIFYING  %10s at level %3d at time=%d"%(fieldname,k,t))
       new_abfile.write_field(flachlnew[k-1,:,:],None,fieldname,k,t)
    elif fieldname == "ECO_diac" :
       print("MODIFYING  %10s at level %3d at time=%d"%(fieldname,k,t))
       new_abfile.write_field(diachlnew[k-1,:,:],None,fieldname,k,t)
    else:
       print("Copying %10s to level %3d at time=%d"%(fieldname,k,t))
       new_abfile.write_field(field,None,fieldname,k,t)

new_abfile.close()

'''
m = Basemap(width=8500000,height=8500000,
resolution='c',projection='stere',\
lat_ts=80,lat_0=80,lon_0=-40.,round='True')
x,y=m(plon,plat)
satlon,satlat = np.meshgrid(lon,lat)
x1,y1 = m(satlon,satlat)


cmin = -1; cmax=1.
cmap = plt.get_cmap('Spectral_r')
figure=plt.figure(figsize=(12,6))
ax=figure.add_subplot(211)
ax.set_position([0.0,0.05,0.5,0.95])
levels = MaxNLocator(nbins=15).tick_values(cmin, cmax)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
pmesh = m.pcolormesh(x,y,diachlnew[0,:,:]+flachlnew[0,:,:]-diachl3D[0,:,:]-flachl3D[0,:,:],cmap=cmap,norm=norm)
m.drawcoastlines(linewidth=0.25)
m.fillcontinents(color='lightgrey')
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='whitesmoke')

cbaxes = figure.add_axes([0.03, 0.2, 0.015, 0.175])
cb = ax.figure.colorbar(pmesh, cax = cbaxes)


cmin = 0; cmax=5.
ax2=figure.add_subplot(212)
ax2.set_position([0.5,0.05,0.5,0.95])

levels = MaxNLocator(nbins=15).tick_values(cmin, cmax)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
pmesh2 = m.pcolormesh(x,y,mapped,cmap=cmap,norm=norm)
m.drawcoastlines(linewidth=0.25)
m.fillcontinents(color='lightgrey')
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='whitesmoke')

cbaxes2 = figure.add_axes([0.95, 0.2, 0.015, 0.175])
cb2 = ax2.figure.colorbar(pmesh2, cax = cbaxes2)

'''






















