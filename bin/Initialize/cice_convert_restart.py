#.!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import sys
from scipy.interpolate import griddata as grd
import matplotlib.pyplot as plt

# testing under /home/sm_alfal/work/sea/TOPAZ5/Forecast14/Forecast14_000/cice/

#-------------------------------------------------------------------------
#   def set_to_val(nc,variable,icemask,val=0.0,n=-999):
#       var    = nc.variables[variable]
#       if n < 0:
#           var[:] = val*icemask
#       else:
#           var[n,:] = val*icemask
#   
# Function that reprojects model into observational grid
#def reproj_mod2obs(X1,Y1,Z1,X2,Y2,method='nearest',mask=None):
def interpolate2points(X1,Y1,Z1,X2,Y2,method='nearest',mask=None):
    #  X1,Y1 are 2D arrays holding input mesh    
    #  Z1 is array to interpolate, has same shape as X1 and Y1
    #  X2,Y2 are 2D arrays holding the target mesh.
    # This function will return an array holding the interplated values 
    #print(f'Interpolate method= {method}')
    minv=np.min(Z1)
    maxv=np.max(Z1)
    X1d = X1.flatten()
    Y1d = Y1.flatten()
    Z1d = Z1.flatten()
    # Interpolation
    # - can be done with other methods ('nearest','linear','cubic'<--doesn't work for our data)
    C  = np.array([X1d, Y1d]).T
    Z2 = grd(points=C,values=Z1d,xi=(X2, Y2), method=method)
   # Perform the interpolation using griddata
    #fine_values = griddata(coarse_points,coarse_values,xi=(fine_lon_mesh, fine_lat_mesh),method='linear' ) # You can change to 'nearest' or 'cubic' if needed
    Z2D = Z2.reshape(X2.shape)
    Z2D = np.clip(Z2D, minv, maxv)
    #print(f'Output array shape= {Z2D.shape} --- Interpolate method= {method}')
    # check output and add mask
    #mask2 = np.isnan(Z2)
    #if mask is not None:
    #    # apply union of mask and model nans
    #    mask2 = np.logical_or(mask1, mask)
    #Z2 = np.ma.array(Z2, mask=mask2)
    return(Z2D)


#def maplev(a,lpp=1):
def landfill(a,lpp=1):
    # gapfilling method
    jm,im=a.shape
    J,I=numpy.where(~numpy.isnan(a))
    with numpy.errstate(invalid='ignore'):
        av=numpy.nansum(a[J,I])/(len(I)*len(J))
    J,I=numpy.where(numpy.isnan(a))
    a[J,I]=av
    b=a
    #lpp=1  # it is better to be set to 100, but for practial reasnon, we keep it very small for now
    i=range(1,im-1)
    j=range(1,jm-1)
    ip1=range(2,im)
    jp1=range(2,jm)
    im1=range(0,im-2)
    jm1=range(0,jm-2)
    
    cc=numpy.zeros(a.shape)
    for k in range(lpp):
        cc[1:-2,1:-2]=b[1:-2,1:-2]+.5/4*( b[1:-2,2:-1]+b[0:-3,1:-2]+b[1:-2,0:-3]+b[2:-1,1:-2]-4.*b[1:-2,1:-2] )
        cc[:,0]=cc[:,1]
        cc[:,im-1]=cc[:,im-2]
        cc[jm-1,:]=cc[jm-2,:]
        b[J,I]=cc[J,I]

    a[J,I]=cc[J,I]
    J,I=numpy.where(numpy.isnan(a))
    a[J,I]=0.
    return a
#-------------------------------------------------------------------------
        

# Make initialization of new cice-file...

if len(sys.argv[:]) < 3:
    print("Usage: python cice_convert_restart.py <input_cice_restartfile> <input_cice_grid_file> <target_cice_gridfile>")
    sys.exit()

#python ../bin/cice_convert_restart.py iced.2024-12-08-00000_mem000.nc cice_grid_coarse.nc cice_grid_fine.nc
#hycomfile = sys.argv[1] #'/disk1/CICE_ini/ocean_ini.nc'
input_cicefile = sys.argv[1] #'iced.1997-01-01-00000.nc'
input_cicegridfile=sys.argv[2] #'cice_grid.nc'
target_cicegridfile=sys.argv[3] #'cice_grid.nc'

output_cicefile = 'interpolated_{}'.format(input_cicefile)
#print('hycom file is {}'.format(hycomfile))
print('input cice restart file is {}'.format(input_cicefile))
print('input cice grid    file is {}'.format(input_cicegridfile))
print('target cice grid   file is {}'.format(target_cicegridfile))
print(' ')
print(f'Output cice restart file is {output_cicefile}')
print(' ')
#nc_hycom   = netCDF4.Dataset(hycomfile)
#print(nc_hycom.file_format)
#nc_cice   = netCDF4.Dataset(cicefile,'r+')
#nc_input_gr = netCDF4.Dataset(input_cicegridfile)
with Dataset(input_cicegridfile,mode='r') as nc_input_gr:
    input_lat = nc_input_gr.variables['ulat'][:]
    input_lon = nc_input_gr.variables['ulon'][:]
    ny_in, nx_in = input_lat.shape
    print('shape of input cice grid is {}x{}'.format(ny_in, nx_in))
#nc_input_gr.close()

#nc_trgt_gr = netCDF4.Dataset(target_cicegridfile)
with Dataset(target_cicegridfile,mode='r') as nc_trgt_gr:
    trgt_lat = nc_trgt_gr.variables['ulat'][:]
    trgt_lon = nc_trgt_gr.variables['ulon'][:]
    trgt_kmt = nc_trgt_gr.variables['kmt'][:]
    ny_trgt, nx_trgt = trgt_lat.shape
    print('shape of target cice grid is {}x{}'.format(ny_trgt, nx_trgt))
#nc_trgt_gr.close()


#mask_h = np.where(aice_hi[:] > 0.1, 1, 0)
#aicen_c= out_nc_cice.variables['aicen']  # ice concentration per category in grid cell
NCAT=5
output_nc_cice = Dataset(output_cicefile,mode='w', format='NETCDF4')
output_nc_cice.createDimension('ni', nx_trgt)
output_nc_cice.createDimension('nj', ny_trgt)
output_nc_cice.createDimension('ncat', NCAT)
# calculate ice mask

#NCAT = [0, 0.6445072, 1.391433, 2.470179, 4.567288, 1e+08] #upper limits of ice categories
varlist2d = ['uvel','vvel','scale_factor','swvdr','swvdf','swidr','swidf','strocnxT','strocnyT',
             'stressp_1','stressp_2','stressp_3','stressp_4','stressm_1','stressm_2','stressm_3',
             'stressm_4','stress12_1','stress12_2','stress12_3','stress12_4','iceumask','frz_onset','fsnow']
#
varlist3d = ['aicen','vicen','vsnon','Tsfcn','iage','FY','alvl', 'vlvl','apnd','hpnd','ipnd','dhs','ffrac']
#varlist3d_sice = ['sice001', 'sice002', 'sice003', 'sice004', 'sice005', 'sice006', 'sice007','qsno001']
#sice_val = 5.0
#varlist3d_qice = ['qice001', 'qice002', 'qice003', 'qice004', 'qice005', 'qice006', 'qice007']
#qice_val = -2.0e8
#varlist3d_lvl = ['alvl', 'vlvl']
#lvl_val = 1.0
varlist3d_ice = ['sice001','qice001', 'sice002','qice002', 'sice003','qice003', 'sice004','qice004',
        'sice005','qice005', 'sice006','qice006', 'sice007','qice007','qsno001']


input_nc_cice  = Dataset(input_cicefile,mode='r')
global_attributes = {attr: input_nc_cice.getncattr(attr) for attr in input_nc_cice.ncattrs()}

icemask_in =input_nc_cice.variables['iceumask'][:]
icemask_interp=interpolate2points(input_lon,input_lat,icemask_in,trgt_lon,trgt_lat,method='nearest')
icemask_interp=icemask_interp*trgt_kmt
# loop over 2Dlist
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loop over 2D variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ')
for ss in varlist2d:
    print(f'-------Reading and interpolating varaiable= {ss}------')
    var2d =input_nc_cice.variables[ss]
    var2d_in=var2d[:]
    var2d_in_interp=interpolate2points(input_lon,input_lat,var2d_in,trgt_lon,trgt_lat,method='nearest')
    #var2d_in_interp=interpolate2points(input_lon,input_lat,var2d_in,trgt_lon,trgt_lat,method='linear')
    print(f'var2d_data_shape= {var2d_in_interp.shape} ')
    # Create a new variable in the new file
    new_var2d = output_nc_cice.createVariable(ss, var2d.dtype, ('nj', 'ni'))
    new_var2d[:,:] = var2d_in_interp*icemask_interp
    # Optionally, copy attributes of the variable
    #for attr_name in variable.ncattrs():
    #    setattr(new_var, attr_name, getattr(variable, attr_name))
    #set_to_val(nc_cice,s,mask_h,0.0)
    
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loop over 3D variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ')
# loop over 2Dlist
for ss in varlist3d + varlist3d_ice :
    print(f'-------Reading and interpolating 3D varaiable= {ss}------')
    var3d =input_nc_cice.variables[ss]
    var3d_in=var3d[:]
    print('shape of 3D var {}'.format(var3d_in.shape))
    new_var3d = output_nc_cice.createVariable(ss, var3d.dtype, ('ncat','nj', 'ni'))
    # Create a new variable in the new file
    for n in range(NCAT):  #N_category=5
        print('Category: '+str(n+1))
        #print('shape of 3D var ss.shape {}'.format(var3d_in[n,:].shape))
        var2d_in_interp=interpolate2points(input_lon,input_lat,var3d_in[n,:],trgt_lon,trgt_lat,method='nearest')
        #print(f'var2d_data_shape= {var2d_in_interp.shape} ')
        new_var3d[n,:] = var2d_in_interp[:]*icemask_interp
#AAA        # Optionally, copy attributes of the variable
#for n in range(N_category):  #N_category=5
#    print('Category: '+str(n+1))
#        #set_to_val(nc_cice,s,mask_h,0.0,n)
#AAA        for s in varlist3d_sice:
#AAA            set_to_val(nc_cice,s,mask_h,sice_val,n)
#AAA        for s in varlist3d_qice:
#AAA            set_to_val(nc_cice,s,mask_h,qice_val,n)
#AAA        for s in varlist3d_lvl:
#        set_to_val(nc_cice,s,mask_h,lvl_val,n)
    #set_to_val(nc_cice,'Tsfcn',mask_h,-3.0,n)
#nc_cice.istep1 = 0
#nc_cice.time = 0
#nc_cice.time_forc = 0
#nc_cice.nyr = 1
#nc_cice.month = 1
#nc_cice.mday = 1
#nc_cice.sec = 0
for attr_name, attr_value in global_attributes.items():
    setattr(output_nc_cice, attr_name, attr_value)

output_nc_cice.close()
input_nc_cice.close()

print(f'Output cice restart file is {output_cicefile}')
#from mpl_toolkilts.basemap import Basemap
#bnd_lat=50.0
#lon0=-45.

#lev=np.array([0.2,0.5,0.7,1.0])
#plt.figure()
#m = Basemap(projection='npstere', lon_0=lon0, boundinglat=bnd_lat, resolution='l')
#x, y = m(lon_c, lat_c)
#m.drawcoastlines()
#m.fillcontinents(color='lightgrey')
#cs = m.contourf(x,y,aice_hi, levels=lev, cmap=plt.cm.bone)
#cbar=m.colorbar(cs)
#
#plt.figure()
#m = Basemap(projection='npstere', lon_0=lon0, boundinglat=bnd_lat, resolution='l')
#x, y = m(lon_h, lat_h)
#m.drawcoastlines()
#m.fillcontinents(color='lightgrey')
#cs = m.contourf(x,y,aice_h2, levels=lev, cmap=plt.cm.bone)
#cbar=m.colorbar(cs)

if False:
    plt.figure()
    plt.imshow(aice_hi, cmap=plt.cm.bone)
    plt.clim(0,1)
    plt.colorbar()
    plt.title('aice int')
    plt.figure()
    plt.imshow(aice_h2, cmap=plt.cm.bone)
    plt.clim(0,1)
    plt.colorbar()
    plt.title('aice topaz')
    plt.figure()
    plt.imshow(hice_hi, cmap=plt.cm.Blues)
    plt.clim(0,3)
    plt.colorbar()
    plt.title('hice int')
    plt.figure()
    plt.imshow(hice_h2, cmap=plt.cm.Blues)
    plt.clim(0,3)
    plt.colorbar()
    plt.title('hice topaz')

plt.show()

print('done')



#AAA    # convert masked to array to array
#AAA    aice_h = aice_h.filled(fill_value=0.0)
#AAA    hice_h = hice_h.filled(fill_value=0.0)
#AAA    # topaz file needs to be averaged
#AAA    aice_h2 = (0.5*(aice_h[0,:,:]+aice_h[1,:,:]))
#AAA    hice_h2 = (0.5*(hice_h[0,:,:]+hice_h[1,:,:]))
#AAA    ny_h, nx_h = aice_h2.shape
#AAA    print(type(aice_h2))
#AAA    print(type(hice_h2))
#AAA  
#AAA    # cice stuff
#AAA    aicen_c   = nc_cice.variables['aicen']  # ice concentration per category in grid cell
#AAA    vicen_c   = nc_cice.variables['vicen']  # volume per unit area of ice (in category n)
#AAA    tsfcn_c   = nc_cice.variables['Tsfcn']
#AAA    
#AAA    print('shape of topaz is {}x{}'.format(ny_h, nx_h))
#AAA    print aicen_c.shape
#AAA    
#AAA    aice_hi=np.zeros((ny_c, nx_c))
#AAA    hice_hi=np.zeros((ny_c, nx_c))
#AAA    
#AAA    # average over navg adjacent grid points on hycom grid
#AAA    navg=2
#AAA        
#AAA    for j in range(ny_c):
#AAA        for i in range(nx_c):
#AAA            distance = np.sqrt( (lat_c[j,i]-lat_h)**2 + (lon_c[j,i]-lon_h)**2 )
#AAA            j_h, i_h = np.where(distance == distance.min())
#AAA            j_h = int(j_h[0])
#AAA            i_h = int(i_h[0])
#AAA            if j_h > navg-1 and j_h < ny_h-navg+1 and i_h > navg-1 and i_h < nx_h-navg+1 :
#AAA                aice_hi[j,i] = np.mean(aice_h2[j_h-navg:j_h+navg,i_h-navg:i_h+navg])
#AAA                hice_hi[j,i] = np.mean(hice_h2[j_h-navg:j_h+navg,i_h-navg:i_h+navg])
#AAA            else :
#AAA                aice_hi[j,i] = aice_h2[j_h,i_h]
#AAA                hice_hi[j,i] = hice_h2[j_h,i_h]
#AAA    
#AAA            for n in range(len(NCAT[1:])):
#AAA                if hice_hi[j,i] > NCAT[n] and hice_hi[j,i] < NCAT[n+1] :
#AAA                    aicen_c[n,j,i] = aice_hi[j,i]
#AAA                    vicen_c[n,j,i] = aice_hi[j,i]*hice_hi[j,i]
#AAA                else :
#AAA                    aicen_c[n,j,i] = 0.0
#AAA                    vicen_c[n,j,i] = 0.0
#AAA                if aice_hi[j,i] > 0.1 :
#AAA                    tsfcn_c[n,j,i] = -5.0
#AAA    

  
