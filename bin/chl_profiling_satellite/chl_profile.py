import numpy as np
from netCDF4 import Dataset as NetCDFFile
import pickle
import abfile
import modeltools.hycom
import datetime
import warnings
import argparse
import getpass
import _BETZY_the_mapping_loop_with_ice
from numpy import dtype
import os
########################
warnings.filterwarnings('ignore')
user = getpass.getuser()

'''
how to use:

The code works with mandatory inputs. 
You need to specify your own following these examples:

--domain=/cluster/work/users/cagyum/TP5a0.06
--nerscdir=/cluster/home/cagyum/NERSC-HYCOM-CICE
--satdir=/cluster/work/users/cagyum/OCCCI/L3/copernicus_globcolor_daily_rep/chl/

--opendap_rep # this overwrites user satdir and gets the Copernicus L3 reprocessed data through opendap
              # you need to have a Copernicus account.
              # the code will use your username and password.
              # you can provide your user and password as a system variable 
              # e.g. in bash: export copernicus_user=XXX  && export copernicus_pass=XXX  

'''
def main(experiment,year,day,satdir,nerscdir,domain,debug,opendap_rep,opendap_nrt):

    global PARAMS_Ardyna, PARAMS_Caglar, depthm, plat, plon, regions, REGmask, profdep
    experiment = 'expt_'+experiment[0:2]+'.'+experiment[2]

    # get domain dimensions
    abgrid = abfile.ABFileGrid(domain + "/" + experiment + "/data/regional.grid","r")
    plon=abgrid.read_field("plon")
    plat=abgrid.read_field("plat")
    scpx=abgrid.read_field("scpx")
    scpy=abgrid.read_field("scpy")
    jdm,idm=plon.shape

    # get domain depth
    abdepth = abfile.ABFileBathy(domain  + "/" + experiment +  "/data/regional.depth.b", \
            "r",idm=idm,jdm=jdm)
    depthm=abdepth.read_field("depth")

    # get the profiling parameters
    f = open(nerscdir+'/bin/chl_profiling_satellite/Ardyna_parameters.pckl','rb')
    PARAMS_Ardyna = pickle.load(f)  
    f.close()

    f = open(nerscdir+'/bin/chl_profiling_satellite/Caglar_parameters.pckl','rb')
    PARAMS_Caglar = pickle.load(f)  
    f.close()

    # get Ardyna regions
    regions = assign_regions()

    # get HYCOM TP5 regions
    f=open(nerscdir + '/input/region_masks_TP5.pckl','rb')
    REGmask = pickle.load(f)
    f.close()

    profdep = 201 # estimated profile level number 


    # read in satellite data
    timecal = ( datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(day)-1) )
    if opendap_rep or opendap_nrt:
       if opendap_rep:
          satfile = 'https://{}:{}@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D'.format(str(os.environ['copernicus_user']),str(os.environ['copernicus_pass']))
       if opendap_nrt:
          satfile = 'https://{}:{}@nrt.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_nrt_l3-multi-4km_P1D'.format(str(os.environ['copernicus_user']),str(os.environ['copernicus_pass']))
       nc = NetCDFFile(satfile)
       satlat = nc.variables['lat'][:]
       satlon = nc.variables['lon'][:]
       time = nc.variables['time'][:]
       time_index = np.abs(time-(timecal - datetime.datetime(1900,1,1)).days).argmin()
       satchl = nc.variables['CHL'][time_index,:,:]
    else:
       globchl = satdir+year+'/'+year+'-'+\
                 str(timecal.month).zfill(2)+'-'+str(timecal.day).zfill(2)+'.nc'
       nc = NetCDFFile(globchl)
       satlat = nc.variables['lat'][:]
       satlon = nc.variables['lon'][:]
       satchl = nc.variables['CHL'][0,:,:]
    #
    satlon_2d,satlat_2d=np.meshgrid(satlon,satlat)
#    if debug:
#       kdfile = satdir[:-4]+'kd/'+year+'/'+year+'-'+\
#                 str(timecal.month).zfill(2)+'-'+str(timecal.day).zfill(2)+'.nc'
#       nckd = NetCDFFile(kdfile)
#       kd = nckd.variables['KD490'][0,:,:]
    
    # mask out high latitudes ( >=70 ) from mid September to avoid artificial high CHL due to angle of the sun
    if int(day)>=259:
        satchl = np.ma.masked_where(satlat_2d>=70.,satchl)
#        if debug:
#           kd  = np.ma.masked_where(satlat_2d>=70.,kd)
        
    # load restart file
    oldfile = domain + "/" + experiment + "/data/restart."+year+"_"+day.zfill(3)+"_00_0000.a"
    f = abfile.ABFileRestart(oldfile,"r",\
        idm=idm,jdm=jdm)
    kdm  = max(f.fieldlevels)

    dia3D = np.zeros((kdm,jdm,idm))
    fla3D = np.zeros((kdm,jdm,idm))
    diachl3D = np.zeros((kdm,jdm,idm))
    flachl3D = np.zeros((kdm,jdm,idm))
    dp3D     = np.zeros((kdm,jdm,idm))
    depth3D  = np.zeros((kdm,jdm,idm))

    dianew    = np.zeros((kdm,jdm,idm))
    diachlnew = np.zeros((kdm,jdm,idm))
    flanew    = np.zeros((kdm,jdm,idm))
    flachlnew = np.zeros((kdm,jdm,idm))
    estimated = np.zeros((profdep,jdm,idm))
    d_estimated = np.zeros((profdep,jdm,idm))
    for k in range(kdm):
       dia3D[k,:,:]    = f.read_field('ECO_dia',k+1)
       fla3D[k,:,:]    = f.read_field('ECO_fla',k+1)
       diachl3D[k,:,:] = f.read_field('ECO_diac',k+1)
       flachl3D[k,:,:] = f.read_field('ECO_flac',k+1)
       dp3D[k,:,:]     = f.read_field('dp',k+1)/modeltools.hycom.onem
       if k == 0:
          depth3D[k,:,:] = dp3D[k,:,:]/2.
       else:
          depth3D[k,:,:] = depth3D[k-1,:,:] + dp3D[k,:,:]/2.

       dianew[k,:,:] = dia3D[k,:,:]
       flanew[k,:,:] = fla3D[k,:,:]
       diachlnew[k,:,:] = diachl3D[k,:,:]
       flachlnew[k,:,:] = flachl3D[k,:,:]

    depthmf = np.asfortranarray(depthm)
    satchlf = np.asfortranarray(satchl)
    platf   = np.asfortranarray(plat)
    plonf   = np.asfortranarray(plon)
    satlatf = np.asfortranarray(satlat)
    satlonf = np.asfortranarray(satlon)
    scpxf   = np.asfortranarray(scpx)
    scpyf   = np.asfortranarray(scpy)
#    if debug:
#       kdf     = np.asfortranarray(kd)
    icefile = domain + "/" + experiment + "/data/cice/iced." + year + "-"+ \
            str(timecal.month).zfill(2)+'-'+str(timecal.day).zfill(2)+'-00000.nc'
    ic = NetCDFFile(icefile)
    icemask = ic.variables["iceumask"][:,:]
    icf = np.asfortranarray(icemask)

    satout = _BETZY_the_mapping_loop_with_ice.main(depthmf,scpxf,scpyf,platf,plonf,icf,satlatf,satlonf,satchlf)
    satout = np.ma.masked_where(satout>1000.,satout)
    satout = np.ma.masked_where(satout<0.,satout)
    satout = np.ma.masked_array(satout,depthm.mask)
#    if debug:
#       kdout = _the_mapping_loop.main(depthmf,scpxf,scpyf,platf,plonf,icf,satlatf,satlonf,kdf)
#       kdout = np.ma.masked_where(kdout>1E4,kdout)
#       kdout = np.ma.masked_array(kdout,depthm.mask) 
#       kdout = 1./kdout
##########################################################################################
# NOW BEGINS THE CHL PROFILING
#
##########################################################################################

    print('adjusting profiles')
    for j in range(jdm):
       for i in range(idm):

          surf = satout[j,i] # copy the local satellite chl
          if np.ma.is_masked(surf):
            pass
          else:
            # get the shape pf the profile
            cz = profiler( surf, plat[j,i], plon[j,i], year, day, plon,use_caglar_parameters=True)
            # convert it to chlorophyll a
            chl_profile       = cz * (surf/(cz[0]))
            chl_profile_depth = np.arange(profdep)

            ################################################
            # now we fit the model towards the estimated profile
            # I used a 30% fit towards the estimated profile

            # take model data in 1D
            dia    = dia3D[:,j,i]
            fla    = fla3D[:,j,i]
            diachl = diachl3D[:,j,i]
            flachl = flachl3D[:,j,i]
            chl    = diachl + flachl
            d1D    = depth3D[:,j,i]
    
            # interpolate estimated profile to model depth
            chl_estimated = np.interp(d1D,chl_profile_depth,chl_profile)
    
            chl_adjusted = np.zeros((kdm))
            flachl_adjsuted = np.zeros((kdm))
            diachl_adjsuted = np.zeros((kdm))
            fla_adjsuted = np.zeros((kdm))
            dia_adjsuted = np.zeros((kdm))
            for k in range(kdm):
                if d1D[k] <= np.max(chl_profile_depth) : # profile estimation is limited to the profdep depth
                                                         # we need to keep model values below as they are
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
            estimated[:,j,i] = chl_profile
            d_estimated[:,j,i] = chl_profile_depth

    ################################################
    # now copy the new information to a new restart file
    newfile = domain + "/" + experiment + "/data/restart."+year+"_"+day.zfill(3)+"_00_0000_NEW"
    new_abfile = abfile.ABFileRestart(newfile,"w",idm=idm,jdm=jdm)
    new_abfile.write_header(f._iexpt,f._iversn,f._yrflag,f._sigver,f._nstep,f._dtime,f._thbase)
    
    for keys in sorted( f.fields.keys() ) :
        fieldname = f.fields[keys]["field"]
        k         = f.fields[keys]["k"]
        t         = f.fields[keys]["tlevel"]
        field     = f.read_field(fieldname,k,t)
        if fieldname == "ECO_dia" :
           print("MODIFYING  %10s at level %3d at time=%d"%(fieldname,k,t))
           dummy = dianew[k-1,:,:]
           field[~depthm.mask] = dummy[~depthm.mask]
#           new_abfile.write_field(dianew[k-1,:,:],None,fieldname,k,t)
           new_abfile.write_field(field,None,fieldname,k,t)
        elif fieldname == "ECO_fla" :
           print("MODIFYING  %10s at level %3d at time=%d"%(fieldname,k,t))
#           new_abfile.write_field(flanew[k-1,:,:],None,fieldname,k,t)
           dummy = flanew[k-1,:,:]
           field[~depthm.mask] = dummy[~depthm.mask]
           new_abfile.write_field(field,None,fieldname,k,t)
        elif fieldname == "ECO_flac" :
           print("MODIFYING  %10s at level %3d at time=%d"%(fieldname,k,t))
#           new_abfile.write_field(flachlnew[k-1,:,:],None,fieldname,k,t)
           dummy = flachlnew[k-1,:,:]
           field[~depthm.mask] = dummy[~depthm.mask]
           new_abfile.write_field(field,None,fieldname,k,t)
        elif fieldname == "ECO_diac" :
           print("MODIFYING  %10s at level %3d at time=%d"%(fieldname,k,t))
#           new_abfile.write_field(diachlnew[k-1,:,:],None,fieldname,k,t)
           dummy = diachlnew[k-1,:,:]
           field[~depthm.mask] = dummy[~depthm.mask]
           new_abfile.write_field(field,None,fieldname,k,t)
        else:
           print("Copying %10s to level %3d at time=%d"%(fieldname,k,t))
           new_abfile.write_field(field,None,fieldname,k,t)
    f.close() 
    new_abfile.close()

    ################################################
    # now replace restart files
    os.rename(oldfile[:-2]+'.a',oldfile[:-2]+'_OLD.a')
    os.rename(oldfile[:-2]+'.b',oldfile[:-2]+'_OLD.b')
    os.rename(newfile+'.a',newfile[:-4]+'.a')
    os.rename(newfile+'.b',newfile[:-4]+'.b')    
   
    if debug:

       namencout = domain + "/" + experiment + "/data/restart."+year+"_"+day.zfill(3)+"_00_0000_NEW.nc" 
       ncout = NetCDFFile(namencout, "w", format="NETCDF4")

       ncout.createDimension("JJ",jdm)
       ncout.createDimension("II",idm)
       ncout.createDimension("dpth",kdm)
       ncout.createDimension("dpro",profdep)
       sat_mapped    = ncout.createVariable('sat_mapped', dtype('double').char, ("JJ","II"))
#       kd_mapped     = ncout.createVariable('kd_mapped', dtype('double').char, ("JJ","II"))
       latitude      = ncout.createVariable('latitude', dtype('double').char, ("JJ","II"))
       longitude     = ncout.createVariable('longitude', dtype('double').char, ("JJ","II"))
       chlold        = ncout.createVariable('chlold', dtype('double').char, ("dpth","JJ","II"))
       chlnew        = ncout.createVariable('chlnew', dtype('double').char, ("dpth","JJ","II"))
       profile       = ncout.createVariable('profile', dtype('double').char, ("dpro","JJ","II")) 
       depth         = ncout.createVariable('depth', dtype('double').char, ("dpth","JJ","II"))
       profile_depth = ncout.createVariable('profile_depth', dtype('double').char, ("dpro","JJ","II"))

       depth3D     = np.ma.masked_array(depth3D,np.repeat(depthm[np.newaxis,:], kdm, 0).mask)
       old         = diachl3D + flachl3D
       old         = np.ma.masked_array(old,np.repeat(depthm[np.newaxis,:], kdm, 0).mask)
       new         = diachlnew + flachlnew
       new         = np.ma.masked_array(new,np.repeat(depthm[np.newaxis,:], kdm, 0).mask)
       estimated   = np.ma.masked_array(estimated,np.repeat(depthm[np.newaxis,:], profdep, 0).mask)
       d_estimated = np.ma.masked_array(d_estimated,np.repeat(depthm[np.newaxis,:], profdep, 0).mask)

       ncout.variables["longitude"][:]  = plon
       ncout.variables["latitude"][:]   = plat
       ncout.variables["sat_mapped"][:] = satout
       ncout.variables["chlold"][:]     = old
       ncout.variables["chlnew"][:]     = new
       ncout.variables["depth"][:]      = depth3D
       ncout.variables["profile"][:]    = estimated
       ncout.variables["profile_depth"][:] = d_estimated 
#       ncout.variables["kd_mapped"][:] = kdout

       ncout.sync()
       ncout.close()
##########################################################################################
# NOW ENDS THE CHL PROFILING
#
##########################################################################################

def select_model_region(lat,lon):
    cooINDEX = abs( plat-lat ) + abs( plon-lon )
    JJ,II = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)

    location = 'NONE'
    for key in REGmask.keys():
        if REGmask[key][JJ,II] == 1: location = key

    return location

def profiler(surf,lat,lon,year,day, *args, **kwargs):

    timecal = ( datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(day)-1) )
    month = timecal.month

    cooINDEX = abs( plat-lat ) + abs( plon-lon )
    JJ,II = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)

    if regions[JJ,II] == 0 : location = 'CENTRAL'
    if regions[JJ,II] == 1 : location = 'BAFFIN'
    if regions[JJ,II] == 2 : location = 'HUDSON'
    if regions[JJ,II] == 3 : location = 'CANADA'
    if regions[JJ,II] == 4 : location = 'BEAUFORT'
    if regions[JJ,II] == 5 : location = 'CHUKCHI'
    if regions[JJ,II] == 6 : location = 'BERING'
    if regions[JJ,II] == 7 : location = 'RUSSIAN'
    if regions[JJ,II] == 8 : location = 'BARENTS'
    if regions[JJ,II] == 9 : location = 'NORDIC'    

    if np.ma.is_masked(regions[JJ,II]):
       print('Water depth below 50 metes, no profiling will be done')
    else:
      if surf < 0.1:
         bins = 'C1'
      elif surf >= 0.1 and surf < 0.3:
         bins = 'C2'
      elif surf >= 0.3 and surf < 0.5:
         bins = 'C3'
      elif surf >= 0.5 and surf < 0.7:
         bins = 'C4'
      elif surf >= 0.7 and surf < 1.0:
         bins = 'C5'
      elif surf >= 1.0 and surf < 3.0:
         bins = 'C6'
      elif surf >= 3.0 and surf < 8.0:
         bins = 'C7'
      elif surf >= 8.0:
         bins = 'C8'

      if surf >= 0.7 :
         season = 'ALL'
         area = 'FULL'
      else:
         if month >= 2 and month <= 4:
            season = 'PRE'
            area = 'FULL'
         elif month >= 5 and month <= 9:
            season = 'POST'
            if bins == 'C3' or bins == 'C4':
               area = 'FULL'
            if bins == 'C1' or bins == 'C2':
               area = location
         else:
            season = 'WINTER'
            area = 'FULL'
      
      cb    = PARAMS_Ardyna[season][area][bins]['cb']
      s     = PARAMS_Ardyna[season][area][bins]['s']
      cmax  = PARAMS_Ardyna[season][area][bins]['cmax']
      zmax  = PARAMS_Ardyna[season][area][bins]['zmax']
      delta = PARAMS_Ardyna[season][area][bins]['delta']

      if kwargs.get('use_caglar_parameters',None):
         region = 'outside'
         if season == 'ALL' or season == 'POST' or season == 'PRE' or season == 'WINTER':
            region = select_model_region(lat,lon)         
            if region == 'Barents' or region == 'NorwegianN' or region == 'NorwegianS':

                  try:
                      PARAM_caglar_1[region][season][bins]['cb']
                      # retrieve caglar's modified parameters 
                      cb    = PARAM_Caglar[region][season][bins]['cb']
                      s     = PARAM_Caglar[region][season][bins]['s']
                      cmax  = PARAM_Caglar[region][season][bins]['cmax']
                      zmax  = PARAM_Caglar[region][season][bins]['zmax']
                      delta = PARAM_Caglar[region][season][bins]['delta']
                  except Exception:
                      pass
        

         print(timecal,season,area,bins,region,'   model J:',JJ,' model I:',II)
      else:
         print(timecal,season,area,bins)
            

      z = np.arange(profdep)
      cz = cb - s*z + cmax * np.exp(-((z - zmax)/delta)**2)

      return cz

def assign_regions():
   # DEFINE REGIONS (www.biogeosciences.net/10/4383/2013/)
   '''
   the order of assigning regions is important. Do not rearrange.
   '''
   regions = depthm*0.

     # BAFFIN BAY - assigned number = 1
   regions[ np.logical_and( np.logical_and(plon>=-80.,plon<=-50),np.logical_and(plat>=70.,plat<=82.) )  ] = 1.
   regions[ np.logical_and( np.logical_and(plon>=-70.,plon<=-50),np.logical_and(plat>=63.,plat<=70.) )  ] = 1.
   regions[ np.logical_and( np.logical_and(plon>=-65.,plon<=-50),np.logical_and(plat>=45.,plat<=63.) )  ] = 1.
   regions = np.ma.masked_where(depthm<50.,regions)

     # HUDSON BAY - assigned number = 2
   regions[ np.logical_and( np.logical_and(plon>=-100.,plon<=-76.5),np.logical_and(plat>=45.,plat<=64.) )  ] = 2.
   regions = np.ma.masked_where(depthm<50.,regions)

     # CANADIAN ARCHIPELAGO - assigned number = 3
   regions[ np.logical_and( np.logical_and( np.logical_and(plon>=-120.,plon<=-60),np.logical_and(plat>=45.,plat<=79.) ),regions==0.)  ] = 3.
   regions[ np.logical_and( np.logical_and(plon>=-128.,plon<=-120),np.logical_and(plat>=45.,plat<=72.) )  ] = 3.
   regions = np.ma.masked_where(depthm<50.,regions)

     # BEAUFORT SEA - assigned number = 4
   regions[ np.logical_and( np.logical_and(plon>=-150.,plon<=-128),np.logical_and(plat>=60.,plat<=72.) )  ] = 4.
   regions = np.ma.masked_where(depthm<50.,regions)

     # CHUKCHI SEA - assigned number = 5
   regions[ np.logical_and( np.logical_and(plon>=-180.,plon<=-150),np.logical_and(plat>=67.,plat<=80.) )  ] = 5.
   regions = np.ma.masked_where(depthm<50.,regions)

     # BERING SEA - assigned number = 6
   regions[ np.logical_and( plon>=160.,np.logical_and(plat>=45.,plat<=67.) )  ] = 6.
   regions[ np.logical_and( plon<=-150.,np.logical_and(plat>=45.,plat<=67.) )  ] = 6.
   regions = np.ma.masked_where(depthm<50.,regions)

     # RUSSIAN SEAS - assigned number = 7
   regions[ np.logical_and( np.logical_and(plon>=69.,plon<=180),np.logical_and(plat>=70.,plat<=80.) )  ] = 7.
   regions[ np.logical_and( np.logical_and(plon>=66.,plon<=180),np.logical_and(plat>=70.,plat<=76.5) )  ] = 7.
   regions[ np.logical_and( np.logical_and(plon>=64.,plon<=180),np.logical_and(plat>=65.,plat<=76.) )  ] = 7.
   regions[ np.logical_and( np.logical_and(plon>=62.,plon<=180),np.logical_and(plat>=65.,plat<=76.) )  ] = 7.
   regions[ np.logical_and( np.logical_and(plon>=60.,plon<=180),np.logical_and(plat>=65.,plat<=76.) )  ] = 7.
   regions[ np.logical_and( np.logical_and(plon>=58.,plon<=180),np.logical_and(plat>=65.,plat<=75.) )  ] = 7.
   regions[ np.logical_and( np.logical_and(plon>=56.,plon<=180),np.logical_and(plat>=65.,plat<=75.) )  ] = 7.
   regions[ np.logical_and( np.logical_and(plon>=56.,plon<=180),np.logical_and(plat>=65.,plat<=74.) )  ] = 7.
   regions = np.ma.masked_where(depthm<50.,regions)

     # BARENTS SEA - assigned number = 8
   regions[ np.logical_and( np.logical_and( np.logical_and(plon>=17.,plon<=70),np.logical_and(plat>=65.5,plat<=80.) ),regions==0.)  ] = 8.
   regions[ np.logical_and( np.logical_and( np.logical_and(plon>=31.,plon<=45),np.logical_and(plat>=63,plat<=66.) ),regions==0.)  ] = 8.
   regions = np.ma.masked_where(depthm<50.,regions)

     # NORWEGIAN & GREENLAND SEA + NORTH ATLANTIC - assigned number = 9
   regions[ np.logical_and( np.logical_and( np.logical_and(plon>=-50.,plon<=17.),np.logical_and(plat>=45.,plat<=80.) ),regions==0.)  ] = 9.
   regions = np.ma.masked_where(depthm<50.,regions)

   return regions


if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('domain',  type=str)
    parser.add_argument('experiment',  type=str)
    parser.add_argument('year',  type=str)
    parser.add_argument('day',  type=str)
    parser.add_argument('--satdir', type=str, \
        default="/cluster/work/users/cagyum/OCCCI/L3/copernicus_globcolor_daily_rep/chl/", \
            help="satellite file directory, should have separate year folders")
    parser.add_argument('--nerscdir', type=str, \
        default="/cluster/home/cagyum/NERSC-HYCOM-CICE/", \
            help="user specific hycom directory")
    parser.add_argument('--debug', action="store_true",default=False)
    parser.add_argument('--opendap_rep', action="store_true",default=False)
    parser.add_argument('--opendap_nrt', action="store_true",default=False)
    args = parser.parse_args()

    main(args.experiment,args.year,args.day,satdir=args.satdir,nerscdir=args.nerscdir,domain=args.domain,debug=args.debug,opendap_rep=args.opendap_rep,opendap_nrt=args.opendap_nrt)
