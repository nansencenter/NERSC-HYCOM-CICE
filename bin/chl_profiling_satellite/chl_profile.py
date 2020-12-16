import numpy as np
from netCDF4 import Dataset as NetCDFFile
import pickle
import abfile
import modeltools.hycom
import datetime
import warnings
import argparse
import getpass
import _FRAM_the_mapping_loop_with_ice
from numpy import dtype
import os
########################
warnings.filterwarnings('ignore')
user = getpass.getuser()

def main(region,experiment,year,day,workdir,satdir,debug):

    experiment = 'expt_'+experiment[0:2]+'.'+experiment[2]
    
    if region == "TP0":
        region = "TP0a1.00"
    if region == "TP2":
        region = "TP2a0.10"
#    if region == "TP4":
#        region = "TP4a0.12"
    if region == "TP5":
        region = "TP5a0.06"
    if region == "NAT":
        region = "NATa1.00"
#    if region == "NA2":
#        region = "NA2a0.80"

    # get domain dimensions
    abgrid = abfile.ABFileGrid(workdir + user + "/" + \
        region + "/topo/regional.grid","r")
    plon=abgrid.read_field("plon")
    plat=abgrid.read_field("plat")
    scpx=abgrid.read_field("scpx")
    scpy=abgrid.read_field("scpy")
    jdm,idm=plon.shape

    # get domain depth
    abdepth = abfile.ABFileBathy(workdir + user + "/" + \
        region + "/" + experiment + "/data/regional.depth.b", \
            "r",idm=idm,jdm=jdm)
    depthm=abdepth.read_field("depth")

    # read in satellite data
    timecal = ( datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(day)-1) )
    globchl = satdir+year+'/'+year+'-'+\
                 str(timecal.month).zfill(2)+'-'+str(timecal.day).zfill(2)+'.nc'
    nc = NetCDFFile(globchl)
    satlat = nc.variables['lat'][:]
    satlon = nc.variables['lon'][:]
    satchl = nc.variables['CHL'][0,:,:]
    #
    satlon_2d,satlat_2d=np.meshgrid(satlon,satlat)
    if debug:
       kdfile = satdir[:-4]+'kd/'+year+'/'+year+'-'+\
                 str(timecal.month).zfill(2)+'-'+str(timecal.day).zfill(2)+'.nc'
       nckd = NetCDFFile(kdfile)
       kd = nckd.variables['KD490'][0,:,:]
    
    # mask out high latitudes ( >=70 ) from mid September to avoid artificial high CHL due to angle of the sun
    if int(day)>=259:
        satchl = np.ma.masked_where(satlat_2d>=70.,satchl)
        if debug:
           kd  = np.ma.masked_where(satlat_2d>=70.,kd)
        
    # load restart file
    oldfile = workdir + user + "/" + \
        region + "/" + experiment + "/data/"+region[0:3]+"restart."+year+"_"+day.zfill(3)+"_00_0000.a"
    f = abfile.ABFileRestart(oldfile,"r",\
        idm=idm,jdm=jdm)
    kdm  = max(f.fieldlevels)
    prok = 21
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
    estimated = np.zeros((prok,jdm,idm))
    d_estimated = np.zeros((prok,jdm,idm))
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
    if debug:
       kdf     = np.asfortranarray(kd)
    icefile = workdir + user + "/" + \
        region + "/" + experiment + "/data/cice/iced." + year + "-"+ \
            str(timecal.month).zfill(2)+'-'+str(timecal.day).zfill(2)+'-00000.nc'
    ic = NetCDFFile(icefile)
    icemask = ic.variables["iceumask"][:,:]
    icf = np.asfortranarray(icemask)

    satout = _FRAM_the_mapping_loop_with_ice.main(depthmf,scpxf,scpyf,platf,plonf,icf,satlatf,satlonf,satchlf)
    satout = np.ma.masked_where(satout>1000.,satout)
    satout = np.ma.masked_where(satout<0.,satout)
    satout = np.ma.masked_array(satout,depthm.mask)
    if debug:
       kdout = _the_mapping_loop.main(depthmf,scpxf,scpyf,platf,plonf,icf,satlatf,satlonf,kdf)
       kdout = np.ma.masked_where(kdout>1E4,kdout)
       kdout = np.ma.masked_array(kdout,depthm.mask) 
       kdout = 1./kdout
##############################
# NOW BEGINS THE CHL PROFILING
#
##############################

    # read mld climatology #############################################
    ncmld = NetCDFFile('./mld_DR003_c1m_reg2.0.nc')
    ncmld.set_auto_mask(False)
    mldlat = ncmld.variables['lat'][:]
    mldlon = ncmld.variables['lon'][:]
    mld3D = ncmld.variables['mld'][:,:,:]
    ncmld.close()
    for l in range(mldlon.shape[0]):  # make longitudes compatible with hycom model
        if mldlon[l]>180. :
           mldlon[l] = -( 360. - mldlon[l] )
    # following two lines will be used to interpolate mld between months
    dayofyear = int(timecal.strftime('%-j'))                                
    mldt = np.array([-15,15,46,74,105,135,166,196,227,258,288,319,349,380]) # day of year every 15th 
    ####################################################################

    # read mixed water profile curves #################################
    # below are parameters taken from Uitz et al., 2006
    # they are stored in a text file and here we read them
    # Uitz classify profiles are 'mixed' and 'stratified'
    # this chunk simply reads in for 'mixed' profiles
    # we will determine whether the profile is mixed later
    Mchlsigma = np.zeros((20,5)) ; Msigma = np.zeros((20,5))
    k=0
    with open("./M_profiles", "r") as filestream:
         for line in filestream:
             currentline = line.split(",")
    
             Msigma[k,0] = float(currentline[1]); 
             Msigma[k,1] = float(currentline[3]);
             Msigma[k,2] = float(currentline[5]); 
             Msigma[k,3] = float(currentline[7]); 
             Msigma[k,4] = float(currentline[9])
             Mchlsigma[k,0] = float(currentline[0]); 
             Mchlsigma[k,1] = float(currentline[2]);
             Mchlsigma[k,2] = float(currentline[4]); 
             Mchlsigma[k,3] = float(currentline[6]); 
             Mchlsigma[k,4] = float(currentline[8]);
             k=k+1
    
    Mchlave = np.array([0.244,0.592,0.885,1.881,6.32])
    zeu_ave = np.array([77.1,53.2,44.,31.5,16.9]) # average EZ depth for the mixing conditions - taken from Uitz et al., 2006
    ####################################################################


# NOW WE BEGIN THE LOOP FOR EACH MODEL POINT

    print 'adjusting profiles'
    for j in range(jdm):
       for i in range(idm):
    #      print j,i
    ##################################################################
    # get mld from climatology
          JJ = np.abs(plon[j,i]-mldlon).argmin()
          II = np.abs(plat[j,i]-mldlat).argmin()
          # for each model point, interpolate MLD values to the exact dayofyear
          # for continuity, last month is concatenated to the front of the year
          # and first month to the back of the year
          mld = np.interp(dayofyear,mldt,\
                np.concatenate((mld3D[11,II,JJ],mld3D[:,II,JJ],mld3D[0,II,JJ]),axis=None))
          # masked mld have very high values, I keep it that way, 
          # so later water column is assigned "mixed" to avoid code malfunction
    ###################################################################

    ##################################################################
    # true/false whether mixed or stratified
          surf = satout[j,i] # copy the local satellite chl
          if np.ma.is_masked(surf):
            pass
          else:
            eup  = 4.61 / (0.041 + 0.04*surf) # taken from ECOSMO - kw=0.041 kb=0.04 ln(1%)=4.61
    
            stratified = 'false'
            # initial estimation of euphotic depth to determine whether mixed or not
            # later on in the code, Morel and Maritorena,2001 estimation will be used
            if eup > mld :
               stratified = 'true'

            if stratified : # if stratified, we need to assign different parameters (from Uitz et al., 2006)
                            # I suggest at this point, it is worth to look at the paper and have an idea on 
                            # what is happening below. Basically a different gaussian shaped curve is assigned
                            # based on the surface CHL value. From what I understand, it is safe to interpolate
                            # between curves, so any surface CHL value from satellite will be fit to an individual
                            # interpolated profile     
              if surf <= 1.0 :
                 A10 = 36.1 ; B10 = 0.357
                 A15 = 42.0 ; B15 = 0.248
              else :
                 A10 = 37.7 ; B10 = 0.615
                 A15 = 43.5 ; B15 = 0.847
    
              Cbc = np.zeros((9)); sc = np.zeros((9)); Cmaxc = np.zeros((9)); 
              sigma_maxc = np.zeros((9)); deltac = np.zeros((9)); sigma_avec = np.zeros((9)); zeu_ave = np.zeros((9))
              Cbc[0] = 0.471; sc[0] = 0.135; Cmaxc[0] = 1.572; sigma_maxc[0] = 0.969; deltac[0] = 0.393; sigma_avec[0] = 0.032; zeu_ave[0] = 119.1
              Cbc[1] = 0.533; sc[1] = 0.172; Cmaxc[1] = 1.194; sigma_maxc[1] = 0.921; deltac[1] = 0.435; sigma_avec[1] = 0.062; zeu_ave[1] = 99.9
              Cbc[2] = 0.428; sc[2] = 0.138; Cmaxc[2] = 1.015; sigma_maxc[2] = 0.905; deltac[2] = 0.630; sigma_avec[2] = 0.098; zeu_ave[2] = 91.0
              Cbc[3] = 0.570; sc[3] = 0.173; Cmaxc[3] = 0.766; sigma_maxc[3] = 0.814; deltac[3] = 0.586; sigma_avec[3] = 0.158; zeu_ave[3] = 80.2
              Cbc[4] = 0.611; sc[4] = 0.214; Cmaxc[4] = 0.676; sigma_maxc[4] = 0.663; deltac[4] = 0.539; sigma_avec[4] = 0.244; zeu_ave[4] = 70.3
              Cbc[5] = 0.390; sc[5] = 0.109; Cmaxc[5] = 0.788; sigma_maxc[5] = 0.521; deltac[5] = 0.681; sigma_avec[5] = 0.347; zeu_ave[5] = 63.4
              Cbc[6] = 0.569; sc[6] = 0.183; Cmaxc[6] = 0.608; sigma_maxc[6] = 0.452; deltac[6] = 0.744; sigma_avec[6] = 0.540; zeu_ave[6] = 54.4
              Cbc[7] = 0.835; sc[7] = 0.298; Cmaxc[7] = 0.382; sigma_maxc[7] = 0.512; deltac[7] = 0.625; sigma_avec[7] = 1.235; zeu_ave[7] = 39.8
              Cbc[8] = 0.188; sc[8] = 0.   ; Cmaxc[8] = 0.885; sigma_maxc[8] = 0.378; deltac[8] = 1.081; sigma_avec[8] = 2.953; zeu_ave[8] = 26.1
    
              if surf <= sigma_avec[0] :
                 Cb = Cbc[0]; s = sc[0]; Cmax = Cmaxc[0]; sigma_max = sigma_maxc[0]; delta = deltac[0]; sigma_ave = sigma_avec[0]; zeu = zeu_ave[0]
              elif surf >= sigma_avec[8] :
                 Cb = Cbc[8]; s = sc[8]; Cmax = Cmaxc[8]; sigma_max = sigma_maxc[8]; delta = deltac[8]; sigma_ave = sigma_avec[8]; zeu = zeu_ave[8]
              else :
                 Cb        = np.interp( surf,sigma_avec,Cbc )
                 s         = np.interp( surf,sigma_avec,sc )
                 sigma_max = np.interp( surf,sigma_avec,sigma_maxc)
                 Cmax      = np.interp( surf,sigma_avec,Cmaxc )
                 delta     = np.interp( surf,sigma_avec,deltac )
                 zeu       = np.interp( surf,sigma_avec,zeu_ave )
              sigma   = np.arange(prok)/10.
              # the shape of the dimensionless profile for stratified waters based on surface chl
              c_sigma = Cb - s*sigma + Cmax*np.exp( -( (sigma - sigma_max)/delta )**2 ) 
              '''
              plt.plot(c_sigma,-sigma)
              plt.scatter(c_sigma,-sigma)
              '''
              chl_int_eup       = A10 * surf**B10 # integrated chl within the euphotic depth for stratified waters
              #chl_ave_eup       = chl_int_eup / eup # mean chl within the euphotic depth for stratified waters
              chl_ave_eup       = chl_int_eup / zeu 
              chl_profile       = chl_ave_eup * c_sigma # chl profile estimated from satellite chl for stratified waters
              #chl_profile_depth = sigma * eup # chl profile sample depth for stratified waters
              chl_profile_depth = sigma * zeu
              '''
              plt.plot(chl_profile,-chl_profile_depth)
              plt.scatter(chl_profile,-chl_profile_depth)
              '''
            else: # if water column is mixed we use below CHL profile shapes 
    
              A10 = 42.1 ; B10 = 0.538
              A15 = 58.5 ; B15 = 0.546
    
              sigma = np.zeros((20)) ; c_sigma = zeros((20))
              if surf <= Mchlave[0] :
                sigma   = Msigma[:,0]
                c_sigma = Mchlsigma[:,0]
                zeu     = zeu_ave[0]
              elif surf >= Mchlave[4] :
                sigma   = Msigma[:,4]
                c_sigma = Mchlsigma[:,4]
                zeu     = zeu_ave[4]
              else:
                for k in range(20):
                    sigma[k]   = np.interp( surf,Mchlave,Msigma[k,:] )
                    c_sigma[k] = np.interp( surf,Mchlave,Mchlsigma[k,:] ) # the shape of the dimensionless profile for mixed waters based on surface chl
                    zeu = np.interp( surf,Mchlave,zeu_ave )
              '''
              plt.plot(c_sigma,-sigma)
              plt.scatter(c_sigma,-sigma)
              '''
              chl_int_eup = A10 * surf**B10 # integrated chl within the euphotic depth for mixed waters
              #chl_ave_eup = chl_int_eup / eup # mean chl within the euphotic depth for stratified waters
              chl_ave_eup       = chl_int_eup / zeu
              chl_profile       = chl_ave_eup * c_sigma # chl profile estimated from satellite chl for mixed waters
              #chl_profile_depth = sigma * eup # chl profile sample depth for mixed waters
              chl_profile_depth = sigma * zeu
              '''
              plt.plot(chl_profile,-chl_profile_depth)
              plt.scatter(chl_profile,-chl_profile_depth)
              '''
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
                if d1D[k] <= np.max(chl_profile_depth) : # profile estimation is limited to the euphotic depth
                                                         # we need to keep model values below
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
    
#    dianew[dianew>1E25] = 0.; dianew[dianew<-1E15] = 0.
#    dianew[diachlnew>1E25] = 0.; dianew[diachlnew<-1E15] = 0.
#    flanew[flanew>1E25] = 0.; flanew[flanew<-1E15] = 0.
#    flanew[flachlnew>1E25] = 0.; flanew[flachlnew<-1E15] = 0.

#    dianew = np.ma.masked_array(dianew,np.repeat(depthm[np.newaxis,:], kdm, 0).mask)
#    flanew = np.ma.masked_array(flanew,np.repeat(depthm[np.newaxis,:], kdm, 0).mask)
#    flachlnew = np.ma.masked_array(flachlnew,np.repeat(depthm[np.newaxis,:], kdm, 0).mask)
#    diachlnew = np.ma.masked_array(diachlnew,np.repeat(depthm[np.newaxis,:], kdm, 0).mask)
    # finally replace the old model values with the estimated one in the restart file
    newfile = workdir + user + "/" + \
        region + "/" + experiment + "/data/"+region[0:3]+"restart."+year+"_"+day.zfill(3)+"_00_0000_NEW"
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

# now replace restart files
    os.rename(oldfile[:-2]+'.a',oldfile[:-2]+'_OLD.a')
    os.rename(oldfile[:-2]+'.b',oldfile[:-2]+'_OLD.b')
    os.rename(newfile+'.a',newfile[:-4]+'.a')
    os.rename(newfile+'.b',newfile[:-4]+'.b')    
   
    if debug:

       namencout = workdir + user + "/" + \
        region + "/" + experiment + "/data/"+region[0:3]+"restart."+year+"_"+day.zfill(3)+"_00_0000_NEW.nc" 
       ncout = NetCDFFile(namencout, "w", format="NETCDF4")

       ncout.createDimension("JJ",jdm)
       ncout.createDimension("II",idm)
       ncout.createDimension("dpth",kdm)
       ncout.createDimension("dpro",prok)
       sat_mapped    = ncout.createVariable('sat_mapped', dtype('double').char, ("JJ","II"))
       kd_mapped     = ncout.createVariable('kd_mapped', dtype('double').char, ("JJ","II"))
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
       estimated   = np.ma.masked_array(estimated,np.repeat(depthm[np.newaxis,:], prok, 0).mask)
       d_estimated = np.ma.masked_array(d_estimated,np.repeat(depthm[np.newaxis,:], prok, 0).mask)

       ncout.variables["longitude"][:]  = plon
       ncout.variables["latitude"][:]   = plat
       ncout.variables["sat_mapped"][:] = satout
       ncout.variables["chlold"][:]     = old
       ncout.variables["chlnew"][:]     = new
       ncout.variables["depth"][:]      = depth3D
       ncout.variables["profile"][:]    = estimated
       ncout.variables["profile_depth"][:] = d_estimated 
       ncout.variables["kd_mapped"][:] = kdout

       ncout.sync()
       ncout.close()
    
if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('region',  type=str)
    parser.add_argument('experiment',  type=str)
    parser.add_argument('year',  type=str)
    parser.add_argument('day',  type=str)
    parser.add_argument('--workdir', type=str, \
        default="/cluster/work/users/", \
            help="machine specific work directory (above user folder)")
    parser.add_argument('--satdir', type=str, \
        default="/cluster/work/users/cagyum/OCCCI/L3/copernicus_globcolor_daily_rep/chl/", \
            help="satellite file directory, should have separate year folders")
    parser.add_argument('--debug', action="store_true",default=False)
    args = parser.parse_args()

    main(args.region,args.experiment,args.year,args.day,workdir=args.workdir,satdir=args.satdir,debug=args.debug)


