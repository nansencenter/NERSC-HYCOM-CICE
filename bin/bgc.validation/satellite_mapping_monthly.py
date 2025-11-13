import numpy as np
from netCDF4 import Dataset as NetCDFFile
import abfile
import getpass
#import _the_loop
import _FRAM_the_mapping_loop
import modeltools.hycom
import warnings
import datetime
import argparse
from numpy import dtype
import netCDF4
from dateutil.relativedelta import relativedelta
warnings.filterwarnings('ignore')
user = getpass.getuser()

def main(region,experiment,date1,date2,workdir):
#    y1 = date1[0:4]
#   m1 = date1[5:7]
#    d1 = date1[8:10]
#
#    y2 = date2[0:4]
#    m2 = date2[5:7]
#    d2 = date2[8:10]
#
    experiment = 'expt_'+experiment[0:2]+'.'+experiment[2]
#    day1=datetime.datetime(int(y1), int(m1), int(d1))-datetime.datetime(int(y1), 1, 1)
#    day2=datetime.datetime(int(y2), int(m2), int(d2))-datetime.datetime(int(y1), 1, 1)

    if region == "TP0":
        region = "TP0a1.00"
    if region == "TP2":
        region = "TP2a0.10"
    if region == "TP4":
        region = "TP4a0.12"
    if region == "TP5":
        region = "TP5a0.06"
    if region == "NAT":
        region = "NATa1.00"
    if region == "NA2":
        region = "NA2a0.80"

    # get domain dimensions
    abgrid = abfile.ABFileGrid(workdir + user + "/" + \
        region + "/topo/regional.grid","r")
    plon=abgrid.read_field("plon")
    plat=abgrid.read_field("plat")
    #ulat=abgrid.read_field("ulat")
    #ulon=abgrid.read_field("ulon")
    #vlat=abgrid.read_field("vlat")
    #vlon=abgrid.read_field("vlon")
    scpx=abgrid.read_field("scpx")
    scpy=abgrid.read_field("scpy")
    jdm,idm=plon.shape

    # get domain depth
    abdepth = abfile.ABFileBathy(workdir + user + "/" + \
        region + "/" + experiment + "/data/regional.depth.b", \
            "r",idm=idm,jdm=jdm)
    depthm=abdepth.read_field("depth")

    nmonths = ( int(date2)+1 - int(date1) ) * 12
    satoutput = np.zeros((nmonths,jdm,idm))
    modeloutput = np.zeros((nmonths,jdm,idm))

    namencout = workdir + user + "/" + \
        region + "/" + experiment + "/data/CHL_monthly_colocated_"+date1+"_"+date2+".nc"
    ncout = NetCDFFile(namencout, "w", format="NETCDF4")

    ncout.createDimension("JJ",jdm)
    ncout.createDimension("II",idm)
    ncout.createDimension("time",nmonths)
#    ncout.createVariable("latitude","f8",("JJ","II"))
#    ncout.createVariable("longitude","f8",("JJ","II"))

    time = ncout.createVariable('time', dtype('double').char, ('time',))
    time.long_name = 'time'
    time.units = 'days since 1970-01-01 00:00'
    time.calendar = 'standard'
    time.axis = 'T'

#    ncout.createVariable("sat_mapped","f8",("time","JJ","II"))
#    ncout.createVariable("model_averaged","f8",("time","JJ","II"))

    sat_mapped = ncout.createVariable('sat_mapped', dtype('double').char, ("time","JJ","II"))
    model_averaged = ncout.createVariable('model_averaged', dtype('double').char, ("time","JJ","II"))
    latitude = ncout.createVariable('latitude', dtype('double').char, ("JJ","II"))
    longitude= ncout.createVariable('longitude', dtype('double').char, ("JJ","II"))

    ncout.variables["longitude"][:]=plon
    ncout.variables["latitude"][:]=plat

    timecal =  datetime.datetime(int(date1), 1, 15)
    for tim in range(nmonths):
        month = int( timecal.strftime('%m') )
        daym = int( timecal.strftime('%d') )
        year = int( timecal.strftime('%Y') )
        print(str(year)+"-"+str(month).zfill(2))
        if region != 'NA2a0.80' :
          f = abfile.ABFileArchv(workdir + user + "/" + \
              region + "/" + experiment + "/data/"+"AVE."+str(month).zfill(2)+"."+str(year)+"."+str(year)+".b","r")
          kdm  = max(f.fieldlevels)
          chl  = np.zeros((kdm,jdm,idm))
          deep = np.zeros((kdm,jdm,idm))
          deepd= np.zeros((kdm,jdm,idm))
          for k in range(kdm) :
              dummy = f.read_field('ECO_diac',k+1)
              chl[k,:,:] = chl[k,:,:] + dummy # first add diatom chl
              dummy = f.read_field('ECO_flac',k+1)
              chl[k,:,:] = chl[k,:,:] + dummy # then add flagellate
              chl[k,:,:] = np.ma.masked_array(chl[k,:,:],dummy.mask)
              dummy = f.read_field('thknss',k+1)/modeltools.hycom.onem
              if k == 0 :
                 deep[k,:,:] = dummy / 2.
                 deepd[k,:,:] = dummy
              else :
                 deepd[k,:,:] = dummy
                 deep[k,:,:] = deepd[k-1,:,:]/2. + deep[k-1,:,:] + dummy / 2.
              deep[k,:,:] = np.ma.masked_array(deep[k,:,:],dummy.mask)

        satdir = 'https://www.oceancolour.org/thredds/dodsC/cci/v5.0-release/geographic/monthly/'
        globchl = satdir+ 'chlor_a/' + str(year) + \
                     '/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-' + \
                            str(year) + str(timecal.month).zfill(2) + '-fv5.0.nc'
        chlname = 'chlor_a'
        #is_it_there=True
        tries = 1000
        for i in range(tries):
         try:
           ncc = NetCDFFile(globchl)
           satlat = ncc.variables['lat'][:]
           satlon = ncc.variables['lon'][:]
    #       time= nc.variables['time'][:]
           satchl = ncc.variables[chlname][0,:,:]
           #satchl = nc.variables['chlor_a'][0,:,:]
           ncc.close()
           break
         except Exception as e :
           if i < tries - 1:
              print("ERROR, try number : ",i+1)
              continue
           else:
              print('File does not exist : chlor_a : ' + \
                 y1 + str(timecal.month).zfill(2) + str(timecal.day).zfill(2))
              #is_it_there=False
              break


#        if is_it_there:
#           ncc = NetCDFFile(globchl)


        globkd = satdir+ 'kd/' + str(year) + \
                    '/ESACCI-OC-L3S-K_490-MERGED-1M_MONTHLY_4km_GEO_PML_KD490_Lee-' + \
                           str(year) + str(timecal.month).zfill(2) + '-fv5.0.nc'
        kdname = 'kd_490'

        for i in range(tries):
           try:
               nck = NetCDFFile(globkd)
               satkd = nck.variables[kdname][0,:,:]
               nck.close()
               break
           except Exception as e :
                if i < tries - 1:
                    print("ERROR, try number : ",i+1)
                    continue
                else:
                    print('File does not exist : chlor_a : ' + \
                       str(year)  + str(timecal.month).zfill(2) + str(timecal.day).zfill(2))
                    satkd = np.copy(satchl)
                    satkd = 1./10.
                    break


    #     satkd[satkd.mask]=1./10. # set 1/kd to 10m if masked. This is effective if chl not masked but kd masked. so chl can compute.

        depthmf = np.asfortranarray(depthm)
        deepf   = np.asfortranarray(deep)
        chlf    = np.asfortranarray(chl)
        satchlf = np.asfortranarray(satchl)
        satkdf   = np.asfortranarray(satkd)
        platf   = np.asfortranarray(plat)
        plonf   = np.asfortranarray(plon)
        satlatf = np.asfortranarray(satlat)
        satlonf = np.asfortranarray(satlon)
        scpxf   = np.asfortranarray(scpx)
        scpyf   = np.asfortranarray(scpy)
 #       maskf   = np.asfortranarray(mask)

#        satout,kdout,modelout = _the_loop.main(depthmf,scpxf,scpyf,platf,plonf,deepf,chlf,satlatf,satlonf,satchlf,satkdf)
        satout,kdout,modelout = _FRAM_the_mapping_loop.main(depthmf,scpxf,scpyf,platf,plonf,deepf,chlf,satlatf,satlonf,satchlf,satkdf)

        satout = np.ma.masked_where(satout>1000.,satout)
        satout = np.ma.masked_where(satout<0.,satout)
        satout = np.ma.masked_array(satout,dummy.mask)
        modelout = np.ma.masked_where(modelout>1000.,modelout)
        modelout = np.ma.masked_where(modelout<0.,modelout)
        modelout = np.ma.masked_array(modelout,dummy.mask)

        ncout.variables["sat_mapped"][tim,:,:] = satout
        ncout.variables["model_averaged"][tim,:,:] = modelout
        ncout.variables["time"][tim] = netCDF4.date2num(timecal, units = time.units, calendar = time.calendar)
        ncout.sync()
  
        timecal = timecal + relativedelta(months=1)

    ncout.close()


if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('region',  type=str)
    parser.add_argument('experiment',  type=str)
    parser.add_argument('date1',  type=str, help='e.g 2010')
    parser.add_argument('date2',  type=str, help='e.g 2010')
    parser.add_argument('--workdir', type=str, \
        default="/cluster/work/users/", \
            help="machine specific work directory (above user folder)")
    args = parser.parse_args()

    main(args.region,args.experiment,args.date1,args.date2,workdir=args.workdir)
