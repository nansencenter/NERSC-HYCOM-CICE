import numpy as np
from netCDF4 import Dataset as NetCDFFile
import abfile
import abplot
import getpass
#import _the_loop
import _the_monthly_mapping_loop_
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
        region + "/" + experiment + "/data/NPP_monthly_colocated_"+date1+"_"+date2+".nc"
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

    sat_mapped_eppley = ncout.createVariable('sat_mapped_eppley', dtype('double').char, ("time","JJ","II"))
    sat_mapped_cafe = ncout.createVariable('sat_mapped_cafe', dtype('double').char, ("time","JJ","II"))
    sat_mapped_vgpm = ncout.createVariable('sat_mapped_vgpm', dtype('double').char, ("time","JJ","II"))
    sat_mapped_cbpm = ncout.createVariable('sat_mapped_cbpm', dtype('double').char, ("time","JJ","II"))

    model_averaged_eppley = ncout.createVariable('model_averaged_eppley', dtype('double').char, ("time","JJ","II"))
    model_averaged_cafe = ncout.createVariable('model_averaged_cafe', dtype('double').char, ("time","JJ","II"))
    model_averaged_vgpm = ncout.createVariable('model_averaged_vgpm', dtype('double').char, ("time","JJ","II"))
    model_averaged_cbpm = ncout.createVariable('model_averaged_cbpm', dtype('double').char, ("time","JJ","II"))
    model_averaged_full = ncout.createVariable('model_averaged_full', dtype('double').char, ("time","JJ","II"))
    latitude = ncout.createVariable('latitude', dtype('double').char, ("JJ","II"))
    longitude= ncout.createVariable('longitude', dtype('double').char, ("JJ","II"))

    ncout.variables["longitude"][:]=plon
    ncout.variables["latitude"][:]=plat

#    # get mask from a random netcdf
#    ncc = NetCDFFile("/cluster/projects/nn9481k/NPP_Satellite/vgpm.2003.11.nc")
#    satnpp = ncc.variables['npp'][:]

    timecal =  datetime.datetime(int(date1), 1, 15)
    for tim in range(nmonths):
        month = int( timecal.strftime('%m') )
        daym = int( timecal.strftime('%d') )
        year = int( timecal.strftime('%Y') )
        print(str(year)+"-"+str(month).zfill(2))
        fname = workdir + user + "/" + \
              region + "/" + experiment + "/data/"+"AVE."+str(month).zfill(2)+"."+str(year)+"."+str(year)+".b"
        arc  = abplot.openfile(fname)
#          kdm  = max(f.fieldlevels)
        npp  = abplot.integr(arc,'ECO_prim',200) * 60. * 60. * 24.
        dummy = abplot.getvarib(arc,'thknss',1)
#          deep = np.zeros((kdm,jdm,idm))
#          deepd= np.zeros((kdm,jdm,idm))
#
#          for k in range(kdm) :
#              dummy = f.read_field('thknss',k+1)/modeltools.hycom.onem
#              if k == 0 :
#                 deep[k,:,:] = dummy / 2.
#                 deepd[k,:,:] = dummy
#              else :
#                 deepd[k,:,:] = dummy
#                 deep[k,:,:] = deepd[k-1,:,:]/2. + deep[k-1,:,:] + dummy / 2.
#              deep[k,:,:] = np.ma.masked_array(deep[k,:,:],dummy.mask)

        modelout = np.ma.masked_array(npp,depthm.mask)
        ncout.variables["model_averaged_full"][tim,:,:] = modelout # store the full model domain without satellite masks

        for prefix in ["vgpm","eppley","cbpm","cafe"]:

            satdir = '/cluster/projects/nn9481k/BGC.Validation/NPP_Satellite/'
            globnpp = satdir+prefix+'.' + str(year) + \
                     '.'+str(timecal.month).zfill(2)+'.nc'
#        file = SD(globnpp, SDC.READ)
#        satnpp = file.select('npp')[:]
#        with h5py.File(globnpp, 'r') as ncc:
#                print('... vgpm in ' + str(year) + " " + str(month).zfill(2) + " available")
#                satnpp = ncc['npp'][:]
#        ncc = netCDF4.Dataset(globnpp, 'r')
            try:
            # with h5py.File(globnpp, 'r') as ncc:
            #     print('... vgpm in ' + str(year) + " " + str(month).zfill(2) + " available")
            #     satnpp = ncc['npp'][:]

            # #ncc = netCDF4.Dataset(globnpp, 'r')
            # #'... vgpm in '+str(year)+" "+str(month).zfill(2)+" available"
            #     satnpp = ncc.variables['npp'][:]
            #     lat_dim, lon_dim = satnpp.shape

            #     satlat = np.linspace(-90, 90, lat_dim)
            #     satlon = np.linspace(-180, 180, lon_dim)
            # #ncc.close()

                ncc = NetCDFFile(globnpp)
                print(globnpp)
                satnpp = ncc.variables['npp'][:,:]
#            satnpp = np.ma.masked_where(satnpp < 0.0, satnpp)
                lat_dim, lon_dim = satnpp.shape

                satlat = -np.linspace(-90, 90, lat_dim)
                satlon = np.linspace(-180, 180, lon_dim)

                depthmf = np.asfortranarray(depthm)
                satnppf = np.asfortranarray(satnpp)
                platf   = np.asfortranarray(plat)
                plonf   = np.asfortranarray(plon)
                satlatf = np.asfortranarray(satlat)
                satlonf = np.asfortranarray(satlon)
                scpxf   = np.asfortranarray(scpx)
                scpyf   = np.asfortranarray(scpy)

                satout = _the_monthly_mapping_loop_.main(depthmf,scpxf,scpyf,platf,plonf,satlatf,satlonf,satnppf) # copied from sat.map.TS folder
                satout = np.ma.masked_where(satout>100000.,satout)
                satout = np.ma.masked_where(satout<0.,satout)
                satout = np.ma.masked_array(satout,dummy.mask)
        
#            modelout = np.ma.masked_where(modelout>1000.,modelout)
#            modelout = np.ma.masked_where(modelout<0.,modelout)
                modelout = np.ma.masked_array(modelout,satout.mask)

                ncout.variables["sat_mapped_"+prefix][tim,:,:] = satout
            except:
                pass
            #ncout.variables["sat_mapped"][tim,:,:] = fill_value
            ncout.variables["model_averaged_"+prefix][tim,:,:] = modelout
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
