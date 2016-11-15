#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime

import datetime
import cfunits
import numpy
import modeltools.tools
import modeltools.forcing.atmosphere
import modeltools.tools
#from mpl_toolkits.basemap import Basemap, shiftgrid
import logging
import abfile


_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False # Dont propagate to parent in hierarchy (determined by "." in __name__)


# Plot scalar field on model grid
def plot_fig(fld,dt,varname,filename) :
   from matplotlib.backends.backend_agg import FigureCanvasAgg
   from matplotlib.figure import Figure
   fig = Figure(figsize=(5,4), dpi=100)
   ax = fig.add_subplot(111)
   canvas = FigureCanvasAgg(fig)
   P=ax.pcolor(fld)
   fig.colorbar(P)
   ax.set_title("%s at %s"%(varname,str(dt)))
   canvas.print_figure(filename)




def atmfor(start,end,af,grid_file="regional.grid",blkdat_file="blkdat.input",plot_diag=False,
      nersc_forcing=False) :

   if nersc_forcing :
      logger.info("Using old NERSC-HYCOM forcing fields")
      
      # Modify names used by hycom
      mynames = dict(modeltools.hycom.variable_names)
      mynames["10u"] = "uwind"
      mynames["10v"] = "vwind"
      mynames["msl"] = "slp"
      mynames["tcc"] = "clouds"
      #mynames["wspd"] = "wndspd"
      mynames["relhum"] = "relhum"
      mynames["taux"] = "tauewd"
      mynames["tauy"] = "taunwd"

      # Modify output units needed
      myunits = dict(modeltools.hycom.variable_units)
      myunits["msl"]  = "hPa"
      myunits["tcc"]  = "1"
      #myunits["wspd"] = "m s-1"
      myunits["relhum"] = "1"

      mylimits = dict(modeltools.hycom.variable_limits)
      mylimits["tcc"] = [0.,1.]
      mylimits["relhum"] = [0.,1.]

      forcingpropertyset = modeltools.forcing.atmosphere.ForcingPropertySet(
            mynames,myunits, mylimits,
            modeltools.forcing.atmosphere.known_vectors)

   # Standard case, use names and units directly from hycom module
   else :

      forcingpropertyset = modeltools.forcing.atmosphere.ForcingPropertySet(
         modeltools.hycom.variable_names,
         modeltools.hycom.variable_units,
         modeltools.hycom.variable_limits,
         modeltools.forcing.atmosphere.known_vectors)

   # Open hycom grid file, read longitude and latitude@
   # TODO: HYCOM-specific
   #za = modeltools.hycom.io.ABFileRegionalGrid(grid_file,"r")
   za = abfile.ABFileGrid(grid_file,"r")
   mlon = za.read_field("plon")
   mlat = za.read_field("plat")
   Nx=mlon.shape[1]
   Ny=mlon.shape[0]
   za.close()

   # parse blkdat to get yearflag
   # TODO: HYCOM-specific
   blkd = modeltools.hycom.BlkdatParser(blkdat_file)
   yrflag = blkd["yrflag"]
   wndflg = blkd["wndflg"]
   lwflag  = blkd["lwflag"]

   # Main loop 
   always_calculate_interpolator = False
   always_calculate_rotator = False
   field_interpolator={}
   vector_rotator={}
   ffiles={}
   dt = start
   while dt <= end :
       
       logger.info("Reading at %s"%str(dt))
       #print af.known_names

       # Read variables
       af.get_timestep(dt)

       # Estimate dependent variable on native grid
       # radflx is downwelling longwave radiation
       # TODO: HYCOM-specific
       if wndflg in [1,2,3] :
           af.calculate_windstress()
           af.calculate_windspeed()
           af.calculate_ustar()

       #  Forcing used by old NERSC-HYCOM
       if nersc_forcing :
          if "relhum" not in af.known_names_explicit : af.calculate_relhum()
       #  Forcing used by new version 
       else :
          if "vapmix" not in af.known_names_explicit : af.calculate_vapmix()
          if "ssrd"   not in af.known_names_explicit : af.calculate_ssrd()
          if lwflag == -1 :
              if "strd"   not in af.known_names_explicit : af.calculate_strd()
          else :
              raise ValueError,"TODO: lwflag<>-1 not supported"

       # Open output files. Dict uses "known name" when mapping to file object
       # TODO: HYCOM-specific
       if dt == start :
          # Open files
          for k,v in forcingpropertyset.items() :
              if k in af.known_names :
                 ffiles[k]=abfile.ABFileForcing(
                       "forcing.%s"%v.name,"w",idm=Nx, jdm=Ny, 
                       cline1=af.name,
                       cline2="%s (%s)"%(v.name,v.cfunit))

       # Interpolation of all fields and unit conversion
       newfld={}
       for kn in [elem for elem in af.known_names if elem in ffiles.keys()] :

          # Read and convert field to units used by HYCOM
          # TODO: HYCOM-specific
          fld=af[kn].data_to_unit(forcingpropertyset[kn].cfunit)

          # Calculate fieldintepolator object
          if kn not in field_interpolator.keys() or always_calculate_interpolator:
             lo,la=af[kn].coords
             field_interpolator[kn]=modeltools.tools.FieldInterpolatorBilinear(lo,la,mlon,mlat)

          #Actual interpolation
          newfld[kn]=field_interpolator[kn].interpolate(fld)

       # Do rotation of u and v components if this the first component of a vector field
       for kn in af.known_names :
           if kn in modeltools.forcing.atmosphere.known_vectors.keys() and kn in ffiles.keys() :
               knu,knv = modeltools.forcing.atmosphere.known_vectors[kn]
               logger.info("Rotating %s,%s "%(knu,knv))

               # Calculate rotateVector object
               if kn not in vector_rotator.keys() or always_calculate_rotator:
                  vector_rotator[kn] = modeltools.tools.rotateVector(mlon,mlat)
               ur,vr=vector_rotator[kn].rotate(newfld[knu],newfld[knv])
               newfld[knu]=ur
               newfld[knv]=vr

       # Apply limits if specified
       for kn in [elem for elem in af.known_names if elem in ffiles.keys()] :
          newfld[kn] = forcingpropertyset[kn].apply_limit(newfld[kn])

       # Loop over open files and write
       # TODO: HYCOM specific
       for kn in ffiles.keys() :
               
           # Variable name used by hycom
           vname=forcingpropertyset[kn].name

           # Write to hycom file
           newdt=af[kn].time
           ord_day,hour,isec=modeltools.hycom.datetime_to_ordinal(newdt,yrflag)
           dtime=modeltools.hycom.dayfor(newdt.year,ord_day,hour,yrflag)
           ffiles[kn].write_field(newfld[kn],newfld[kn],vname,dtime,af.timestep_in_days)

           # Write diagnostics, if requested
           if plot_diag :
              tmp="forcing.%s."%vname 
              tmp=tmp+"%Y%m%d%H.png"
              tmp=newdt.strftime(tmp)
              logger.info( "plotting %s"%tmp)
              plot_fig(newfld[kn],newdt,vname,tmp)

       # Increase time
       dt = dt + af.timestep
       
   # CLose files            
   for kn in ffiles.keys() :
       ffiles[kn].close()
   af=[]





if __name__ == "__main__" :
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)
   parser = argparse.ArgumentParser(description='Prepare HYCOM forcing files from a set of input files')
   parser.add_argument('--plot_diag', action="store_true")
   parser.add_argument('--nersc_forcing', action="store_true")
   parser.add_argument('start_time', action=DateTimeParseAction, help='Start time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   parser.add_argument('end_time',   action=DateTimeParseAction, help='Stop  time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   parser.add_argument('xml_file',   type=str, help='xml file containing definition of forcing dataset(s)')
   parser.add_argument('xml_id',   type=str, help='xml id of forcing dataset to use in xml file')

   args = parser.parse_args()


   # Set up AtmosphericForcing object, which keeps track of data to be read
   af=modeltools.forcing.atmosphere.AtmosphericForcing(args.xml_file,args.xml_id)
   atmfor(args.start_time,args.end_time,af,plot_diag=args.plot_diag,nersc_forcing=args.nersc_forcing)#,gridfile="regional.grid",blkdat_file="blkdat.input")
