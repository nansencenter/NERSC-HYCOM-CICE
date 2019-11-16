#import scipy.io.netcdf
import netCDF4 
import numpy
import logging
import re  
import xml.etree.ElementTree
import datetime
import cfunits
import netcdftime
import scipy


# Set up logger
_loglevel=logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False # Dont propagate to parent in hierarchy (determined by "." in __name__)


class FieldReaderError(Exception):
    """Base class for exceptions in this module."""
    pass


class FieldReader(object) :
   def __init__(self,filenametemplate,coord_props={},time_offset=datetime.timedelta(0)) :
      self._filenametemplate = filenametemplate
      self._filename         = None
      self._coord_props      = coord_props
      self._time_offset      = time_offset  


   def file_is_open(self,newfilename) :
      if self._filename is None :
         return False
      elif newfilename == self._filename :
         return True
      else :
         return False

   def find_timestep(self,dt,varname) :
      tmp = self._coordmap[varname]["time"]-dt
      #print self._filename
      #print dt
      #print tmp[0]
      tmp = [ i.days*86400 + i.seconds  for i in tmp]
      #print tmp[0]
      #print numpy.where(numpy.array(tmp)==0)
      I = numpy.where(numpy.array(tmp)==0)[0]
      return I


   def get_grid(self,varname,dt) :
      self.open_if_needed(dt)
      if self._coordrank[varname]["lon"] > self._coordrank[varname]["lat"] :
         return numpy.meshgrid(self._coordvar["lon"],self._coordvar["lat"])
      else:
         return numpy.meshgrid(self._coordvar["lat"],self._coordvar["lon"])


   def get_coords(self,varname,dt) :
      self.open_if_needed(dt)
      return self._coordvar["lon"],self._coordvar["lat"]


   def get_proj4grid(self,varname,dt) :
      raise NotImplementedError,"get_proj4grid not implemented"



   @classmethod
   def get_field_reader(cls,filenametemplate,format,coord_props=None,time_offset=datetime.timedelta(0)) :
      if format == "netcdf" :
         return NetcdfFieldReader(filenametemplate,coord_props=coord_props,time_offset=time_offset)
      else :
         raise FieldReaderError,"Only netcdf supported at the moment"


class NetcdfFieldReader(FieldReader) :
   def __init__(self,filenametemplate,coord_props={},time_offset=datetime.timedelta(0)) :
      super(NetcdfFieldReader,self).__init__(filenametemplate,coord_props=coord_props,time_offset=time_offset)

   def open(self) :
      #self._nc = scipy.io.netcdf.netcdf_file(self._filename,"r")
      #print "open started"
      logger.info("Opening %s"%self._filename)
      self._nc = netCDF4.Dataset(self._filename,"r")

      # Read and map coordinate variables of all input vars
      self._coordvar={}
      self._coordmap={}
      self._coordrank={}
      for varname,var in self._nc.variables.items() :

         #print varname
         self._coordmap[varname] = {}
         self._coordrank[varname] = {}

         for inumber,i in enumerate(var.dimensions) :
            coordvar = self._nc.variables[i]

            # Set coordinate attributes from coord_props if specified
            if "units" not in dir(coordvar) : 
               if i in self._coord_props and "units" in self._coord_props[i] :
                  logger.info("Setting units from explicit coordinate properties for variable %s"%i)
                  unit_string = self._coord_props[i]["units"]
               else  :
                  raise FieldReaderError,"No units specified for variable %s"%i
            else  :
               if i in self._coord_props and "units" in self._coord_props[i] :
                  logger.warning("Overriding units for variable %s"%i)
                  unit_string = self._coord_props[i]["units"]
               else :
                  unit_string = coordvar.units

            unit = cfunits.Units(unit_string)

            if i not in self._coordvar.keys() :
               # Convert to datetime. use netcdftime as handling is better. 
               coordvals = numpy.array(self._nc.variables[i][:])
               if unit.isreftime :

                  if "calendar" not in dir(coordvar) : 
                     if i in self._coord_props and "calendar" in self._coord_props[i] :
                        logger.info("Setting calendar from explicit coordinate properties for variable %s"%i)
                        calendar = self._coord_props[i]["calendar"]
                     else  :
                        raise FieldReaderError,"No calendar specified for variable %s"%i
                  else  :
                     if i in self._coord_props and "calendar" in self._coord_props[i] :
                        logger.warning("Overriding calendar for variable %s"%i)
                        calendar = self._coord_props[i]["calendar"]
                     else :
                        calendar = coordvar.calendar

                  tmp=netcdftime.utime(unit_string,calendar=calendar)
                  self._coordvar["time"]=tmp.num2date(coordvals)

               elif unit.islongitude :
                  self._coordvar["lon"] = coordvar[:]
               elif unit.islatitude :
                  self._coordvar["lat"] = coordvar[:]
               else :
                  raise FieldReaderError,"Dont know how to handle coordinate variable %s"%i

            if unit.isreftime :
               self._coordmap[varname]["time"] = self._coordvar["time"]
               self._coordrank[varname]["time"] = inumber
            elif unit.islongitude :
               self._coordmap[varname]["lon"] = self._coordvar["lon"]
               self._coordrank[varname]["lon"] = inumber
            elif unit.islatitude :
               self._coordmap[varname]["lat"] = self._coordvar["lat"]
               self._coordrank[varname]["lat"] = inumber
            else :
               raise FieldReaderError,"Dont know how to handle coordinate variable %s"%i

         #print self._coordmap[varname].keys()
      #print "open finished"

   def close(self) :
      self._nc.close()
      self._filename=None

   def open_if_needed(self,dt) :
      # Open file if necessary
      tmpdt=dt-self._time_offset
      newfilename=tmpdt.strftime(self._filenametemplate)
      if not self.file_is_open(newfilename) : 
         if self._filename is not None : self._nc.close()
         self._filename = newfilename
         self.open()


   def get_timestep(self,varname,dt) :
      # Open file if necessary
      self.open_if_needed(dt)

      #TODO: Get rank of time variable in field. For now its assumed to be first (which is the normal)

      #Find timestep
      I = self.find_timestep(dt,varname)
      return self._nc.variables[varname][I,:] # Will actually read 3D, 4D etc




# TODO:
# To implement other readers:
# Subclass FieldReader
# on open (or init) : Define self._coordvar (gets coordinate variables in file)
#                     Time coordinate must be defined as datetime
# on open (or init) : Define self._coordmap (links variables to coord variables)
# Implement init, open, get_timestep and close





class ForcingField(object) :

   def __init__(self,name,filenametemplate,varname,unit,format,accumulation_time=None,rootPath=None,
         coord_props={}) :
      self._name             = name             # Variable names known to this module
      self._filenametemplate = filenametemplate # File known to this module
      self._varname          = varname          # Variable name in file
      self._units            = unit
      self._format           = format

      if rootPath is not None :
         self._filenametemplate = self._filenametemplate.replace("[rootPath]",rootPath)

      self._accumulation_time=accumulation_time
      if self._accumulation_time is not None :
         if self._accumulation_time[-1] == "h" :
            self._accumulation_time = datetime.timedelta(hours=int(self._accumulation_time[:-1]))
         else :
            raise AtmosphericForcingError,"time step must be specified in hours (hours + letter 'h')"

      if self._accumulation_time is not None  :
         # Convert from accumulated to flux
         logger.info("Converting accumulated field (varname=%s) to flux"%self._varname)
         self._units        = self._units+ " s**-1"
         tmp= self._accumulation_time
         self._accumulation_scale_factor = 1./(tmp.days*86400. + tmp.seconds)
      else :
         self._accumulation_scale_factor = 1.
         self._accumulation_time=datetime.timedelta(0)
      #print name,accumulation_time


      self._cfunit             = cfunits.Units(units=self._units)
      self._format             = format
      self._fieldreader = FieldReader.get_field_reader(self._filenametemplate,format,coord_props=coord_props,time_offset=self._accumulation_time) 



   def get_timestep(self,dt,unit=None) : 

      outdt=dt

      # Do unit conversion to get correct output unit
      if unit is not None :
         mycfunit = cfunits.Units(unit)

      logger.debug("Varname %s : imposed unit=%s, my unit=%s"%(self._varname,str(mycfunit), self._units))

      #tmp = numpy.squeeze(self._fieldreader.get_timestep(self._varname,dt))
      tmp = numpy.squeeze(self._fieldreader.get_timestep(self._varname,dt))*self._accumulation_scale_factor
      if not self._cfunit.equals(mycfunit) :
         #print "Unit conversion:",self.varname,"unit=",self._cfunit, "targetunit=", mycfunit
         #print "Unit conversion:max=",tmp.max()
         #print self._accumulation_time
         tmp=cfunits.Units.conform(tmp,self._cfunit,mycfunit)
         #print "Unit conversion:max after=",tmp.max()

#Approach 2: Calculate average at this time
      # If this is an accumulated field, we need to get next field and interpolate
      # TODO: Check if this is really necessary
      if self._accumulation_time <> datetime.timedelta(0):
         dt2 = dt + self._accumulation_time
         logger.info("Computing interpolated value for accumulated field %s (Reading additional field at %s)"%(self._varname,str(dt2)))
         tmp2 = numpy.squeeze(self._fieldreader.get_timestep(self._varname,dt2))*self._accumulation_scale_factor
         if not self._cfunit.equals(mycfunit) :
            #print "Unit conversion:",self.varname,"unit=",self._cfunit, "targetunit=", mycfunit
            #print "Unit conversion:max=",tmp.max()
            tmp2=cfunits.Units.conform(tmp2,self._cfunit,mycfunit)
            #print "Unit conversion:max after=",tmp.max()
         tmp = 0.5*(tmp + tmp2)
         outdt = dt


      #TODO - may be optimized out
      # Sets grid and coords explicitly upon read
      logger.debug("NB: Implicit read of coordinates and grid when reading variable %10s at %s"%(self._varname,str(dt)))
      self.get_coords(dt)
      self.get_grid(dt)

      logger.info("Reading name %20s, varname=%20s"%(self._name,self._varname))
      self._data=tmp
      self._time=outdt

   def get_coords(self,dt) : 
      self._coordx,self._coordy =  self._fieldreader.get_coords(self._varname,dt)

   def get_grid(self,dt) : 
      self._gridx ,self._gridy  = self._fieldreader.get_grid(self._varname,dt)


   @property
   def filenametemplate(self) : 
      return self._filenametemplate

   @property
   def varname(self) : 
      return self._varname

   @property
   def units(self) : 
      return self._units

   @property
   def cfunit(self) : 
      return self._cfunit

   @property
   def format(self) : 
      return self._format

   @property
   def data(self) : 
      return self._data

   def data_to_unit(self,newunit):
      return cfunits.Units.conform(self.data ,self._cfunit,newunit)
           
   @property
   def time(self) : 
      return self._time

   @property
   def grid(self) : 
      return self._gridx,self._gridy

   @property
   def coords(self) : 
      return self._coordx,self._coordy

   @property 
   def is_readable(self) :
      return self._fieldreader is not None




class ForcingFieldFromXml(ForcingField) :
   def __init__(self,xml_element,format,rootPath=None,coord_props={}) :
      name             = xml_element.attrib["known_name"]    # Variable names known to this module
      filenametemplate = xml_element.attrib["file"]          # File known to this module
      varname          = xml_element.attrib["varname"]       # Variable name in file
      units            = xml_element.attrib["units"]         # 

      accumulation_time=None
      if "accumulated" in xml_element.attrib.keys() :
         accumulation_time=xml_element.attrib["accumulated"]
      super(ForcingFieldFromXml,self).__init__(name,filenametemplate,varname,units,format,accumulation_time=accumulation_time,rootPath=rootPath,coord_props=coord_props)




class ForcingFieldCopy(ForcingField) :
   def __init__(self,new_name,instance,unit) :
      self._name             = new_name         # Variable names known to this module
      self._filenametemplate = None             # File known to this module
      self._varname          = None             # Variable name in file
      self._format           = None
      self._fieldreader      = None
      self._units            = unit
      self._cfunit             = cfunits.Units(units=self._units)

      # TODO: use ref or copy?
      self._gridx,self._gridy   = instance.grid
      self._coordx,self._coordy = instance.coords
      self._time = instance.time
      self._data = instance.data

   def set_data(self,data) :
      self._data = data

   def get_timestep(self,dt,unit=None) : 
      self._data   = None
      self._time   = None
      self._gridx  = None
      self._gridy  = None
      self._coordx = None
      self._coordy = None

      
