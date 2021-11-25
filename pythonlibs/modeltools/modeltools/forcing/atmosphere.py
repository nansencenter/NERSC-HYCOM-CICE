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
import modeltools.tools

# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False # TODO: :Not sure why this is needed...

_stefanb=5.67e-8
_s0=1365.                 # w/m^2  solar constant
_absh2o=0.09              # ---    absorption of water and ozone
_airdns=1.2

# Variable names known by this module
_all_known_names = [
      "10u",
      "10v",
      "2t",
      "2d",
      "msl",
      "ci",
      "tcc",
      "tp",
      "ro",
      "ssrd",
      "strd",
      "ssr",
      "str",
      "taux",
      "tauy",
      "wspd",
      "ustar",
      "sradtop",
      "vapmix",
      "relhum"
      ]

# Units used by internal calculations in this module
_assumed_units = {
      "10u":"m s**-1",
      "10v":"m s**-1",
      "2t":"K",
      "2d":"K",
      "msl":"Pa",
      "ci":"1",
      "tcc":"1",
      "tp":"m s**-1",
      "ro":"1",
      "ssrd":"W m**-2",
      "strd":"W m**-2",
      "ssr":"W m**-2",
      "str":"W m**-2",
      "taux":"N m**-2",
      "tauy":"N m**-2",
      "wspd":"m s**-1",
      "ustar":"m s**-1",
      "sradtop":"W m**-2",
      "vapmix":"kg kg**-1",
      "relhum":"1"
      }



known_vectors = {
      "taux" : ("taux","tauy"),
      "10u"  : ("10u" ,"10v" )
      }

class AtmosphericForcingError(Exception):
    """Base class for exceptions in this module."""
    pass





class AtmosphericForcing(object) :
   # These are the fileds this routine knows about, and can use to calculate new fields

   def __init__(self,configfile,forcing_dataset,rootpath="") :
      self._configfile = configfile
      self._forcing_dataset=forcing_dataset
      self._tree=xml.etree.ElementTree.ElementTree(file=configfile)

      elements = self._tree.findall('forcing_datasets/forcing_dataset[@name="%s"]'%forcing_dataset)
      if len(elements) > 1 : 
         msg = "Could not find unique dataset %s"%forcing_dataset
         raise AtmosphericForcingError(msg)
      if len(elements) == 0 : 
         msg = "Could not find dataset %s"%forcing_dataset
         raise AtmosphericForcingError(msg)

      # Parse the  Self._element attributes 
      self._rootPath = None
      self._element=elements[0]
      if not rootpath :
          self._rootPath = self._element.attrib["rootPath"]
      else :
          self._rootPath = rootpath
      self._timestep = self._element.attrib["timestep"]
      #if "rootPath" in self._element.attrib.keys(): self._rootPath = self._element.attrib["rootPath"]

     
      if self._timestep[-1] == "h" :
         self._timestep = datetime.timedelta(hours=int(self._timestep[:-1]))
      else :
         raise AtmosphericForcingError("time step must be specified in hours (hours + letter 'h')")

      # Get format - only netcdf currently supported 
      if "format" in self._element.attrib.keys(): 
         self._format=self._element.attrib["format"]
         if self._format != "netcdf" :
            raise AtmosphericForcingError("Only netcdf supported at the moment")
      else :
         raise AtmosphericForcingError("Format must be specified")

      # Parse the available fields and create forcingfield class
      elements = self._element.findall('field')
      self._fields   = {} 
      for xml_element in elements :

         # Name must be "known" to this routine
         name             = xml_element.attrib["known_name"]    # Variable names known to this module
         if name not in _all_known_names :
            msg = "Unknown field with name %s"%name
            raise AtmosphericForcingError(msg)

         # We can specify coordinate properties (some times these are wrongly specified or missing).
         # Here we treat the coords as a dict
         tmp = xml_element.findall('coordinate')
         coord_props={}
         for el2 in tmp :
            if "varname" in el2.attrib.keys() :
               coord_props[el2.attrib["varname"]] = dict([(elem[0],elem[1]) for elem in el2.attrib.items() if elem[0] != "varname"])
         #print coord_props
         #print name,coord_props

         self._fields[name] = modeltools.tools.ForcingFieldFromXml(xml_element,self._format,rootPath=self._rootPath,
               coord_props=coord_props)


   def get_timestep(self,dt,varnames=None) :
      for k,v in self._fields.items() :
         if varnames is None or k in varnames :
            #if v.is_readable :
            #   v.get_timestep(dt,unit=_assumed_units[k])
            v.get_timestep(dt,unit=_assumed_units[k])


   def get_grid(self,dt,varnames=None) :
      flddict={}
      for k,v in self._fields.items() :
         if varnames is None or k in varnames :
            v.get_grid(dt)


   def get_coords(self,dt,varnames=None) :
      flddict={}
      for k,v in self._fields.items() :
         if varnames is None or k in varnames :
            v.get_coords(dt)


#   def get_proj4grid(self,dt,varnames=None) :
#      flddict={}
#      for k,v in self._fields.items() :
#         if varnames is None or k in varnames :
#            flddict[k]=v.get_proj4grid(self._varnames[k],dt)
#      return flddict

#   @property 
#   def field(self) : 
#      return self._fields

   def __getitem__(self,i) :
      return self._fields[i]


   @property
   def timestep(self) :
      return self._timestep

   @property
   def timestep_in_days(self) :
      t=self.timestep
      t=t.days+t.seconds/86400.
      return t

   @property
   def varnames(self) :
      return [ elem.varname for elem in self._fields.values ]

   # Lists all field "known" names  (readable from file and calculated)
   @property
   def known_names(self) :
      return self._fields.keys()

   # Lists all fields (readable from file)
   @property
   def known_names_explicit(self) :
      return [elem[0] for elem in self._fields.items() if elem[1].is_readable]


   def calculate_windstress(self) :
      logger.info("Calculating wind stress (taux, tauy) from wind fields")
      if "10u" in self.known_names and "10v" in self.known_names :
         self._fields["taux"] =  modeltools.tools.ForcingFieldCopy("taux",self["10u"],_assumed_units["taux"])
         self._fields["tauy"] =  modeltools.tools.ForcingFieldCopy("tauy",self["10v"],_assumed_units["tauy"])
         tmp1,tmp2 = windstress(self["10u"].data,self["10v"].data)
         self["taux"].set_data(tmp1)
         self["tauy"].set_data(tmp2)
      else :
         raise AtmosphericForcingError("Can not calculate wind stress without 10 meter winds")


   def calculate_windspeed(self) :
      logger.info("Calculating wind speed (wspd) from wind fields")
      if "10u" in self.known_names and "10v" in self.known_names :
         self._fields["wspd"]    = modeltools.tools.ForcingFieldCopy("wspd",self._fields["10u"],_assumed_units["wspd"])
         self["wspd"].set_data(numpy.sqrt(self["10u"].data**2+self["10v"].data**2))
      else :
         raise AtmosphericForcingError("Can not calculate wind speed without 10 meter winds")


   def calculate_ustar(self) :
      logger.info("Calculating ustar from taux, tauy")
      if "taux" in self.known_names and "taux" in self.known_names :
         self._fields["ustar"]    = modeltools.tools.ForcingFieldCopy("ustar",self["taux"],_assumed_units["ustar"])
         self["ustar"].set_data(numpy.sqrt((self["taux"].data**2+self["tauy"].data**2)*1e-3))
      else :
         raise AtmosphericForcingError("Can not calculate wind stress without 10 meter winds")


   def calculate_vapmix(self) :
      logger.info("Calculating vapmix")
      if "msl" in self.known_names and "2d" in self.known_names:
         e = satvap(self["2d"].data) 
         self._fields["vapmix"]    = modeltools.tools.ForcingFieldCopy("vapmix",self["2t"],_assumed_units["vapmix"])
         self["vapmix"].set_data(vapmix(e,self["msl"].data))
      else :
         raise AtmosphericForcingError("Can not calculate wind stress without 10 meter winds")
     
     
   def calculate_strd(self) :
      logger.info("Calculating strd")
      # Calculates downwelling longwave radiation
      if "tcc" in self.known_names and "2t" in self.known_names and "2d" in self.known_names :
         e = satvap(self["2d"].data)
         self._fields["strd"]    = modeltools.tools.ForcingFieldCopy("strd",self["2d"],_assumed_units["strd"])
         self._fields["strd"].set_data(strd_efimova_jacobs(self["2t"].data,e,self["tcc"].data))
      # Estimate from surface parameters
      elif "tcc" in self.known_names and "2t" in self.known_names :
         self._fields["strd"]          = modeltools.tools.ForcingFieldCopy("strd",self["2t"],_assumed_units["strd"])
         self["strd"].set_data(strd_maykut_jacobs(self["2t"].data,e,self["tcc"].data))
      else :
         raise AtmosphericForcingError("Can not calculate STRD")

     
   def calculate_ssrd(self) :
      logger.info("Calculating ssrd")
      # Calculates downwelling shortwave radiation
      #if "ssrd" in self.known_names :
      #   pass
      ## Estimate from surface parameters
      #elif "tcc" in self.known_names  :
      if "tcc" in self.known_names  :
         lo,la= self["tcc"].grid
         srad_top,cosz,cosz_noon =  qsw_et(self["tcc"].time,lo,la)
         ssrd =  qsw_allsky_rosato(srad_top,cosz,cosz_noon,self["tcc"].data) 
         self._fields["ssrd"] = modeltools.tools.ForcingFieldCopy("ssrd",self["tcc"],_assumed_units["ssrd"])
         self["ssrd"].set_data(ssrd)
         self._fields["sradtop"] = modeltools.tools.ForcingFieldCopy("sradtop",self["tcc"],_assumed_units["sradtop"])
         self["sradtop"].set_data(srad_top)
      else :
         raise AtmosphericForcingError("Can not calculate SSRD")


#MOSTAFA: BEGIN

   def calculate_strd_bignami(self) :
      logger.info("Calculating nersc downwelling strd (Bignami 1995)")
      # Calculates downwelling longwave radiation
      if "tcc" in self.known_names and "2t" in self.known_names and "2d" in self.known_names :
         e = satvap(self["2d"].data)
         self._fields["strd"]    = modeltools.tools.ForcingFieldCopy("strd",self["2d"],_assumed_units["strd"])
         self._fields["strd"].set_data(strd_bignami(self["2t"].data,e,self["tcc"].data))
      else :
         raise AtmosphericForcingError("Can not calculate STRD")
   # Here I decompose the net lonwave radiation formulation into two parts: (1) term with no SST effect; (2) term with effects from SST
   # The SST here is from either observation or atmospheric model result.
   def calculate_lwrad_budyko(self) :
      logger.info("Calculating nersc downwelling strd (Budyko 1974)")
      # Calculates downwelling longwave radiation
      if "tcc" in self.known_names and "2t" in self.known_names and "2d" in self.known_names :
         e = satvap(self["2d"].data)
         lo,la= self["tcc"].grid
         self._fields["lwrad"]    = modeltools.tools.ForcingFieldCopy("strd",self["2d"],_assumed_units["strd"])
         self._fields["lwrad"].set_data(lwrad_budyko(la,self["2t"].data,e,self["tcc"].data))
      else :
         raise AtmosphericForcingError("Can not calculate STRD")


   def calculate_lwrad_berliand(self) :
      logger.info("Calculating downwelling lwrad (Berliand (1952)")
      # Calculates downwelling longwave radiation
      if "tcc" in self.known_names and "2t" in self.known_names and "2d" in self.known_names :
         e = satvap(self["2d"].data)
         lo,la= self["tcc"].grid
         self._fields["lwrad"]    = modeltools.tools.ForcingFieldCopy("strd",self["2d"],_assumed_units["strd"])
         self._fields["lwrad"].set_data(lwrad_berliand(self["2t"].data,e,self["tcc"].data))
      else :
         raise AtmosphericForcingError("Can not calculate STRD")

#
#     TODO: more options will be appeared here
#
#MOSTAFA: END


   def calculate_slp(self) :
      logger.info("Calculating slp")
      if "msl" in self.known_names :
         self._fields["slp"]    = modeltools.tools.ForcingFieldCopy("slp",self["msl"],_assumed_units["msl"])
         self["slp"].set_data( self["msl"].data * 1e-2)
      else :
         raise AtmosphericForcingError("Can not calculate slp fields")


   def calculate_relhum(self) :
      logger.info("Calculating relhum")
      if "2t" in self.known_names and "msl" in self.known_names and "2d" in self.known_names:
         e = satvap(self["2t"].data)
         ed= satvap(self["2d"].data)
         #print self["2d"].data.max(),self["2t"].data.max(),self["msl"].data.max()
         self._fields["relhum"]    = modeltools.tools.ForcingFieldCopy("relhum",self["2t"],_assumed_units["relhum"])
         self["relhum"].set_data(relhumid(e,ed,self["msl"].data)/100.)
      else :
         raise AtmosphericForcingError("Can not calculate wind stress without 10 meter winds")


   @property
   def name(self) : return self._forcing_dataset




class ForcingProperty(object) :
   """Class that keeps some metadata about a forcing dataset"""
   def __init__(self,known_name,variable_name,variable_unit,variable_limits,vector_info) :
      self._known_name      = known_name
      self._variable_name   = variable_name
      self._variable_unit   = variable_unit
      self._vector_info     = vector_info
      self._variable_limits = variable_limits

   def apply_limit(self,data) :
      """ Apply specified limits to input data """
      mydata=numpy.copy(data)
      #print self._variable_limits
      if self._variable_limits[0] is not None :
         mydata = numpy.maximum(self._variable_limits[0],mydata)
      if self._variable_limits[1] is not None :
         mydata = numpy.minimum(self._variable_limits[1],mydata)
      #print "jau",mydata
      return mydata

   @property
   def name(self) :
      return self._variable_name

   @property
   def cfunit(self) :
      return cfunits.Units(self._variable_unit)

   @property
   def unit(self) :
      return self._variable_unit


class ForcingPropertySet(object) :
   """Class that keeps some metadata about a dictionary of forcing dataset (uses ForcingProperty)"""
   def __init__(self,variable_names,variable_units,variable_limits,vector_info) :
      self._forcingproperties={}
      for key in variable_names.keys() :
         unit  = variable_units[key]
         vname = variable_names[key]
         vinfo =None
         lims  =[None,None]
         if key in vector_info.keys()     : vinfo = vector_info[key]
         if key in variable_limits.keys() : lims  = variable_limits[key]
         self._forcingproperties[key] = ForcingProperty(key,vname,unit,lims,vinfo)


   def __getitem__(self,kn) : return self._forcingproperties[kn]
   def keys(self) : return self._forcingproperties.keys()
   def items(self) : return self._forcingproperties.items()




      
#MOSTAFA: BEGIN


def lwrad_berliand(tair,e,cc) :
   # downwelling longwave radiation
   #tair : air temperature [K]
   #e    : near surface vapor pressure [Pa]
   #cc   : cloud cover (0-1)

   # Below formula assumes pressure in mBar
   e_mbar = e * 0.01
   emiss=0.97
   # Clear sky downwelling longwave flux from Eimova(1961)
   strd = emiss*_stefanb * tair**4 * (0.39-0.05*numpy.sqrt(e_mbar))

   # Cloud correction by Berliand (1952)
   strd = strd * (1. - 0.6823 * cc*cc)

   return strd


def strd_bignami(tair,e,cc) :
   # downwelling longwave radiation
   #tair : air temperature [K]
   #e    : near surface vapor pressure [Pa]
   #cc   : cloud cover (0-1)

   # Below formula assumes pressure in mBar (1hPa=1 mbar; 1Pas=0.01mbar)
   e_mbar = e * 0.01
   emiss=0.97
   # Clear sky downwelling longwave flux from Bignami(1995)
   strd = _stefanb * tair**4 * (0.653+0.00535*e_mbar)

   # Cloud correction by Berliand (1952)
   strd = strd * (1. + 0.1762 * cc*cc)

   return strd


def lwrad_budyko(lat,tair,e,cc) :
   # downwelling longwave radiation
   #tair : air temperature [K]
   #e    : near surface vapor pressure [Pa]
   #cc   : cloud cover (0-1)

   et=611.*10.**(7.5*(tair-273.16)/(tair-35.86))
#   print numpy.max(tair),numpy.max(et),numpy.max(e)
#   exit(0)
   # Below formula assumes pressure in mBar
   emiss=0.97
   deg2rad = numpy.pi/180.
   # Clear sky downwelling longwave flux from Budyko(1961)
   chi = 0.5+0.246*numpy.abs(lat*deg2rad)
   cc_cliped=numpy.minimum(1.0,numpy.maximum(cc,0.))
   term1=(0.254-4.95e-5*et)*(1.0-chi*cc_cliped**(1.2) )
   strd = emiss*_stefanb * tair**4 * (term1)


   return strd


#MOSTAFA: END


#http://www.nersc.no/biblio/formulation-air-sea-fluxes-esop2-version-micom
#http://regclim.met.no/rapport_4/presentation16/presentation16.htm
#def qlwd(tair,plat,cc,td)
#   fqlw  =emiss*stefanb*tair**3
#   fqlwcc=1.-(.5+.246*abs(plat(i,j)*radian))*cc**1.2
#c
#      fqlwi1=fqlw*tair*((.254-4.95e-5*vpair_i)*fqlwcc-4.)
#      fqlwi2=fqlw*4.
#c
#      fqlww1=fqlw*tair*((.254-4.95e-5*vpair_w)*fqlwcc-4.)
#      fqlww2=fqlwi2

def strd_efimova_jacobs(tair,e,cc) :
   #tair : air temperature [K]
   #e    : near surface vapor pressure [Pa]
   #cc   : cloud cover (0-1)

   # Below formula assumes pressure in mBar
   e_mbar = e * 0.01

   # Clear sky downwelling longwave flux from Eimova(1961)
   strd = _stefanb * tair**4 * (0.746+0.0066*e_mbar) 
   
   # Cloud correction by Jacobs(1978)
   strd = strd * (1. + 0.26 * cc) 

   return strd

   
def strd_maykut_jacobs(tair,cc) :
   #tair : air temperature [K]
   #cc   : cloud cover (0-1)

   # Below formula assumes pressure in mBar
   e_mbar = e * 0.01

   # Clear sky downwelling longwave flux from Maykut and Church (1973). 
   strd = _stefanb * tair**4 * 0.7855
   
   # Cloud correction by Jacobs(1978)
   strd = strd * (1. + 0.26 * cc) 

   return strd



def windstress(uwind,vwind) :
   karalight=True
   ws=numpy.sqrt(uwind**2+vwind**2)
   if karalight :
      wndfac=numpy.maximum(2.5,numpy.minimum(32.5,ws))
      cd_new = 1.0E-3*(.692 + .0710*wndfac - .000700*wndfac**2)
   else :
      wndfac=(1.+sign(1.,ws-11.))*.5
      cd_new=(0.49+0.065*ws)*1.0e-3*wndfac+cd*(1.-wndfac)
   wfact=ws*_airdns*cd_new
   taux = uwind*wfact
   tauy = vwind*wfact
   ustar = numpy.sqrt((taux**2+tauy**2)*1e-3)
   return taux, tauy



def vapmix(e,p) :
   # Input is :
   # e = vapour pressure (partial pressure of vapour)
   # p = air pressure
   vapmix = 0.622 * e / (p-e)
   return vapmix


def satvap(t) :
   # This function calculates the saturation vapour pressure
   # [Pa] from the temperature [deg K].
   # Modified: Anita Jacob, June '97
   #
   # Input: t: temperature [deg K]
   # Output: satvap: saturation vapour pressure at temp. t
   #
   # es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual
   #data c1/610.78/,t00/273.16/
   c1=610.78
   t00=273.16

   #if (t < t00) then
   #   c3 = 21.875
   #   c4 = 7.66
   #else
   #   c3 = 17.269
   #   c4 = 35.86
   #endif
   #KAL !!! c3 = numpy.where(t < t00,21.875,7.66)
   #KAL !!! c4 = numpy.where(t < t00,17.269,35.86)

   # Old hycom
   #c3 = numpy.where(t < t00, 21.875,17.269)
   #c4 = numpy.where(t < t00,  7.66, 35.86)

   # From newest IFS (CY41R2)
   c3 = numpy.where(t < t00, 22.587,17.502)
   c4 = numpy.where(t < t00, -0.7  ,32.19)

   aa = c3 * (t - t00)
   bb = t - c4
   cc=aa/bb

   #if (cc < -20.0) then
   #   satvap=0.0
   #else
   #   satvap = c1 * exp(aa/bb)
   satvap=numpy.where(cc<-20.0,0.0,c1 * numpy.exp(aa/bb))
   return satvap

def  relhumid(sva,svd,msl) :
   # This routine calculates the relative humidity by the 
   # dew point temperature and the mean sea level pressure.
   # Modified: Anita Jacob, June '97
   # Input:
   #    sva: saturatn vapour press at air temp [K]
   #    svd: saturatn vapour press at dew pt temp [K]
   #    msl: pressure at mean sea level [Pa]
   # Output: 
   #   relhumid: Relative Humidity

   # We use the Tetens formula:
   # es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual
   #              es(Tdew)        p - es(Tair)
   # RH = 100 *  -----------   *  ------------
   #             p - es(tdew)       es(Tair)
   aaa=msl - svd
   aaa = svd/aaa
   bbb = (msl - sva)/sva
   relhumid = 100. * aaa * bbb
   return relhumid



def qsw_et(dtime,plon,plat) :
   # BAsed on equations in 
   # Fourier series representation of the position of the sun
   # J. W. Spencer
   # CSIRO Division of Building Research
   # Melbourne, Victoria
   #
   # dtime is datetime object
   # plon in degrees
   # plat in degrees
   # NB: Only suitable for "present day climate"

   radian = numpy.pi/180.

   tmp  = dtime-datetime.datetime(dtime.year,1,1,0,0,0)
   tmp2 = datetime.datetime(dtime.year+1,1,1,0,0,0)-datetime.datetime(dtime.year,1,1,0,0,0)
   tmp = tmp.days+tmp.seconds/86400.
   tmp2 = tmp2.days+tmp2.seconds/86400.
   logger.debug("Day of year - 31. dec of year: %.4f %.4f"%(tmp,tmp2))
   dangle = 2*numpy.pi * float(tmp) / tmp2
   hangle = dtime.hour/24. + dtime.minute/(24.*60.) + dtime.second/3600.
   #hangle = numpy.mod(hangle-0.5,1.)*2*numpy.pi   # Solar hour angle, 0 at noon
   hangle = (hangle-0.5)*2*numpy.pi   # Solar hour angle, 0 at noon

   logger.debug("day angle=%.4f"%dangle)
   logger.debug("time hour angle at Greenwich=%.4f"%hangle)

   if abs(dtime.year - 2000) > 3000. :
      raise AtmosphericForcingError("qsw_et only suitable for present day climate")


   # Solar declination in radians
   decli=.006918+.070257*numpy.sin(dangle)   -.399912*numpy.cos(dangle)      \
                +.000907*numpy.sin(2.*dangle)-.006758*numpy.cos(2.*dangle)   \
                +.001480*numpy.sin(3.*dangle)-.002697*numpy.cos(3.*dangle)

   # Equation of time in radians
   eot = 0.0000075 + 0.001868*numpy.cos(dangle) - 0.032077*numpy.sin(dangle) \
        -0.014615*numpy.cos(2*dangle)           - 0.040849*numpy.sin(2*dangle)

   # Inverse square distance from sun (1/r**2)
   isqd =   1.000110 + 0.034221*numpy.cos(dangle)   + 0.001280*numpy.sin(dangle)  \
                     + 0.000719*numpy.cos(2*dangle) + 0.000077*numpy.sin(2*dangle)

   # eot indicates offset from UT time. Negative means that UT time is faster than  solar time
   hangle = hangle + eot

   # Local solar hour angle
   loc_hangle=hangle+plon*radian

   # Solar Zenith angle
   cosz_noon = numpy.sin(plat*radian)*numpy.sin(decli) + numpy.cos(plat*radian)*numpy.cos(decli)
   cosz      = numpy.sin(plat*radian)*numpy.sin(decli) + numpy.cos(plat*radian)*numpy.cos(decli)*numpy.cos(loc_hangle)

   cosz     =numpy.maximum(0.,numpy.minimum(1.,cosz     ))
   cosz_noon=numpy.maximum(0.,numpy.minimum(1.,cosz_noon))

   srad =_s0*cosz

   logger.debug("declination in degrees      : %.4f"%(decli*180./numpy.pi))
   logger.debug("equation of time  in degrees: %.4f"%(eot*180./numpy.pi))
   logger.debug("equation of time  in minutes: %.3f"%(1440 * eot/(2*numpy.pi)))
   #logger.debug("Zenith at solar noon        : ",numpy.arccos(cosz)*180./numpy.pi)
   logger.debug("time hour angle at Greenwich, corrected for eot=%.4f"%hangle)

   return srad,cosz,cosz_noon


def qsw_allsky_rosato(srad_top,cosz,cosz_noon,cc) :
   # Follows Rosato and Miyakoda[1988]
   # srad = cloud-top incident radiation
   # cosz = cosine of solar zenith angle

   # direct component
   sdir=srad_top*0.7**(1./(cosz+1e-2))     #direct radiation component
   sdif=((1.-_absh2o)*srad_top-sdir)*.5        #diffusive radiation component

   # Solar altitude
   altdeg=numpy.maximum(0.,numpy.arcsin(cosz_noon))*180./numpy.pi #solar noon altitude in degrees

   cfac=(1.-0.62*cc+0.0019*altdeg)               #cloudiness correction by Reed(1977)
   ssurf=(sdir+sdif)*cfac

   return ssurf









def qsw0(qswtime,daysinyear,cc,plat,plon) :
   #
   # --- -------------------------------------------------------------------
   # --- compute 24 hrs mean solar irrradiance at the marine surface layer
   # --- (unit: w/m^2)
   # --- -------------------------------------------------------------------
   #
   # --- Average number of days in year over a 400-year cycle (Gregorian Calendar)
   daysinyear400=365.2425
   #c --- set various quantities
   pi2=8.*numpy.arctan(1.)          #        2 times pi
   deg=360./pi2             #        convert from radians to degrees
   rad=pi2/360.             #        convert from degrees to radians
   eepsil=1.e-9             #        small number
   ifrac=24                 #        split each 12 hrs day into ifrac parts
   fraci=1./ifrac           #        1 over ifrac
   absh2o=0.09              # ---    absorption of water and ozone
   s0=1365.                 # w/m^2  solar constant
   radian=rad
#c
#c --- -------------------------------------------------------------------
#c --- compute 24 hrs mean solar radiation at the marine surface layer 
#c --- -------------------------------------------------------------------
#C --- KAL: TODO - adhere to hycom time setup
   day=numpy.mod(qswtime,daysinyear)    #0 < day < 364
   day=numpy.floor(day)
#c
   dangle=pi2*day/float(daysinyear)   #day-number-angle, in radians 
   if day<0. or day>daysinyear+1 :
      print('qsw0: Error in day for day angle')
      print('Day angle is ',day,daysinyear,qswtime)
      raise NameError("test")
      
# --- compute astronomic quantities -- 
   decli=.006918+.070257*numpy.sin(dangle)   -.399912*numpy.cos(dangle)      \
                +.000907*numpy.sin(2.*dangle)-.006758*numpy.cos(2.*dangle)   \
                +.001480*numpy.sin(3.*dangle)-.002697*numpy.cos(3.*dangle)

   sundv=1.00011+.001280*numpy.sin(dangle)   +.034221*numpy.cos(dangle)      \
                +.000077*numpy.sin(2.*dangle)+.000719*numpy.cos(2.*dangle)

# --- compute astronomic quantities

   sin2=numpy.sin(plat*radian)*numpy.sin(decli)
   cos2=numpy.cos(plat*radian)*numpy.cos(decli)
#
# --- split each day into ifrac parts, and compute the solar radiance for 
# --- each part. by assuming symmetry of the irradiance about noon, it
# --- is sufficient to compute the irradiance for the first 12 hrs of
# --- the (24 hrs) day (mean for the first 12 hrs equals then the mean
# --- for the last 12 hrs)
#
# --- TODO - This routine can also return daily varying solar heat flux
   scosz=0.
   stot=0.
   for npart in range(1,25) :
      bioday=day+(npart-.5)*fraci*.5
      biohr=bioday*86400.                #hour of day in seconds
      biohr=numpy.mod(biohr+43200.,86400.)    #hour of day;  biohr=0  at noon
      hangle=pi2*biohr/86400.            #hour angle, in radians
#
      cosz=numpy.maximum(0.,sin2+cos2*numpy.cos(hangle)) #cosine of the zenith angle
      scosz=scosz+cosz                     #  ..accumulated..
      srad =s0*sundv*cosz                  #extraterrestrial radiation
#
#         sdir=srad*0.7**(1./(cosz+eepsil))    #direct radiation component
#         sdir=srad * exp(-0.356674943938732447/(cosz+eepsil))         
# ---    KAL prevent underflow - .7^100 = 3x10^-16 
      sdir=srad*0.7**(numpy.minimum(100.,1./(cosz+eepsil)))    #direct radiation component
#
      sdif=((1.-absh2o)*srad-sdir)*.5               #diffusive radiation component
      altdeg=numpy.maximum(0.,numpy.arcsin(numpy.minimum(1.0,sin2+cos2)))*deg #solar noon altitude in degrees
      cfac=(1.-0.62*cc+0.0019*altdeg)               #cloudiness correction 
      ssurf=(sdir+sdif)*cfac
      stot=stot+ssurf

#     enddo
   scosz=scosz*fraci               #24-hrs mean of  cosz
   radfl0=stot*fraci               #24-hrs mean shortw rad in w/m^2
#
# --  Original formula was wrong ...
#     !cawdir(i,j)=1.-numpy.maximum(0.15,0.05/(scosz+0.15)) 
#     !cawdir(i,j)=1.-numpy.maximum(0.03,0.05/(scosz+0.15))  !Correction   - Mats
   cawdir=1.-numpy.minimum(0.15,0.05/(scosz+0.15))   #Correction 2 - KAL
#     enddo
#     enddo
#     enddo
#$OMP END PARALLEL DO
#
#     end subroutine qsw0

   return radfl0,cawdir





