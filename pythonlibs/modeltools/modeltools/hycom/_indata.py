import cfunits

# Our "standard" variable names, and how they map to hycom names
variable_names = {
      "10u"     : "wndewd",
      "10v"     : "wndnwd",
      "2t"      : "airtmp",
      "vapmix"  : "vapmix",
      "msl"     : "mslprs",
      "tp"      : "precip",
      "ssrd"    : "shwflx",
      "strd"    : "radflx",
      "ssr"     : "nswrad",
      "str"     : "nlwrad",
      "taux"    : "tauewd",
      "tauy"    : "taunwd",
      "wspd"    : "wndspd",
      "2d"      : "dewpt"
      }


# Units that are required by hycom (udunits)
variable_units = {
      "10u"     : "m s**-1",
      "10v"     : "m s**-1",
      "wspd"    : "m s**-1",
      "2t"      : "degree_Celsius",
      "vapmix"  : "kg kg**-1",
      "msl"     : "Pa",
      "tp"      : "m s**-1",
      "ssrd"    : "W m**-2",
      "strd"    : "W m**-2",
      "ssr"     : "W m**-2",
      "str"     : "W m**-2",
      "taux"    : "N m**-2",
      "tauy"    : "N m**-2",
      "2d"      : "K"
      }


variable_limits = { 
      "vapmix" : [0.,None],
      "tp"     : [0.,None]
      }




#atmosphere_variable_cfunit = dict(
#      [(elem[0],cfunits.Units(elem[1])) for elem in atmosphere_variable_units.items() ]
#      )
#
#def variable_cfunit(kn) :
#   return cfunits.Units(atmosphere_variable_units[kn])
#

