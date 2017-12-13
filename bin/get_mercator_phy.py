#!/bin/env python

# Name: get_mercator.py
# Purpose: download ocean model data from Mercator
# Author: Mostafa Bakhoday-Paskyabi (Mostafa.Bakhoday@nersc.no)
# Created: Sep 2017
# Copyright: (c) NERSC Norway 2017
# Licence:
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.



import os
from calendar import monthrange

# Please download "motu-client" and set its path here 
# The path of the motu-client package & principle python code
motu_client_name = '/Users/mba098/Documents/Mostafa/Software/models/python_packages/motu-client-python-master/src/python/motu-client.py'

http_address= 'http://purl.org/myocean/ontology/service/database#GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS'

mercator_motu='http://nrtcmems.mercator-ocean.fr/mis-gateway-servlet/Motu'

product_name ='global-analysis-forecast-phy-001-024'

usr      = 'YOUR MERCATOR ACCOUNT'
password = 'YOUR MERCATOR PASSWORD'


# Get latitude and longitude of region of interest as [lonmin latmin lonmax latmax]
region = [-180, 30, 179.91667175293, 90];
depth=[0.494, 5727.9171]
StartYear = 2011
StopYear  = 2012

StandardNames = ['sea_surface_height_above_geoid',
                      'sea_ice_thickness',
                      'sea_ice_area_fraction',
                      'eastward_sea_ice_velocity',
                      'northward_sea_ice_velocity',
                      'sea_water_potential_temperature',
                      'sea_water_salinity',
                      'eastward_sea_water_velocity',
                      'northward_sea_water_velocity',
                      'ocean_mixed_layer_thickness_defined_by_sigma_theta',
                      'sea_water_potential_temperature_at_sea_floor']

VarNames     =  ['zos',
                 'sithick',
                 'siconc',
                 'usi',
                 'vsi',
                 'thetao',
                 'so',
                 'uo',
                 'vo',
                 'mlotst',
                 'bottomT']

# Set true for example the first element if you would like to download 'zos'.
# Please pick up your list of variables from "VarNames" and set through in the corresponding
# element in "Flags".
Flags       = [True,
               True,
               True,
               True,
               True,
               True,
               True,
               True,
               True,
               True,
               True]

DownloadString=''
VarString=''
for i in range(1,len(Flags)+1,1):
    if Flags[i-1]==True:
        VarString = VarString + ' -v %s  ' % VarNames[i-1]

VarString = ' -z '+str(depth[0])+' -Z '+str(depth[1])+VarString
reg=' -x '+str(region[0])+' -X '+str(region[2])+' -y '+str(region[1])+' -Y '+str(region[3])

DownloadString = motu_client_name+ ' -u '+usr+' -p '+password+' -m '+mercator_motu+' -s '+http_address+' -d '+product_name+reg


def extract_data(DownloadString,VarString,StartDate,StopDate,SaveFileName):
    extract_data = 'python '+DownloadString+' -t "'+StartDate+' 12:00:00" -T "'+StopDate+ ' 12:00:00" '+VarString+' -o ./ -f '+SaveFileName
    return extract_data


for iyear in range(StartYear,StopYear,1):
    Year = str(iyear)
    if not os.path.exists(opath+str(iyear)):
        os.makedirs(opath+str(iyear))
    for imonth in range(1,13,1):
        ndays = monthrange(iyear,imonth)[1]
        Month = "%02d" % imonth
        for iday in range(1,ndays+1,1):
            Day = "%02d" % iday
            StartDate = Year+'-'+Month+'-'+Day
            StopDate  = Year+'-'+Month+'-'+Day
            print StartDate+'-'+StopDate
            SaveFileName=opath+str(iyear)+'/'+'MERCATOR-PHY-24-'+Year+'-'+Month+'-'+Day+'-12.nc'
            print SaveFileName
            fname = extract_data(DownloadString,VarString,StartDate,StopDate,SaveFileName)
            print '==============================================='
            print '               ',fname
            os.system(fname)
            print '==============================================='

