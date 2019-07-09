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

# USAGE: python ../bin/get_mercator_bio.py 

import os
from calendar import monthrange

# The path of the motu-client package & principle python code
motu_client_name = '/cluster/home/mostafab/NERSC-HYCOM-CICE/motu-client-python-master/src/python/motu-client.py'

#http_address= 'http://purl.org/myocean/ontology/service/dataset-global-analysis-forecast-bio-001-014'

mercator_motu='http://nrtcmems.mercator-ocean.fr/mis-gateway-servlet/Motu'

#http_address='GLOBAL_ANALYSIS_FORECAST_BIO_001_014-FILE'
http_address='GLOBAL_ANALYSIS_FORECAST_BIO_001_014-TDS'

product_name ='dataset-global-analysis-forecast-bio-001-014'

#
#/wrk/project_bhao/nersc_shared/MERCATOR_DATA/scripts> python /wrk/project_bhao/nersc_shared/MOSTAFA/packages/motu-client-python-1.0.3/src/python/motu-client.py -u mbakhodaypaskya -p MostafaCMEMS2017 -m http://nrtcmems.mercator-ocean.fr/mis-gateway-servlet/Motu -s GLOBAL_ANALYSIS_FORECAST_BIO_001_014-TDS -d dataset-global-analysis-forecast-bio-001-014  -x -180 -X 179.5 -y -89 -Y 90 -t "2013-02-02 12:00:00" -T "2013-02-02 12:00:00" -z 0.493 -Z 0.4942 -v NO3 -v PO4 -v Si  -o ./ -f ./MERCATOR-BIO-14-2012-01-03-12.nc

#usr      = 'INSERT YOUR MERCATOR USER NAME'
#password = 'INSERT YOUR MERCATOR PASSWORD'

usr      = ''
password = ''



# output folder path directory
opath='../'
# Get latitude and longitude of region of interest as [lonmin latmin lonmax latmax]
region = [-180, 30, 179.5, 90];
depth=[0.494, 5727.9171]
#depth=[0.494, 8000.00]
StartYear = 2011
StopYear  = 2019


StandardNames = ['mole_concentration_of_nitrate_in_sea_water (NO3)',
                      'net_primary_production_of_biomass_expressed_as_carbon_per_unit_volume_in_sea_water (PP)',
                      'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water (O2)',
                      'mole_concentration_of_phosphate_in_sea_water (PO4)',
                      'mole_concentration_of_silicate_in_sea_water (SI)',
                      'mass_concentration_of_chlorophyll_a_in_sea_water (CHL)',
                      'mole_concentration_of_dissolved_iron_in_sea_water (FE)',
                      'mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water (PHYC)',
                      'cell_thickness ()',
                      'model_level_number_at_sea_floor ()']

VarNames     =  ['NO3',
                 'PP',
                 'O2',
                 'PO4',
                 'Si',
                 'CHL',
                 'FE',
                 'PHYC',
                 ' ',
                 ' ']

Flags       = [True,
               False,
               False,
               True,
               True,
               False,
               False,
               False,
               False]

DownloadString=''
VarString=''
for i in range(1,len(Flags)+1,1):
    if Flags[i-1]==True:
        VarString = VarString + ' -v %s' %VarNames[i-1]

VarString = ' -z '+str(depth[0])+' -Z '+str(depth[1])+VarString
reg=' -x '+str(region[0])+' -X '+str(region[2])+' -y '+str(region[1])+' -Y '+str(region[3])

DownloadString = motu_client_name+ ' -u '+usr+' -p '+password+' -m '+mercator_motu+' -s '+http_address+' -d '+product_name+reg


def extract_data(DownloadString,VarString,StartDate,StopDate,SaveFileName):
    extract_data = 'python '+DownloadString+' -t "'+StartDate+' 12:00:00" -T "'+StopDate+ ' 12:00:00" '+VarString+' -o ./ -f '+SaveFileName
    return extract_data


for iyear in range(StartYear,StopYear,1):
    Year = str(iyear)
    if not os.path.exists(opath+str(iyear)+'/BIO'):
        os.makedirs(opath+str(iyear)+'/BIO')
    for imonth in range(1,13,1):
        ndays = monthrange(iyear,imonth)[1]
        Month = "%02d" % imonth
        for iday in range(1,ndays+1,1):
            Day = "%02d" % iday
            StartDate = Year+'-'+Month+'-'+Day
            StopDate  = Year+'-'+Month+'-'+Day
            print StartDate+'-'+StopDate
            SaveFileName=opath+str(iyear)+'/BIO'+'/'+'MERCATOR-BIO-14-'+Year+'-'+Month+'-'+Day+'-00.nc'
            print SaveFileName
            fname = extract_data(DownloadString,VarString,StartDate,StopDate,SaveFileName)
            print '==============================================='
            print '               ',fname
            os.system(fname)
            print '==============================================='
