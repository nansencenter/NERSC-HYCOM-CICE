import sys
from grids_construction import Grids

# data path definition
GloFASdata_path='/cluster/projects/nn9481k/GloFAS_data/TOPAZrunoff_data'
if GloFASdata_path not in sys.path:
   sys.path.append(GloFASdata_path)

#ldd file .nc
ldd_file=GloFASdata_path+'/river/auxiliary_files__glofas_v4_0/ldd_glofas_v4_0.nc'

#uparea file .nc
uparea_file=GloFASdata_path+'/river/auxiliary_files__glofas_v4_0/uparea_glofas_v4_0.nc'

#elevation file .nc
elevation_file=GloFASdata_path+'/river/auxiliary_files__glofas_v4_0/elevation_glofas_v4_0.nc'

#river correction .nc
river_correction_file=GloFASdata_path+'/river/correction_files/largest_river_correction_2001_2019.nc'

#edit_file
estuary_edit_file=GloFASdata_path+'/river/GloFAS_edit_files/estuaries_edit.csv'

#-------------------------------PARAMETERS (no need to change anything here)-----------------------------------

climatology=False #True if file is monthly climatology (keep it False)

month=1 #Month of year in model (don't worry about it)

correct_GloFAS=True #Always True

Edit_estuaries=False #True to take into account estuary edit file and modify the positions for certain estuaries

distance_threshold = 300 #[km] Estuaries further than this from the nearest ocean cell of the model grid are dropped


#--------------------------------CONFIG (don't change anything unless you really know what you're doing)------------------------------------------
#To select a region in particular
region_name='Arctic' #Choose region name (no need to change that here)

#Dict--> 'Region_name':[ min_lon,min_lat,max_lon,max_lat]
Region={'Arctic':[-180,50,180,90], 'Russia':[40,50,180,90], 'Greenland':[-70,55,-10,85],
        'America':[-170,50,-60,90],'Scandinavia':[0,50,40,90], 'East':[0,50,180,90],
        'West':[-180,50,0,90]} #Borders of the region


#--------------------------------Building grids-------------------------------------------

def compute_discharge_noUI(river_data,topaz_grid_file,topaz_depth_file,input_grids=None,lazy_mode=False,climatology=False,month=1,correct_GloFAS=True,Edit_estuaries=False,Propagation_cleaning_step=True,rradius = 80000,alongshoreradius = 200000,sigma=3,pointsource=False,geus=None):
    
    # With point sources the radii and sigma are not used, so they are not printed
    if pointsource:
        spread = "Spreading: pointsource (all water in the mouth's cell)"
    else:
        spread = 'rradius='+str(rradius)+"\n"+'alongshoreradius='+str(alongshoreradius)+"\n"+'sigma='+str(sigma)
    print("-----Parameters-----"+"\n"+"Correct GloFAS:"+ str(correct_GloFAS)+"\n"+ "Edit estuary:"+ str(Edit_estuaries)+"\n"+"Cleaning step in propagation:"+str(Propagation_cleaning_step)+"\n"+spread+"\n"+"---------------------------")

    grids=Grids(topaz_grid_file,topaz_depth_file,river_data,ldd_file,uparea_file,elevation_file,
            region_name=region_name,region_list=Region,river_correction_file=river_correction_file,
            distance_threshold=distance_threshold,correct_GloFAS=correct_GloFAS,
            climatology=climatology,month=month,lazy_mode=lazy_mode,input_grids=input_grids,
            input_edit_file=estuary_edit_file, edit_estuaries=Edit_estuaries,
        Propagation_cleaning_step=Propagation_cleaning_step,rradius=rradius, alongshoreradius=alongshoreradius,
        sigma=sigma,pointsource=pointsource,geus=geus)

    return grids

