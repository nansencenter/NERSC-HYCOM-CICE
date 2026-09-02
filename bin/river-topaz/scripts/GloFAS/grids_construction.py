import xarray as xr
import abfile
from utils.grid_operations import refine_grid
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
# mpl.use('TkAgg')

from Agrid import Agrid
from GloFAS_grid import GloFAS_grid



class Grids:
    """Class that creats an objects containing all the grids (topaz and GloFAS) as well as a Dataframe with the estuaries"""
    def __init__(self,topaz_grid_file,
                 topaz_depth_file,
                 river_disharge_file,
                 ldd_file,
                 uparea_file,
                 elevation_file,
                 river_correction_file,
                 region_list,
                 region_name='Arctic',
                 field_name='dis24',
                 uparea_threshold=1e9,
                 elevation_thresold=10,
                 correct_GloFAS=True,
                 climatology=False,
                 month=None,
                 lazy_mode=False,
                 input_grids=None,
                 input_edit_file=None,
                 edit_estuaries=False,
                 river_name=None,
                 Propagation_cleaning_step=True,
                 rradius = 80000,
                 alongshoreradius = 200000):
        """
        Contructor of the Grids class. It will create the differents grids (Agrids and GloFAS) using associated classes, add the necessary variables int the GloFAS grid, then run the interpolation
        and find the estuaries
        :param topaz_grid_file:
        :param topaz_depth_file:
        :param river_disharge_file:
        :param ldd_file:
        :param uparea_file:
        :param elevation_file:
        :param region_list:
        :param date_time:
        :param region_name:
        :param idm:
        :param jdm:
        :param field_name:
        :param uparea_threshold:
        :param elevation_thresold:
        """

        if lazy_mode==False: #We do complete process:interpolation, finding estuaries, propagation...

            # ---------------------Grid file--------------------------
            self.topaz_grid = Agrid(topaz_grid_file, region_name=region_name, region_list=region_list)

            # -----------------Depth file----------------------
            self.topaz_depth_grid = Agrid(topaz_depth_file, region_name=region_name, region_list=region_list,idm=self.topaz_grid.TOPAZ_idm,jdm=self.topaz_grid.TOPAZ_jdm, type='depth', lon=self.topaz_grid.lon,
                               lat=self.topaz_grid.lat, mask=self.topaz_grid.mask)

            # ------------------------River data--------------------------
            self.dis_grid = GloFAS_grid(river_disharge_file, region_name=region_name, region_list=region_list,field_name=field_name, climatology=climatology,month=month)

            # ------------------------------Local drainage direction, upstream area and elevation-----------------------------------

            self.dis_grid.add_variable(ldd_file, 'ldd', static=True, x_coord='lon',
                                  y_coord='lat')  # Look at NetCDF metadata for name of x_coord and y_coord


            self.dis_grid.add_variable(uparea_file, 'uparea', static=True, x_coord='longitude', y_coord='latitude')

            self.dis_grid.add_variable(elevation_file, 'elevation', static=True, x_coord='longitude', y_coord='latitude')


            #-----------------------------Operations: Interpolation, finding estuaries, correction----------------------------
            self.dis_grid.get_glofas_into_topaz(self.topaz_depth_grid)

            self.dis_grid.find_estuaries(uparea_threshold=uparea_threshold, elevation_threshold=elevation_thresold)

            self.dis_grid.load_modif_file(input_edit_file=input_edit_file,do_edit=edit_estuaries)

            self.dis_grid.correct_largest_russian_rivers(correction_file=river_correction_file,correct_GloFAS=correct_GloFAS,climatology=climatology)

            # self.dis_grid.isolate_river(river_name=river_name)

            self.dis_grid.propagate_river_plume(depth_grid=self.topaz_depth_grid,Propagation_cleaning_step=Propagation_cleaning_step,
                                                rradius = rradius,alongshoreradius = alongshoreradius)

            self.dis_grid.check_total_fluxes()



        #--------------------------------------LAZY MODE where previous grid is used---------------------------

        else: #In lazy mode, we use data from input Grids:

            self.topaz_grid=input_grids.topaz_grid

            self.topaz_depth_grid=input_grids.topaz_depth_grid


            # ------------------------River data--------------------------
            self.dis_grid = GloFAS_grid(river_disharge_file, region_name=region_name, region_list=region_list,field_name=field_name, climatology=climatology,month=month)


            #-----------------------------Operations: Interpolation, finding estuaries, correction----------------------------
            self.dis_grid.get_glofas_into_topaz(self.topaz_depth_grid,lazy_mode=True, input_grids=input_grids)

            self.dis_grid.load_modif_file(input_edit_file=input_edit_file)

            self.dis_grid.find_estuaries(uparea_threshold=uparea_threshold, elevation_threshold=elevation_thresold,lazy_mode=True,input_grids=input_grids)

            self.dis_grid.load_modif_file(input_edit_file=input_edit_file, do_edit=edit_estuaries)

            self.dis_grid.correct_largest_russian_rivers(correction_file=river_correction_file,correct_GloFAS=correct_GloFAS,climatology=climatology)


            self.dis_grid.propagate_river_plume_lazily(depth_grid=self.topaz_depth_grid,input_grids=input_grids)

            self.dis_grid.check_total_fluxes()

