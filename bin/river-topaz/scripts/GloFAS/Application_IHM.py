import pandas as pd

from grids_construction import Grids
from Estuary_editor import Estuary_editor
from interface import Ui
import numpy as np
import xarray as xr

from scipy.signal import convolve2d

from tkinter import *

from matplotlib.widgets import  RectangleSelector
from matplotlib.widgets import Button

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh


class Appli:
    def __init__(self):
        """
        Contructor of the Appli class. It will launch the user interface (Ui class to slect paramters). Once the users click on run, the Grids object will be calculated and the Map Viewer
        will be opened.
        """
        Region = {'Norway': [5, 55, 30, 80], 'Bothnia': [18, 60, 28, 70], 'Vaasa': [21, 28, 59, 66],
                  'Iceland': [-40, -10, 50, 80], 'Arctic': [-180, 50, 180, 90],
                  'Russia':[40,50,180,90], 'Greenland':[-70,55,-10,85],'America':[-170,50,-60,90],
                  'Scandinavia':[0,50,40,90]}
        self.ui=Ui() #We open the user interface where the user can choose the files and the parameters

        #open river data file:
        river_disharge_file = self.ui.p.param('Input files', 'GloFAS daily disharge file').value()
        river_disharge_file=xr.open_dataset(river_disharge_file)


        self.Grids=Grids(topaz_grid_file=self.ui.p.param('Input files', 'TOPAZ grid A file').value(),
                 topaz_depth_file=self.ui.p.param('Input files', 'TOPAZ depth grid A file').value(),
                 river_disharge_file=river_disharge_file,
                 ldd_file=self.ui.p.param('Input files', 'GloFAS LDD file').value(),
                 uparea_file=self.ui.p.param('Input files', 'GloFAS uparea file').value(),
                 elevation_file=self.ui.p.param('Input files', 'GloFAS elevation file').value(),
                 river_correction_file=self.ui.p.param('Input files', 'River correction file').value(),
                 region_list=Region,
                 region_name=self.ui.p.param('Parameters', 'Region').value(),
                 field_name=self.ui.p.param('Parameters', 'field name').value(),
                 uparea_threshold=self.ui.p.param('Parameters', 'uparea threshold [m2]').value(),
                 elevation_thresold=self.ui.p.param('Parameters', 'elevation threshold [m]').value(),
                 correct_GloFAS=self.ui.p.param('Parameters', 'Correct GloFAS data').value(),
                 climatology=self.ui.p.param('Parameters', 'Climatology').value(),
                 month=self.ui.p.param('Parameters', 'Month').value(),
                 input_edit_file=self.ui.p.param('Input files', 'Estuaries edit file').value(),
                 edit_estuaries=self.ui.p.param('Parameters', 'Edit estuaries').value(),
                 lazy_mode=False)


        df_input_edit=pd.read_csv(self.ui.p.param('Input files', 'Estuaries edit file').value())

        self.Grids.df_input_edit=df_input_edit

        #For developpment purposes:
        # Estuary_editor(self.Grids,zone=[18,60,28,70])


        Map_viewer(self.Grids,self.ui) #Opening the map viewer, where the magic happens



class Map_viewer:
    """Class creating the map viewer interface where it is possible to view estuaries and select areas to edit"""

    def __init__(self,Grids,ui):
        self.Grids=Grids
        self.df_edit = pd.DataFrame(
            {'GloFAS_longitude': [], 'GloFAS_latitude': [], 'TOPAZ_idx': [], 'TOPAZ_longitude': [],
             'TOPAZ_latitude': []})  # dataframe containing new modifications


        self.GloFAS_estuary_output_file=ui.p.param('Output files', 'Estuaries file').value()  #Output files name
        self.estuary_edit_output_file=ui.p.param('Output files', 'Modification file').value()

        self.Open_map_viewer()


    def Open_map_viewer(self):
        #----Plot the GloFAS Data
        self.x, self.y = np.meshgrid(self.Grids.dis_grid.ds.longitude, self.Grids.dis_grid.ds.latitude)
        proj = ccrs.Stereographic(central_latitude=90.0, central_longitude=-40.0)
        pxy = proj.transform_points(ccrs.PlateCarree(), self.x, self.y)
        self.px = pxy[:, :, 0]
        self.py = pxy[:, :, 1]

        #We will downsample Z to plot efficiently
        n=2 #downsampling window size, n=None for no downsample
        self.downscaled_x, self.downscaled_y, self.downscaled_dis=\
            self.x[::n,::n],self.y[::n,::n],self.downscale(self.Grids.dis_grid.ds[self.Grids.dis_grid._GloFAS_grid__field_name],n)

        flux_min = self.Grids.dis_grid.df_dis_Estuary.dis_Estuary.min()
        flux_max=1e3
        # flux_max = Grids.dis_grid.df_dis_Estuary.dis_Estuary.max()

        self.fig =plt.figure(figsize=(11,6))
        self.fig.canvas.manager.set_window_title('Map viewer')
        self.ax = self.fig.add_subplot(111, projection=proj)
        plt.tight_layout()
        cmap = plt.get_cmap('rainbow')
        cmap.set_under('black')  # Color for value less than vmin
        self.ax.add_feature(cfeature.GSHHSFeature('auto', edgecolor='black', facecolor='none', linewidth=0.1))
        add_edge = 0
        if add_edge == 1:
            self.ax.add_feature(cfeature.GSHHSFeature('auto', edgecolor='darkgray', facecolor='#f2f2f2', linewidth=0.4))
        self.ax.set_extent((-3300000, 3050000, -3500000, 3550000), crs=ccrs.NorthPolarStereo())

        # im = ax.pcolormesh(px,py,Grids.dis_grid.ds['dis24'])
        # im=ax.contourf(x,y,Grids.dis_grid.ds['dis24'],transform=ccrs.PlateCarree())
        self.im = self.ax.pcolormesh(self.downscaled_x, self.downscaled_y, self.downscaled_dis,transform=ccrs.PlateCarree(),
                           cmap=cmap,vmin=flux_min, vmax=flux_max)

        self.sc_estuaries = self.ax.scatter( self.Grids.dis_grid.df_dis_Estuary.longitude_glofas_on_topaz,
                         self.Grids.dis_grid.df_dis_Estuary.latitude_glofas_on_topaz,
                         c= self.Grids.dis_grid.df_dis_Estuary.dis_Estuary,
                         s=30, cmap=cmap, label='GloFAS Estuary projected on Topaz grid',
                         vmin=flux_min, vmax=flux_max, marker='s', edgecolors='black', transform=ccrs.PlateCarree())


        gl = self.ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=0.5)
        # gl.top_labels = False
        # gl.right_labels = False
        gl.xlocator = mticker.FixedLocator(range(-180, 181, 20))
        gl.ylocator = mticker.FixedLocator(range(30, 91, 10))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'}
        dpi = 200

        self.rect=RectangleSelector(self.ax, self.select_callback,
        useblit=True,
        button=[1, 3],  # disable middle button
        minspanx=5, minspany=5,
        spancoords='pixels',
        interactive=True)

        #Adding button to open editor
        axbutton = self.fig.add_axes([0.85, 0.05, 0.1, 0.075])
        button = Button(axbutton, 'Open editor')
        button.on_clicked(self.open_estuary_editor)

        #Adding a button to save the edits in a txt file
        ax_write_button = self.fig.add_axes([0.05, 0.3, 0.1, 0.075])
        button_write = Button(ax_write_button, 'Save edits')
        button_write.on_clicked(self.write_modications)

        #Adding a button to load the edits from a txt file
        ax_load_button = self.fig.add_axes([0.05, 0.05, 0.1, 0.075])
        button_load = Button(ax_load_button, 'Load edits')
        button_load.on_clicked(self.load_modif_file)


        #Adding a button to save the estuary dataframe in a txt file
        ax_save_button = self.fig.add_axes([0.85, 0.3, 0.1, 0.075])
        button_save = Button(ax_save_button, 'Save estuaries')
        button_save.on_clicked(self.save_estuaries)

        self.fig.canvas.mpl_connect('key_press_event', self.toggle_selector)
        plt.show()


    def get_lon_lat_from_proj(self):
        """
        Finds the minimum and maximum longitudes and latitudes of the points inside the the rectangle selector
        These will be the coordinates of the corners in the estuary editor
        """
        corners=self.rect.corners #Coordinates of the corners of rectangle selector (projected x y coords)
        LON,LAT=np.zeros(4),np.zeros(4)

        for i in range(4):
            p_ll = corners[0][i], corners[1][i] #We select on of the corners
            P_dist = np.sqrt((p_ll[0] - self.px) ** 2 + (p_ll[1] - self.py) ** 2) #We calculate the distance between the corner and the GloFAS grid points (in projected x y coords)
            ind_ll = np.unravel_index(np.argmin(P_dist, axis=None), P_dist.shape) #We select the GloFAS grid point with min distance to corner point
            LON[i],LAT[i]=self.x[ind_ll],self.y[ind_ll] #We find LON and LAT associated with this point

        self.min_lon,self.max_lon=np.min(LON),np.max(LON)
        self.min_lat, self.max_lat = np.min(LAT), np.max(LAT)


    def select_callback(self, eclick, erelease):
        """
        Callback for line selection.

        *eclick* and *erelease* are the press and release events.
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        # print(f"({x1:3.2f}, {y1:3.2f}) --> ({x2:3.2f}, {y2:3.2f})")
        # print(f"The buttons you used were: {eclick.button} {erelease.button}")
        self.get_lon_lat_from_proj()

        print("min_lon,max_lon=",self.min_lon,self.max_lon)
        print("min_lat,max_lat=", self.min_lat, self.max_lat)

    def toggle_selector(self,event):
        print('Key pressed.')
        if event.key == 't':
            selector=self.rect
            name = type(selector).__name__
            if selector.active:
                print(f'{name} deactivated.')
                selector.set_active(False)
            else:
                print(f'{name} activated.')
                selector.set_active(True)

    def open_estuary_editor(self,event):
        """
        Opens the estuary editor when cliking on the button
        If the changes made in the editor are saved, the function will also make the changes in the Dataframe and update the plot
        """
        zone=self.min_lon,self.min_lat,self.max_lon,self.max_lat

        print("Opening the editor:",zone)

        Editor=Estuary_editor(self.Grids,zone)


        print("Editor Closed")
        if Editor.save==True: #If we click on 'Save and exit', then we add the changes made in the editor to the Dataframe
            for i in range(len(Editor.Index)):
                idx=Editor.Index[i]
                Nlon,Nlat=Editor.New_lon[i],Editor.New_lat[i]
                self.Grids.dis_grid.df_dis_Estuary.at[idx, 'latitude_glofas_on_topaz'] = Nlat
                self.Grids.dis_grid.df_dis_Estuary.at[idx, 'longitude_glofas_on_topaz'] = Nlon
                self.Grids.dis_grid.df_dis_Estuary.at[idx, 'idx_glofas_on_topaz'] = Editor.New_idx[i]

                #Add it to the edit dataframe
                new_row={'GloFAS_longitude':self.Grids.dis_grid.df_dis_Estuary.loc[idx].longitude,
                         'GloFAS_latitude':self.Grids.dis_grid.df_dis_Estuary.loc[idx].latitude,
                         'TOPAZ_idx':Editor.New_idx[i],
                         'TOPAZ_longitude':Nlon,'TOPAZ_latitude':Nlat}
                self.df_edit.loc[len(self.df_edit)]=new_row

                #re-propagate the river plume for the modified estuaries
                self.Grids.dis_grid.propagate_river_plume(self.Grids.topaz_depth_grid,river_index_df=idx,re_propagate=True)

            self.update_plot()


    def write_modications(self, event):
        """
        The modifications made to the estuary positions are written if a csv file
        """
        filename = self.estuary_edit_output_file
        self.df_edit.to_csv(filename)
        print("Estuaries edit dataframe saved in " + filename)

    def load_modif_file(self, event):
        """
        The modifications written in the input csv modification file are added in the map viewer and in the Dataframe
        """
        for i in range(len(self.Grids.df_input_edit)):

            # For each line in the file, we look for corresponding point in th estuary dataframe
            GloFAS_lon, GloFAS_lat = self.Grids.df_input_edit.loc[i].GloFAS_longitude, self.Grids.df_input_edit.loc[
                i].GloFAS_latitude
            TOPAZ_longitude, TOPAZ_latitude = self.Grids.df_input_edit.loc[i].TOPAZ_longitude, \
                                              self.Grids.df_input_edit.loc[i].TOPAZ_latitude

            #---This is now done automatically---
            # try:  # If this point is in the estuary dataframe, we make the replacements (changing the point of estuary on topaz grid)
            #     self.Grids.dis_grid.df_dis_Estuary.loc[(self.Grids.dis_grid.df_dis_Estuary.longitude == GloFAS_lon) &
            #                                            (self.Grids.dis_grid.df_dis_Estuary.latitude == GloFAS_lat),
            #                                            ['idx_glofas_on_topaz', 'longitude_glofas_on_topaz',
            #                                             'latitude_glofas_on_topaz']] = \
            #         [self.Grids.df_input_edit.loc[i].TOPAZ_idx,
            #          TOPAZ_longitude, TOPAZ_latitude]
            #
            #
            # except:  # But the point may not be in the estuary dataframe
            #     print("Point at:" + str(GloFAS_lon) + "," + str(GloFAS_lat) + " is not in the list of estuaries")

            # In any case, the loaded edits are also added to edit list (so they will be in output edit file)
            new_row = {'GloFAS_longitude': self.Grids.df_input_edit.loc[i].GloFAS_longitude,
                       'GloFAS_latitude': self.Grids.df_input_edit.loc[i].GloFAS_latitude,
                       'TOPAZ_idx': self.Grids.df_input_edit.loc[i].TOPAZ_idx,
                       'TOPAZ_longitude': TOPAZ_longitude, 'TOPAZ_latitude': TOPAZ_latitude}
            self.df_edit.loc[len(self.df_edit)] = new_row

        print("Edits loaded")
        self.update_plot()  # And of course we update the plot

    def save_estuaries(self, event):
        """
        We save the estuaries position and flux (and all other info) in a csv file
        This file will eventually serve as an input for TOPAZ
        """

        idm,jdm=self.Grids.topaz_grid.idm,self.Grids.topaz_grid.jdm
        self.Grids.dis_grid.df_dis_Estuary.idx_glofas_on_topaz = self.Grids.dis_grid.df_dis_Estuary.idx_glofas_on_topaz.apply(
            np.unravel_index, args=((jdm, idm), 'C'))  #We transform TOPAZ grid index back into a tuple (i,j) for the TOPAZ model

        filename = self.GloFAS_estuary_output_file
        self.Grids.dis_grid.df_dis_Estuary.to_csv(filename)
        print("Estuaries dataframe saved in " + filename)


    def update_plot(self):
        flux_min = self.Grids.dis_grid.df_dis_Estuary.dis_Estuary.min()
        flux_max=1e3
        # flux_max = self.Grids.dis_grid.df_dis_Estuary.dis_Estuary.max()
        cmap = plt.get_cmap('rainbow')
        cmap.set_under('black')  # Color for value less than vmin
        self.sc_estuaries.remove()
        self.sc_estuaries = self.ax.scatter( self.Grids.dis_grid.df_dis_Estuary.longitude_glofas_on_topaz,
                         self.Grids.dis_grid.df_dis_Estuary.latitude_glofas_on_topaz,
                         c= self.Grids.dis_grid.df_dis_Estuary.dis_Estuary,
                         s=30, cmap=cmap, label='GloFAS Estuary projected on Topaz grid',
                         vmin=flux_min, vmax=flux_max, marker='s', edgecolors='black', transform=ccrs.PlateCarree())

        plt.draw()

    def downscale(self,matrix,n=None):
        #We will downsample Z to plot efficiently
        if n==None:
            return matrix

        kernel = np.ones((n, n))
        convolvedZ = convolve2d(matrix, kernel, mode='valid')

        return convolvedZ[::n, ::n] / n










