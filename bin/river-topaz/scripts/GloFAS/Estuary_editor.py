import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtWidgets

import numpy as np
import pandas as pd
import xarray as xr
import copy

import matplotlib.pyplot as plt
from pyqtgraph.parametertree import interact, ParameterTree, Parameter
from pyqtgraph.dockarea.Dock import Dock
from pyqtgraph.dockarea.DockArea import DockArea

class Estuary_editor:
    def __init__(self,Grids,zone):
        self.zone=zone

        self.df_dis_Estuary=copy.deepcopy(Grids.dis_grid.df_dis_Estuary)
        self.ds=copy.deepcopy(Grids.dis_grid.ds)
        self.field_name=Grids.dis_grid._GloFAS_grid__field_name

        self.topaz_depth_grid_flux=copy.deepcopy(Grids.dis_grid.river_flux_TOPAZ)
        self.topaz_depth_grid_depth=copy.deepcopy(Grids.topaz_depth_grid.depth)
        self.topaz_depth_grid_lon = copy.deepcopy(Grids.topaz_depth_grid.lon)
        self.topaz_depth_grid_lat = copy.deepcopy(Grids.topaz_depth_grid.lat)
        self.topaz_depth_grid_idx=np.array([i for i in range(len(self.topaz_depth_grid_depth))]) #indices of points from topaz grid

        #Regional mask
        self.apply_regional_mask()


        #Make app
        self.make_app()



    def apply_regional_mask(self, mask=None):
        """
        Masks as specific region
        :param mask: 2 array corresponding to the region of interest
        """
        min_lon, min_lat, max_lon, max_lat = self.zone

        #Mask depth grid (mask for region and mask points with no depth)
        mask = np.where(
            (self.topaz_depth_grid_lon < max_lon) & (self.topaz_depth_grid_lon > min_lon) &
            (self.topaz_depth_grid_lat < max_lat) & (self.topaz_depth_grid_lat > min_lat) & (np.isnan(self.topaz_depth_grid_depth)==False))
        self.topaz_depth_grid_lon = self.topaz_depth_grid_lon[mask]
        self.topaz_depth_grid_lat = self.topaz_depth_grid_lat[mask]
        self.topaz_depth_grid_depth = self.topaz_depth_grid_depth[mask]
        self.topaz_depth_grid_flux=self.topaz_depth_grid_flux[mask]
        self.topaz_depth_grid_idx=self.topaz_depth_grid_idx[mask]

        #mask Dataframe
        cond=(self.df_dis_Estuary.longitude<max_lon) & (self.df_dis_Estuary.longitude>min_lon) & \
        (self.df_dis_Estuary.latitude < max_lat) & (self.df_dis_Estuary.latitude > min_lat)

        self.df_dis_Estuary=self.df_dis_Estuary[cond]

        #Mask dataset
        mask_lon = (self.ds['longitude'] >= min_lon) & (self.ds['longitude'] <= max_lon)
        mask_lat = (self.ds['latitude'] >= min_lat) & (self.ds['latitude'] <= max_lat)

        self.ds = self.ds.where(mask_lon & mask_lat, drop=True)



    def make_app(self):

        #Setup the editor
        self.app = pg.mkQApp("Estuary editor")
        mw = QtWidgets.QMainWindow()
        mw.resize(800, 800)
        area = DockArea()
        d1 = Dock("Settings", size=(1, 1)) #Dock for buttons
        area.addDock(d1,'left')
        d2 = Dock("Editor", size=(500,500)) #dock for plot
        area.addDock(d2,'right')

        view = pg.GraphicsLayoutWidget()  ## GraphicsView with GraphicsLayout inserted by default
        mw.setCentralWidget(area)
        mw.show()
        mw.setWindowTitle('Estuary editor')
        w1 = view.addPlot()
        d2.addWidget(view)



        #Add parameter tree
        self.tree = ParameterTree()

        params=  [
            {'name': 'Confirm change', 'type':'action'},

            {'name': 'Quit', 'type': 'group', 'children': [
                {'name': 'Save and exit', 'type': 'action'},
                {'name': 'Exit', 'type': 'action'},
            ]}
        ]

        p=Parameter.create(name='params', type='group', children=params)

        p.param('Quit', 'Save and exit').sigActivated.connect(self.save)
        p.param('Quit', 'Exit').sigActivated.connect(self.exit)
        p.param('Confirm change').sigActivated.connect(self.edit_estuary)

        self.tree.setMinimumWidth(150)
        self.tree.addParameters(p, showTop=True)

        self.textbox = Parameter.create(name="text", type="text", readonly=True)
        self.textbox.setValue("select an estuary")
        self.tree.addParameters(self.textbox)
        d1.addWidget(self.tree)

        self.action='1' #if action==1, you should select a topaz grid point, if action=='2', you should select an estuary

        ## Make all plots clickable
        self.clickedPen = pg.mkPen('b', width=2)
        self.lastClicked = []



        #----Plot the GloFAS Data
        x, y = np.meshgrid(self.ds.longitude, self.ds.latitude)


        # Colormap
        colorMap = pg.colormap.get("turbo")
        colorMap.setMappingMode('diverging')
        cm_list=[QtGui.QColor(255,255,255)]
        cmap = plt.get_cmap('rainbow')
        for i in range(cmap.N):
            r,g,b,a = cmap(i)
            cm_list.append(QtGui.QColor(int(r*255), int(g*255), int(b*255)))
        cm=pg.ColorMap(None,cm_list)


        #Pcolormesh for GloFAS data
        pg.setConfigOption('imageAxisOrder', 'row-major')  # Switch default order to Row-major

        Z=xr.where(np.isnan(self.ds[self.field_name]),-9999,self.ds[self.field_name])

        img = pg.PColorMeshItem( x,  y, Z[:-1, :-1], colorMap=cm)
        w1.setAspectLocked()
        w1.addItem(img)
        # w1.invertY(True)
        w1.addColorBar(img, colorMap=cm,values=(-2*np.max(Z.values)/cmap.N,np.max(Z.values)))


        #Scatter plot for Topaz grid points
        #Colormap for topaz grid points

        #Color by depth:
        nPts = self.topaz_depth_grid_depth.shape[0]
        # colormap = pg.colormap.get("ocean",source='matplotlib')
        # valueRange = np.linspace(np.min(self.topaz_depth_grid_depth), np.max(self.topaz_depth_grid_depth), num=nPts)
        # colors = colormap.getLookupTable(0, 1.0, nPts=nPts)
        # self.brushes_depth = colors[np.searchsorted(valueRange, self.topaz_depth_grid_depth)]


        #Color by flux
        colormap = pg.colormap.get("viridis",source='matplotlib')
        valueRange = np.linspace(np.min(self.topaz_depth_grid_flux), np.max(self.topaz_depth_grid_flux), num=nPts)
        colors = colormap.getLookupTable(0, 1.0, nPts=nPts)
        self.brushes_flux=colors[np.searchsorted(valueRange, self.topaz_depth_grid_flux)]

        self.s2 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(None))
        self.s2.addPoints(x=self.topaz_depth_grid_lon,
                     y=self.topaz_depth_grid_lat,
                     data=self.topaz_depth_grid_idx,
                     brush=self.brushes_flux)
        w1.addItem(self.s2)
        self.s2.sigClicked.connect(self.clicked_topaz_grid_point)


        #ScatterPlot for the estuaries

        #Colormap for esuaries
        nPts = self.df_dis_Estuary.longitude_glofas_on_topaz.values.shape[0]
        colormap = pg.colormap.get("rainbow",source='matplotlib')
        valueRange = np.linspace( self.df_dis_Estuary['dis24'].min(), self.df_dis_Estuary['dis24'].max(), num=nPts)
        colors = colormap.getLookupTable(0,  1, nPts=nPts)
        self.brushes = colors[np.searchsorted(valueRange, self.df_dis_Estuary['dis24'])]

        self.s1 = pg.ScatterPlotItem(size=20, pen=pg.mkPen(None))
        self.s1.addPoints(x=self.df_dis_Estuary.longitude_glofas_on_topaz.values,
                     y=self.df_dis_Estuary.latitude_glofas_on_topaz.values,
                     data=self.df_dis_Estuary.index.values,
                     brush=self.brushes,symbol='t')
        w1.addItem(self.s1)
        self.s1.sigClicked.connect(self.clicked_estuary)

        #Add buttons


        #List of modifications
        self.Index=[] #List of index in datraframe to modidy
        self.New_lon=[] # new latitude of estuary at given index
        self.New_lat=[] # new longitude of estuary at given index
        self.New_idx=[] #new index of estuary from topaz grid

        pg.exec()


    def edit_estuary(self):
        """
        Changes the position of the estuary inside the Dataframe. This function is called when clicking on "confirm changes"
        after having selected a new estuary position.
        """
        print("Editor opened")

        self.Index.append(self.idx_estuary)

        self.New_lon.append(self.Nlon)
        self.New_lat.append(self.Nlat)
        self.New_idx.append(self.Nidx)

        # We change the information inside the dataframe
        self.df_dis_Estuary.at[self.idx_estuary,'latitude_glofas_on_topaz']=self.Nlat
        self.df_dis_Estuary.at[self.idx_estuary, 'longitude_glofas_on_topaz'] = self.Nlon
        self.df_dis_Estuary.at[self.idx_estuary, 'idx_glofas_on_topaz'] = self.Nidx

        print("point ", self.idx_estuary, "is now at ",self.Nlon,self.Nlat," (index on topaz grid:",self.Nidx,")")
        self.update_plot()

        self.textbox.setValue("Select an estuary"+"\n"+ "Or select save and exit to leave")


    def update_plot(self):
        self.s1.clear()
        self.s1.addPoints(x=self.df_dis_Estuary.longitude_glofas_on_topaz.values,
                     y=self.df_dis_Estuary.latitude_glofas_on_topaz.values,
                     data=self.df_dis_Estuary.index.values,
                     brush=self.brushes,symbol='t')




    def save(self):

        self.app.quit()
        self.save=True


    def exit(self):

        self.app.quit()
        self.save=False



    def clicked_estuary(self,plot,points):
        """
        Called when clicking on an estuary to be moved. The function saves the index of the estuary (index in the dataframe).
        :param plot:
        :param points: Point clicked
        """
        if self.action=='1':  #if action==1, you should select a topaz grid point, if action=='2', you should select an estuary
            for p in self.lastClicked:
                p.resetPen()
            print("clicked estuary with index ",points[0].data()," at position ",points[0].pos())
            for p in points:
                p.setPen(self.clickedPen)
            self.lastClicked = points
            self.choose_estuary=True
            self.idx_estuary=points[0].data()
            self.action='2'
            self.textbox.setValue('The estuary selected is at position: '+str(points[0].pos())+'\n'+'Now select a point from the topaz grid')
        else:
            print('select a topaz grid point')


    def clicked_topaz_grid_point(self,plot,points):
        """
        Called when clicking on a new TOPAZ grid point after having selected an estuary to be moved. The estuary will be moved on the new TOPAZ grid point (if "confirm changes" is clicked)
        :param plot:
        :param points:
        """
        if self.action=='2': #if action==1, you should select a topaz grid point, if action=='2', you should select an estuary
            for p in self.lastClicked:
                p.resetPen()
            print("clicked topaz grid point at position ",points[0].pos())
            for p in points:
                p.setPen(self.clickedPen)
            self.lastClicked = points
            self.choose_grid=True
            self.Nlon, self.Nlat = points[0].pos()
            self.Nidx=points[0].data()
            self.action = '1'
            self.textbox.setValue('The topaz grid point is at position: ' + str(
                points[0].pos()) + '\n' + 'Click on confirm to make the change')
        else:
            print('select an estuary')