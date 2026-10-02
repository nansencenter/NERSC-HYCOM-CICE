import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.axes import Axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
# import plotly.express as px



# def Plot_plotly(dis_grid,depth_grid):
#     x, y = np.meshgrid(dis_grid.ds.longitude, dis_grid.ds.latitude)
#     fig=px.imshow(dis_grid.ds['dis24'],y=dis_grid.ds.latitude,x=dis_grid.ds.longitude)
#     fig.update_yaxes(autorange=True)
#     fig.show()




def Plot_estuaries(dis_grid,depth_grid):
    land_cmap=cmap = ListedColormap(["black"])

    cmap = plt.get_cmap('rainbow')
    cmap.set_under('black') #Color for value less than vmin
    eps=5 #Flux threshold for plot (cell will appear as black if value under eps)


    fig,ax=plt.subplots(figsize=(15,12))
    # ax.set_aspect('equal')
    # x1,y1,x2,y2=dis_grid._GloFAS_grid__Region[dis_grid._GloFAS_grid__region_name]
    x1,y1,x2,y2=[120,60,140,80]

    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)
    x,y=np.meshgrid(dis_grid.ds.longitude, dis_grid.ds.latitude)
    # im=ax.pcolormesh(x,y,ds_dis['land'],cmap=land_cmap,label='GloFAS grid point')
    im = ax.pcolormesh(x,y, dis_grid.ds['dis24'], vmin=eps,cmap=cmap)


    # fig,ax=plt.subplots()
    # plt.tight_layout()
    # x,y=np.meshgrid(dis_grid.ds.longitude, dis_grid.ds.latitude)
    # # im=ax.pcolormesh(x,y,ds_dis['land'],cmap=land_cmap,label='GloFAS grid point')
    # im = ax.pcolormesh(x,y, dis_grid.ds['uparea'],cmap=cmap)
    # cb = plt.colorbar(im, label='uparea')


    TOPAZ_mask=~depth_grid.depth.mask

    plume_topaz_mask=np.where((depth_grid.lon <x2) & (depth_grid.lon > x1) & (depth_grid.lat < y2) & (depth_grid.lat > y1))
    TOPAZ_min_flux,TOPAZ_max_flux=np.min(dis_grid.river_flux_TOPAZ[TOPAZ_mask]),np.max(dis_grid.river_flux_TOPAZ[TOPAZ_mask])
    sc2=ax.scatter(depth_grid.lon[TOPAZ_mask],depth_grid.lat[TOPAZ_mask],s=30,c=dis_grid.river_flux_TOPAZ[TOPAZ_mask],cmap='rainbow',
                   label='TOPAZ depth grid point',vmin=0,vmax=np.max(dis_grid.river_flux_TOPAZ[plume_topaz_mask]))

    # sc2=ax.scatter(depth_grid.lon[TOPAZ_mask],depth_grid.lat[TOPAZ_mask],s=40,c=depth_grid.depth[TOPAZ_mask],cmap='ocean',
    #                label='TOPAZ depth grid point')


    mask_Estuary=(dis_grid.df_dis_Estuary.longitude>x1) & (dis_grid.df_dis_Estuary.longitude<x2) \
         & (dis_grid.df_dis_Estuary.latitude>y1) & (dis_grid.df_dis_Estuary.latitude<y2)
    flux_min=dis_grid.df_dis_Estuary.dis_Estuary[mask_Estuary].min()
    flux_max = dis_grid.df_dis_Estuary.dis_Estuary[mask_Estuary].max()
    sc1=ax.scatter(dis_grid.df_dis_Estuary.longitude,dis_grid.df_dis_Estuary.latitude,
                   c=dis_grid.df_dis_Estuary.dis_Estuary,
                   s=60,cmap=cmap,label='GloFAS Estuary',
                   vmin=flux_min,vmax=flux_max,marker='v',edgecolors='black')



    sc3=ax.scatter(dis_grid.df_dis_Estuary.longitude_glofas_on_topaz,dis_grid.df_dis_Estuary.latitude_glofas_on_topaz,
                   c=dis_grid.df_dis_Estuary.dis_Estuary,
                   s=60,cmap=cmap,label='GloFAS Estuary projected on Topaz grid',
                   vmin=flux_min,vmax=flux_max,marker='s',edgecolors='black')


    ax.set_aspect('equal')
    ax.set_title('Estuaries '+dis_grid._GloFAS_grid__region_name)
    # Print figure and remove wite space.
    cax1 = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    cax2 = fig.add_axes([ax.get_position().x1 + 0.09, ax.get_position().y0, 0.02, ax.get_position().height])

    cba=plt.colorbar(sc1,label='Total estuary flux '+r'$m^3/s$',cax=cax1)
    cbb=plt.colorbar(sc2,label='Propagated flux '+r'$m^3/s$',cax=cax2)
    # cbc=plt.colorbar(sc3,label='climatology flux '+r'$m^3/s$')


    topaz_grid_point = mlines.Line2D([], [], color='green', marker='o', ls='', label='TOPAZ depth grid point')
    GloFAS_estuary = mlines.Line2D([], [], color='red', marker='v', ls='', label='GloFAS estuary')
    GloFAS_estuary_on_topaz_grid = mlines.Line2D([], [], color='red', marker='s', ls='', label='GloFAS estuary on TOPAZ depth grid')

    plt.legend(handles=[topaz_grid_point, GloFAS_estuary,GloFAS_estuary_on_topaz_grid],loc='upper right')

    dpi = 200
    fig.canvas.print_figure("../../fig/Estuaries_"+dis_grid._GloFAS_grid__region_name, bbox_inches='tight', dpi=dpi)
    plt.show()


def Plot_statistics(dis_grid):
    Uparea_Thresholds=np.linspace(1,1e11,5000)
    Nb_estuaries=np.zeros_like(Uparea_Thresholds)
    Flux_sum = np.zeros_like(Uparea_Thresholds)

    for idx,t in enumerate(Uparea_Thresholds):
        dis_grid.find_estuaries(t)

        Nb_estuaries[idx]=dis_grid.df_dis_Estuary.dis_Estuary.count() #Number of estuaries
        Flux_sum[idx]=dis_grid.df_dis_Estuary.dis_Estuary.sum() #Sum of all fluxes from estuaries

    fig1,ax1=plt.subplots()
    color='blue'
    ax1.semilogx(Uparea_Thresholds,Nb_estuaries,'-',label='Number of estuaries',color=color)
    ax1.set_xlabel('uparea threshold '+r'$m^2$')
    ax1.set_ylabel('Number of estuaries',color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2=ax1.twinx()
    color='red'
    ax2.semilogx(Uparea_Thresholds,Flux_sum,'-',label='Sum of estuaries flux',color=color)
    ax2.set_ylabel('Flux'+r'$m^3/s$',color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    plt.show()


def Plot_maps(dis_grid):


    cmap = plt.get_cmap('rainbow')
    cmap.set_under('black') #Color for value less than vmin
    eps=10
    x,y=np.meshgrid(dis_grid.ds.longitude, dis_grid.ds.latitude)

    fig,ax=plt.subplots()
    plt.tight_layout()
    im = ax.pcolormesh(x,y, dis_grid.ds[dis_grid._GloFAS_grid__field_name],cmap=cmap,vmin=10,vmax=1e3)
    cb=plt.colorbar(im,label='flux '+r'$[m^3/s]$')
    ax.set_aspect('equal')

    fig,ax=plt.subplots()
    plt.tight_layout()
    im = ax.pcolormesh(x,y, dis_grid.ds['uparea'],cmap=cmap,vmax=1e10)
    cb = plt.colorbar(im, label='uparea '+r'$[m^2]$')
    ax.set_aspect('equal')

    fig,ax=plt.subplots()
    plt.tight_layout()
    im = ax.pcolormesh(x,y, dis_grid.ds['elevation'],cmap=cmap)
    cb = plt.colorbar(im, label='elevation '+r'$[m]$')
    ax.set_aspect('equal')

    plt.show()




def Plot_arctic(dis_grid,depth_grid,region_name='Arctic'):
    dis_grid.find_estuaries(uparea_threshold=2e9, elevation_threshold=10)
    flux_min=dis_grid.df_dis_Estuary.dis_Estuary.min()
    flux_max = dis_grid.df_dis_Estuary.dis_Estuary.max()

    cmap = plt.get_cmap('rainbow')
    cmap.set_under('black') #Color for value less than vmin
    eps=10
    x,y=np.meshgrid(dis_grid.ds.longitude, dis_grid.ds.latitude)


    fig_titl=region_name
    filename=region_name+'_estuaries_topaz'
    ###
    proj=ccrs.Stereographic(central_latitude=90.0,central_longitude=-40.0)
    pxy = proj.transform_points(ccrs.PlateCarree(), x, y)
    px=pxy[:,:,0]
    py=pxy[:,:,1]
    figure =plt.figure(figsize=(8,8))
    ax=figure.add_subplot(111, projection=proj)
    plt.tight_layout()
    ax.add_feature(cfeature.GSHHSFeature('auto', edgecolor='black', facecolor='none', linewidth=0.1))
    add_edge=0
    if add_edge == 1:
        ax.add_feature(cfeature.GSHHSFeature('auto', edgecolor='darkgray', facecolor='#f2f2f2', linewidth=0.4))
    ax.set_extent((-2400000, 2050000, -2800000, 2500000), crs=ccrs.NorthPolarStereo())
    P=plt.pcolormesh(px,py,dis_grid.ds[dis_grid._GloFAS_grid__field_name],cmap=cmap,shading='auto',vmin=10,vmax=1e3)
    #-----------------------------
    # # Plot the points
    # ax.scatter(lons_v, lats_v, color='red', marker='o', transform=ccrs.PlateCarree(), s=14)

    # sc2=ax.scatter(depth_grid.lon,depth_grid.lat,s=30,c=depth_grid.depth,cmap='ocean', #Plot TOPAZ grid points
    #                label='TOPAZ depth grid point',transform=ccrs.PlateCarree())

    TOPAZ_mask=~depth_grid.depth.mask
    TOPAZ_min_flux,TOPAZ_max_flux=np.min(dis_grid.river_flux_TOPAZ[TOPAZ_mask]),np.max(dis_grid.river_flux_TOPAZ[TOPAZ_mask])
    sc2=ax.scatter(depth_grid.lon[TOPAZ_mask],depth_grid.lat[TOPAZ_mask],s=30,c=dis_grid.river_flux_TOPAZ[TOPAZ_mask],cmap='rainbow', vmin=TOPAZ_min_flux, vmax=TOPAZ_max_flux,
                   label='TOPAZ depth grid point',transform=ccrs.PlateCarree())

    sc1=ax.scatter(dis_grid.df_dis_Estuary.longitude,dis_grid.df_dis_Estuary.latitude,
                   c=dis_grid.df_dis_Estuary.dis_Estuary,
                    cmap=cmap,label='GloFAS Estuary',
                   vmin=flux_min,vmax=flux_max,marker='v',edgecolors='black',
                   transform=ccrs.PlateCarree(), s=30)

    sc3=ax.scatter(dis_grid.df_dis_Estuary.longitude_glofas_on_topaz,dis_grid.df_dis_Estuary.latitude_glofas_on_topaz,
                   c=dis_grid.df_dis_Estuary.dis_Estuary,
                   s=30,cmap=cmap,label='GloFAS Estuary projected on Topaz grid',
                   vmin=flux_min,vmax=flux_max,marker='s',edgecolors='black',transform=ccrs.PlateCarree())




    aspect = 40; pad_fraction = 0.25
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)


    topaz_grid_point = mlines.Line2D([], [], color='green', marker='o', ls='', label='TOPAZ depth grid point')
    GloFAS_estuary = mlines.Line2D([], [], color='red', marker='v', ls='', label='GloFAS estuary')
    GloFAS_estuary_on_topaz_grid = mlines.Line2D([], [], color='red', marker='s', ls='', label='GloFAS estuary on TOPAZ depth grid')

    plt.legend(handles=[topaz_grid_point, GloFAS_estuary,GloFAS_estuary_on_topaz_grid],loc='upper right')


    # Print figure and remove wite space.
    cax1 = figure.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    cax2 = figure.add_axes([ax.get_position().x1 + 0.2, ax.get_position().y0, 0.02, ax.get_position().height])

    cba=plt.colorbar(sc1,label='flux '+r'$m^3/s$',cax=cax1)
    cbb=plt.colorbar(sc2,label='depth [m]',cax=cax2)
    #if clim is not None : P.set_clim(clim)
    #P.set_clim([1.0,1.0])
    ax.set_title(fig_titl)
    #Print figure.
    #ax.gridlines(linewidth = 0.5)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,linewidth = 0.5)
    #gl.top_labels = False
    #gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 20))
    gl.ylocator = mticker.FixedLocator(range(30, 91, 10))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'gray'}
    gl.ylabel_style = {'size': 10, 'color': 'gray'}
    dpi=200

    figure.canvas.print_figure("../fig/Fig_"+filename,bbox_inches='tight',dpi=dpi)
    ax.clear()

    plt.close(figure)
