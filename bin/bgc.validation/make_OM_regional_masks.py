import numpy as np
import abfile
import warnings
import matplotlib.pyplot as plt
import pickle
import argparse
warnings.filterwarnings('ignore')

# CAGLAR, 11.Nov.2025
# need to be done once to create masks for further use
# usage: python make_OM_regional_masks.py TP5a0.06 /cluster/work/users/cagyum/TP5a0.06/topo/regional.grid

# ---------
def regmask(masks,lon,lat):
    # you need to run this interactively. This is only for visual checking the regions.
    # command:
    # regmask(maskall,plon,plat)

    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    #cmap = plt.get_cmap('Spectral_r')
    cmap = plt.get_cmap('tab20')

    figsize=(8,8)
    ig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.NorthPolarStereo()})
    ax.set_extent((-180, 180, 50, 90), crs=ccrs.PlateCarree())
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    ax.add_feature(cfeature.COASTLINE, edgecolor='black',linewidth=0.5)
    ax.pcolormesh(lon,lat,masks,transform = ccrs.PlateCarree(),cmap=cmap)
# ---------

def main(region,grid):

  # get domain dimensions
  #abgrid = abfile.ABFileGrid("/cluster/work/users/cagyum/" + region + "/topo/regional.grid","r")
  #abgrid = abfile.ABFileGrid("./regional.grid","r")
  abgrid = abfile.ABFileGrid(grid,"r")
  plon=abgrid.read_field("plon")
  plat=abgrid.read_field("plat")

  x=np.copy(plon)
  y=np.copy(plat)

  maskall = np.zeros((y.shape))

  # barents (BRTS)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>16.), x<=58.)] = 1
  brts=masks.copy()
  maskall = maskall + masks*10.0

  # lofoten (LFTB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>7.), x<=16.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<73.), x>-10.), x<=7.)] = 1
  lftb=masks.copy()
  maskall = maskall + masks*1.0

  # norwegian (NRWB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=60 , y<68.), x>-10.), x<=16.)] = 1
  nrwb=masks.copy()
  maskall = maskall + masks*6.0

  # irminger (IBIB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=55 , y<66.), x>-46.), x<=-10.)] = 1
  ibib=masks.copy()
  maskall = maskall + masks*3.0

  # greenland (GSIS)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=73 , y<81.), x>-25.), x<=7.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=66 , y<73.), x>-37.), x<=-10.)] = 1
  gsis=masks.copy()
  maskall = maskall + masks*13.0

  # kara (KARA)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>58.), x<=100.)] = 1
  kara=masks.copy()
  maskall = maskall + masks*12.0

  # eurasia (ERSB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=78 , y<90.), x>100.), x<=140.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=81 , y<90.), x>-25.), x<=100.)] = 1
  ersb=masks.copy()
  maskall = maskall + masks*2.0

  # laptev (LPES)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<78.), x>100.), x<=180.)] = 1
  lpes=masks.copy()
  maskall = maskall + masks*8.0

  # amerasia (ARSB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=81 , y<90.), x>-123.), x<=-25.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=78 , y<90.), x>140.), x<=180.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=76 , y<90.), x>-181.), x<=-123.)] = 1

  arsb=masks.copy()
  maskall = maskall + masks*14.0

  # chukchi (CHBF)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=66 , y<76.), x>-180.), x<=-158.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<76.), x>-158.), x<=-123.)] = 1

  chbf=masks.copy()
  maskall = maskall + masks*4.0

  # bering (BRGS)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=56 , y<66.), x>=-180.), x<=-158.)] = 1

  brgs=masks.copy()
  maskall = maskall + masks*9.0

  # canadian (CACH)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>-123.), x<=-78.)] = 1

  cach=masks.copy()
  maskall = maskall + masks*5.0

  # baffin (BFFB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>-78.), x<=-51.)] = 1

  bffb=masks.copy()
  maskall = maskall + masks*11.0

  # labrador (LBRD)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=55 , y<68.), x>-64.), x<=-46.)] = 1

  lbrd=masks.copy()
  maskall = maskall + masks*7.0


  masks = {'LBRD':lbrd, 'BFFB':bffb, 'CACH':cach,\
         'BRGS':brgs, 'CHBF':chbf, 'ARSB':arsb,\
          'LPES':lpes, 'ERSB':ersb, 'KARA':kara,\
          'GSIS':gsis, 'IBIB':ibib, 'NRWB':nrwb,\
          'LFTB':lftb, 'BRTS':brts}

  outname = region[0:3]
  #f=open('./OMmask'+outname+'.pckl','wb')
  #pickle.dump(masks,f)
  #f.close()

  with open('/nird/projects/NS9481K/BGC.Validation/OMmask'+outname+'.pckl', 'wb') as f:
    pickle.dump(masks, f, protocol=pickle.HIGHEST_PROTOCOL)


# To load the mask from a different script (modify path and filename if necessary in your own copy):
# import pickle
# f=open('./OMmask'+outname+'.pckl','rb')
# masks = pickle.load(f)
# f.close()

# To retrieve the subregion mask, e.g. NRWB:
# masks['NRWB'] is the same array as plon and plat, but with Norwegian Sea points = 1, the rest = 0.

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('region',  type=str, help='e.g. TP5a0.06')
    parser.add_argument('grid',  type=str, help='e.g. /cluster/work/users/cagyum/TP5a0.06/topo/regional.grid')
    args = parser.parse_args()

    main(args.region,args.grid)
