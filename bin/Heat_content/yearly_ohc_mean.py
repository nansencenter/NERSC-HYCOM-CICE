import argparse
import configparser
import numpy as np
import numpy.ma as ma
import datetime
from datetime import date
import netCDF4
import os
from glob import glob
import xarray as xr

from cal_ohc import ohc
from diag_errinfo import errinfo 
from diag_errinfo import write_onefld 
from plot_ohc import plot_2Dfld
from plot_ohc import plot_myline
from plot_ohc import plot_2lines
from plot_ohc import plot_3lines
from plot_ohc import plot_3lines_anom

#from ngallery_utils import DATASETS

def fld_extraNC(filename,fldname,Nrec):
   import netCDF4 as nc
   ds=nc.Dataset(filename,'r')
   varall=ds.variables.keys()
   if fldname in varall:
      print(fldname)
      if Nrec>=0:
         fld=ds.variables[fldname][Nrec]
      else:
         fld=ds.variables[fldname][:]
      #print(np.shape(fld))
      return fld 
   else:
      #print(varall)
      print("No data ...")
      print('')


# load grid information
# regional.grid.a for the grid cell area
def load_TP2grid(Modpath):
   import abfile.abfile as abf
   Fgrid=Modpath+"regional.grid.a"
   afile=abf.AFile(400,380,Fgrid,"r")
   px=afile.read_record(9)
   py=afile.read_record(10)

   Fdepth=Modpath+"regional.depth.a"
   afile=abf.AFile(400,380,Fdepth,"r")
   pdep=afile.read_record(0)

   return px,py,pdep 


def load_run(yr,filetype,sourdir):
   if filetype=='yearly':
      tempfile=sourdir+"TP2thetao_archv%i.nc"%yr
      salnfile=sourdir+"TP2so_archv%i.nc"%yr
   else:
      print("Under developing...")
      quit()
   dsT=xr.open_dataset(tempfile)
   Rtemp=dsT['thetao'][0,:,:,:]
   dsS=xr.open_dataset(salnfile)
   Rsaln=dsS['so'][0,:,:,:]
   Rdep=dsT['depth'][:]
   return Rtemp,Rsaln,Rdep

cfg=configparser.ConfigParser()
cfg_file=os.path.join(os.path.dirname(os.path.abspath(__file__)),'config.ini')
if not os.path.exists(cfg_file):
   print("Error: config.ini not found at %s"%cfg_file)
   quit()
cfg.read(cfg_file)

TopoDir  =cfg['paths']['topo_dir']
wrkdrt   =cfg['paths']['work_dir']
GlorysDir=cfg['paths']['glorys_dir']
NewrunDir=cfg['paths']['new_run_dir']
OldrunDir=cfg['paths']['old_run_dir']

def load_nemoTP2(yr,filetype,sourdir):
   if filetype=='yearly':
      tempfile=sourdir+"cmems_mod_glo_phy_my_0.083deg_P1M-m_thetao_180.00W-179.92E_30.00N-90.00N_0.49-5727.92m_%i_mean_TP2.nc"%yr
      salnfile=sourdir+"cmems_mod_glo_phy_my_0.083deg_P1M-m_so_180.00W-179.92E_30.00N-90.00N_0.49-5727.92m_%i_mean_TP2.nc"%yr
   else:
      print("Under developing...")
      quit()
   dsT=xr.open_dataset(tempfile)
   #print(dsT['thetao'].shape)
   Gtemp=dsT['thetao'][0,:,:,:]
   dsS=xr.open_dataset(salnfile)
   Gsaln=dsS['so'][0,:,:,:]
   Gdep=dsT['depth'][:]
   return Gtemp,Gsaln,Gdep


parser=argparse.ArgumentParser(description='Calculate Ocean Heat Content')
parser.add_argument('--region', type=int, default=0,
                    help='Region index to compute OHC for (3-8 from NAtlantic_Arctic_regions file). '
                         'Regions 1 and 2 are outside the TOPAZ domain. Default 0 = whole domain.')
parser.add_argument('--depth', type=int, default=0,
                    help='Depth range to compute OHC for: '
                         '1=total, 2=0-300m, 3=0-700m, 4=0-2000m, 5=700-2000m, 6=2000m_plus. '
                         'Default 0 = all depth ranges.')
parser.add_argument('--anomaly', action='store_true',
                    help='Also plot anomaly time series (deviation from each dataset\'s temporal mean).')
args=parser.parse_args()

region_file=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'NAtlantic_Arctic_regions_generic12_TP2.nc')
if args.region in (1, 2):
   print("Error: regions 1 and 2 are outside the TOPAZ domain. Use regions 3-8.")
   quit()
elif args.region > 0:
   ds_reg=netCDF4.Dataset(region_file,'r')
   regionMask=(ds_reg.variables['regions'][:]==args.region)
   ds_reg.close()
   region_tag='-reg%i'%args.region
   print("Applying region mask: region %i"%args.region)
else:
   regionMask=None
   region_tag=''
   print("Computing for whole domain")

Px,Py,Pdepth=load_TP2grid(TopoDir)
Parea0=Px*Py       # unit: m2

depth_ranges = ['total', '0-300m', '0-700m', '0-2000m','700-2000m','2000m_plus']

depth_dict = { 1: 'total', 2: '0-300m', 3: '0-700m',
               4: '0-2000m', 5: '700-2000m', 6: '2000m_plus'}

axes_texts ={'full': 'Full depth', '0-300': '0 - 300 m', '0-700': '0 - 700 m',
              '0-20': '0 - 2000 m', '7-20': '700 - 2000 m', '2plus':  '> 2000 m',}

#filenc='/nird/datapeak/NS9481K/SEACLIM/GLORYS_annual/cmems_mod_glo_phy_my_0.083deg_P1M-m_thetao_180.00W-179.92E_30.00N-90.00N_0.49-5727.92m_1993_mean_TP2.nc'




# NB: the depth levels are different for these two model products
depth_keys=([args.depth] if args.depth in depth_dict else list(depth_dict.keys()))
for isub in depth_keys:
   #if isub<5: 
   #   continue
   #print(depth_dict[2])
   match isub:
      case 1:
          levmin=0.; levmax=4000.
          
      case 2:
          levmin=0.; levmax=300.
      case 3:
          levmin=0.; levmax=700.
      case 4:
          levmin=0.; levmax=2000.
      case 5:
          levmin=700.; levmax=2000.
      case _:
          levmin=2000.; levmax=4000.
   depthMask=(Pdepth>levmin)&(Pdepth<10000.)
   if regionMask is not None:
      maskMod=np.where(depthMask&regionMask)
   else:
      maskMod=np.where(depthMask)
   # global valid water dots
   Parea=Parea0[maskMod]
   Pdeep=Pdepth[maskMod]
   NG=np.size(Pdeep)
   Atotal=np.nansum(Parea)   # total ocean area for this mask [m2]
 
   Osurf=depth_dict[isub]
   # yearly cycle
   Years=[]
   Gmean=[]
   Nmean=[]
   Omean=[]
   for iyr in range(1993,2025,1):
      print("preprocessing %d ..."%iyr)
      Years.append(iyr)
      FileG=wrkdrt+"GOHC_%i-%s%s.nc"%(iyr,Osurf,region_tag)
      FileN=wrkdrt+"NOHC_%i-%s%s.nc"%(iyr,Osurf,region_tag)
      FileO=wrkdrt+"OOHC_%i-%s%s.nc"%(iyr,Osurf,region_tag)

      print("")

      # --- GLORYS ---
      if os.path.exists(FileG):
         print("extract GLORYS OHC from cached file")
         tmpG=fld_extraNC(FileG,"ohc_G%s"%Osurf,-1)
         Gohc=tmpG[maskMod]
      else:
         GT,GS,Glev=load_nemoTP2(iyr,'yearly',GlorysDir)
         Ngk=GT.shape[0]
         print("load data from GLORYS ")
         tmpGT=GT.values[:,maskMod[0],maskMod[1]]
         tmpGS=GS.values[:,maskMod[0],maskMod[1]]
         Gohc=ohc(Ngk,NG,tmpGT,tmpGS,Glev.values,levmin,levmax,Parea,Pdeep)
         tmp0=np.zeros(Pdepth.shape,dtype=float)
         tmp0[maskMod]=Gohc
         tmp0[np.where(tmp0==0.)]=np.nan
         write_onefld(FileG,tmp0,"ohc_G%s"%Osurf,"J","Ocean heat content in %s"%Osurf,"Derived from GLORYS")

      # --- New_ref_run ---
      if os.path.exists(FileN):
         print("extract New_ref_run OHC from cached file")
         tmpN=fld_extraNC(FileN,"ohc_N%s"%Osurf,-1)
         Nohc=tmpN[maskMod]
      else:
         NT,NS,Nlev=load_run(iyr,'yearly',NewrunDir)
         Nnk=NT.shape[0]
         print("load data from New_ref_run ")
         tmpNT=NT.values[:,maskMod[0],maskMod[1]]
         tmpNS=NS.values[:,maskMod[0],maskMod[1]]
         Nohc=ohc(Nnk,NG,tmpNT,tmpNS,Nlev.values,levmin,levmax,Parea,Pdeep)
         tmp0=np.zeros(Pdepth.shape,dtype=float)
         tmp0[maskMod]=Nohc
         tmp0[np.where(tmp0==0.)]=np.nan
         write_onefld(FileN,tmp0,"ohc_N%s"%Osurf,"J","Ocean heat content in %s"%Osurf,"Derived from New_ref_run")

      # --- Old_ref_run ---
      if os.path.exists(FileO):
         print("extract Old_ref_run OHC from cached file")
         tmpO=fld_extraNC(FileO,"ohc_O%s"%Osurf,-1)
         Oohc=tmpO[maskMod]
      else:
         OT,OS,Olev=load_run(iyr,'yearly',OldrunDir)
         Nok=OT.shape[0]
         print("load data from Old_ref_run ")
         tmpOT=OT.values[:,maskMod[0],maskMod[1]]
         tmpOS=OS.values[:,maskMod[0],maskMod[1]]
         Oohc=ohc(Nok,NG,tmpOT,tmpOS,Olev.values,levmin,levmax,Parea,Pdeep)
         tmp0=np.zeros(Pdepth.shape,dtype=float)
         tmp0[maskMod]=Oohc
         tmp0[np.where(tmp0==0.)]=np.nan
         write_onefld(FileO,tmp0,"ohc_O%s"%Osurf,"J","Ocean heat content in %s"%Osurf,"Derived from Old_ref_run")

      Gmean.append(np.nansum(Gohc)/Atotal)
      Nmean.append(np.nansum(Nohc)/Atotal)
      Omean.append(np.nansum(Oohc)/Atotal)
      print("")

   print("Plotting the time series for %s"%depth_dict[isub])
   ytext="yearly_"+depth_dict[isub]+region_tag
   print(ytext)
   Pxx=[]
   for ik in list(Years):
      Pxx.append(date(ik,1,1))

   plot_3lines(Pxx,Gmean,Nmean,Omean,"GLORYS","New_ref_run","Old_ref_run",[],ytext)

   if args.anomaly:
      Ganom=np.array(Gmean)-np.nanmean(Gmean)
      Nanom=np.array(Nmean)-np.nanmean(Nmean)
      Onanom=np.array(Omean)-np.nanmean(Omean)
      plot_3lines_anom(Pxx,Ganom,Nanom,Onanom,"GLORYS","New_ref_run","Old_ref_run",ytext)



#   plot_myline(Pxx,Rbias,'','BIAS',ytext)
#   print("bias done!")
#
#   plot_myline(Pxx,Rrmse,'','RMSE',ytext)
#   print("rmse done!")
#   print(Pxx)
#   plot_myline(Pxx,Rcorr,[.95, 1.],'Coeff',ytext)
#   print("coeff done!")
   

#print("done!")

#fig.text(0.03, 0.5, 'OHC (OM)', va='center', rotation='vertical')











