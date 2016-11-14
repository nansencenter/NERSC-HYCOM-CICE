## test_pivots.py
## Author: Timothy Williams
## Date:   20151111
## test the pivots obtained from get_pivots.sh etc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

# routines from $SWARP_ROUTINES/py_funs
# - add this dir to PYTHONPATH or use sys.path.append
import mod_reading as Mr
import fns_plotting as Fplt

SHOW  = 0
odir  = 'out'
if not os.path.exists(odir):
   os.mkdir(odir)

###############################################################################################
##work out filename for test date, and get record number (irec) if necessary:
d_path   = '/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC/2015/analysis_m2/'
ncfil    = d_path+'/SWARP_WW3_ARCTIC-12K_20151109.nc' #full name of input file with data lon/lat 
ncfil2   = ncfil                                      #full name of input file with test variable
irec     = 0

var_name    = 'hs';
dname       = 'ww3_arctic'
piv_path    = '/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC/' # path of pivot point binaries
M_INP       = '/work/shared/nersc/msc/ModelInput/'

test_xtra_vars = 1
if test_xtra_vars:
   var_name2   = ['fp','dir']
   Nv          = len(var_name2)
   ncfil3      = []
   for j in range(Nv):#check values of other variables at (itest,jtest)
      # names of their input data files
      ncfil3.append(ncfil)
###############################################################################################

###############################################################################################
print('\nReading '+ncfil)
nci         = Mr.nc_getinfo(ncfil)
VBL         = nci.get_var(var_name,time_index=irec)   # get original variable
dlon,dlat   = nci.get_lonlat()
Vmin0       = VBL.values.min()
Vmax0       = VBL.values.max()

###############################################################################################


# loop over model domains:
# for MOD in ['TP4']:
# for MOD in ['BS1']:
# for MOD in ['FR1']:
for MOD in ['TP4','BS1','FR1']:

   print('\n*********************************************************')
   print('Testing pivot points for '+MOD+' model domain')
   print('*********************************************************\n')

   ###############################################################################################
   # model specifics
   if MOD=='TP4':
         Model    = 'TP4a0.12'
         topo_dir = M_INP+'/TOPAZ4/'+Model+'/topo/'

         # test one point
         itest = 430
         jtest = 500

   elif MOD=='BS1':
         Model    = 'BS1a0.045'
         topo_dir = M_INP+'/Barents_Hyc2.2.12/topo/'

         # test one point
         itest = 350
         jtest = 350

   elif MOD=='FR1':
         Model    = 'FR1a0.03'
         topo_dir = M_INP+'/FramStrait_Hyc2.2.12/'+Model+'/topo/'

         # test one point
         itest = 200
         jtest = 200
   ###############################################################################################

   ###############################################################################################
   # read binary files
   rfil  = topo_dir+'/regional.grid.a'
   print('\nReading '+rfil)
   plon  = Mr.get_array_from_HYCOM_binary(rfil,1)
   plat  = Mr.get_array_from_HYCOM_binary(rfil,2)
   dims  = plon.shape

   afil     = piv_path+'/pivots/'+Model+'_pivots_'+dname+'.a'
   print('\nReading '+afil)
   ipiv  = Mr.get_array_from_HYCOM_binary(afil,1,dims=dims)
   jpiv  = Mr.get_array_from_HYCOM_binary(afil,2,dims=dims)
   dist  = Mr.get_array_from_HYCOM_binary(afil,3,dims=dims)

   tlon  = plon[itest-1,jtest-1]
   tlat  = plat[itest-1,jtest-1]
   ###############################################################################################

   ###############################################################################################
   # interpolate variable
   print('Interpolating variable')
   vinterp  = np.nan*np.zeros(dims)
   idm,jdm  = dims

   for i in range(idm):
      for j in range(jdm):
         ii = ipiv[i,j]
         jj = jpiv[i,j]
         if ii*jj>0:
            mask  = VBL.values.mask[jj-1,ii-1]
            if not mask:
               vv             = VBL.values.data[jj-1,ii-1]
               vinterp[i,j]   = vv
               # print(i,j,ii-1,jj-1,vv)

            # else:
            #    print('\nWARNING! '+var_name+' is masked var at:')
            #    print('Data grid:')
            #    print(ii,jj)
            #    print(dlon[jj-1,ii-1],dlat[jj-1,ii-1])
            #    print('Model grid')
            #    print(i,j)
            #    print(plon[i,j],plat[i,j])
            #    print('\n')
   ###############################################################################################

   ###############################################################################################
   # plot interpolated variable
   bmap  = Fplt.start_HYCOM_map(MOD)
   Vi    = np.ma.array(vinterp,mask=np.isnan(vinterp))
   Vmin  = Vi.min()
   Vmax  = Vi.max()
   cmap  = cm.jet
   cmap.set_bad(color='w')

   print('Plotting interpolated variable')
   print('Range in '+var_name+':')
   print(Vmin,Vmax)
   fig2  = plt.figure()
   ax2   = fig2.add_subplot(1,1,1)
   PC    = bmap.pcolor(plon,plat,Vi,ax=ax2,latlon=True,\
                        vmin=Vmin,vmax=Vmax)
   fig2.colorbar(PC)
   Fplt.finish_map(bmap,ax=ax2)
   if SHOW:
      fig2.show()
   else:
      figname  = odir+'/'+MOD+'_'+var_name+'_interp.png'
      print('Saving '+figname+'\n')
      fig2.savefig(figname)
      ax2.cla()
      plt.close(fig2)
   ###############################################################################################

   ###############################################################################################
   if 1:
      print('Plotting original variable')
      print('Range in '+var_name+':')
      print(Vmin0,Vmax0)

      if SHOW:
         fig,ax,bmap = nci.plot_var(var_name,HYCOMreg=MOD,time_index=irec,show=SHOW,\
                                    clim=[Vmin,Vmax])
      else:
         fig,ax,bmap = nci.plot_var(var_name,HYCOMreg=MOD,time_index=irec,show=SHOW,\
                                    clim=[Vmin,Vmax])
         #
         figname  = odir+'/'+MOD+'_'+var_name+'.png'
         print('Saving '+figname+'\n')
         fig.savefig(figname)
         #
         ax.cla()
         plt.close(fig)
   ###############################################################################################


# 
# ip    = ipiv(itest,jtest);
# jp    = jpiv(itest,jtest);
# if ip==0
#    disp('itest,jtest are out of the data grid');
#    return;
# end
# 
# figure(3);
# glon  = lon(jp,ip);
# glat  = lat(jp,ip);
# 
# %%
# disp('Single point test:')
# disp('(can also use this to check HYCOM is reading/interpolating data correctly)')
# itest,jtest
# piv_test = [ip jp]
# 
# ll_test  = [tlon tlat;glon glat]
# dtest_km = [dist(itest,jtest) GEN_great_circle_dist(tlon,tlat,glon,glat)]/1e3
# disp(['test <<' var_name '>> at itest/jtest:']);
# vtest    = [vinterp(itest,jtest) vbl(jp,ip)]
# %%
# if test_xtra_vars
#    for j=1:Nv
#       %% get variable;
#       if rank_var==3%%3d (time-dep) variable;
#          [vbl0,sf,ao]   = NCget_v2(ncfil3{j},{var_name2{j},3});
#          vbl            = squeeze(vbl0(irec,:,:));
#          clear vbl0;
#       else
#          [vbl,sf,ao] = NCget_v2(ncfil3{j},{var_name2{j},2});
#       end
# 
#       %%scale etc
#       if ~isempty(sf)
#          vbl   = vbl*sf;
#       end
#       if ~isempty(ao)
#          vbl   = vbl+ao;
#       end
# 
#       disp(['test <<' var_name2{j} '>> at itest/jtest:']);
#       vtest    = vbl(jp,ip)
#    end
# end
# 
# disp('Plotting grids:')
# m_proj('stereographic','lon',tlon,'lat',tlat,'rad',rad2,'rec','on','rot',0);
# 
# [X,Y] = m_ll2xy(plon,plat);
# plot(X,Y,'or');
# disp([Model,' = red']);
# hold on;
# 
# [X,Y] = m_ll2xy(lon,lat);
# plot(X,Y,'.g');
# 
# for j=1:length(Model)
#    ss(j) = ' ';
# end
# ss(1:4)  = 'Data';
# disp([ss ' = green']);
# 
# [X,Y] = m_ll2xy(tlon,tlat);
# plot(X,Y,'^');
# [X,Y] = m_ll2xy(glon,glat);
# plot(X,Y,'v');
# 
# m_grid('fontname','times','fontsize',fs,'linewidth',lw);
# m_coast('patch',[ .2 .2 .2]);
# box off;
# hold off;
