%% test_pivots.m
%% Author: Timothy Williams
%% Date:   20130527, 15:17:50 CEST
%% test the pivots obtained from get_pivots.sh etc
clear
addpath ../../../../matlab_testprogs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%work out filename for test date, and get record number (irec) if necessary:
test_year   = {2012,'2012'};
test_mon    = {5,'05'}; 
wave_dt0    = 6;%h
wave_dt     = wave_dt0/24;
test_day    = 2+3*wave_dt;%%NB starts from 0
%%
test_day2   = {1+floor(test_day)};%%calendar day of month
if test_day2{1}>=10
   test_day2{2}   = num2str(test_day2{1});
else
   test_day2{2}   = ['0' num2str(test_day2{1})];
end
cal_date = [test_year{2} test_mon{2} test_day2{2}];
irec     = 1+round((test_day-floor(test_day))/wave_dt);

var_name    = 'significant_wave_height';
dname       = 'wamnsea';
dpath       = '/work/shared/nersc/msc/WAMNSEA/2012/analysis';%%path of input data file
%dpath       = '/work/shared/nersc/msc/WAVES_INPUT/WAMNSEA/';%%path of input data file
ncfil       = [dpath '/wam_nsea.an.' cal_date '.nc'];%%full name of input file with data lon/lat
ncfil2      = ncfil;%%full name of input file with data

test_xtra_vars = 1;
if test_xtra_vars
   var_name2   = {'peak_wave_period' 'wave_direction'};
   Nv          = length(var_name2);
   for j=1:Nv%%check values of other variables at (itest,jtest)
      %%local names of their input data files
      ncfil3{j}   = ncfil2;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_lon_name  = 'longitude';
d_lat_name  = 'latitude';
rank_lonlat = 2;%%rank of lon/lat fields (1 or 2);
rank_var    = 3;%%rank of the variable (3 if space-time-dep, 2 if space-dep);
d_missing   = 9.969209968386869e+36;%value given to missing values

M_INP = '/work/shared/nersc/msc/ModelInput/';

%%choose model;
mod_no   = 3;
switch mod_no
   case 1
      Model    = 'TP4a0.12' 
      topo_dir = [M_INP 'TOPAZ4/' Model '/topo/'];
      idm      = 800;
      jdm      = 880;

      %%for plotting:
      clon  = -15;
      clat  = 80;
      rad   = 40;
      rad2  = 2;

      %%test one point
      if 0
         [itest,jtest]  = find(ipiv==1);
         if 0%~isempty(itest)
            itest = itest(1);
            jtest = jtest(1);
         else
            itest = round(idm*.4);
            jtest = round(jdm*.5);
         end
      else
         itest = 430;
         jtest = 500;
      end

   case 2
      Model    = 'BS1a0.045'
      topo_dir = [M_INP 'Barents_Hyc2.2.12/topo/'];
      idm      = 510;
      jdm      = 450;

      %%for plotting:
      clon  = 50;
      clat  = 74;
      rad   = 18;
      rad2  = 1;

      %%test one point
      if 1
         [itest,jtest]  = find(ipiv==1);
         if 0%~isempty(itest)
            itest = itest(1);
            jtest = jtest(1);
         else
            itest = round(idm*.3);
            jtest = round(jdm*.5);
         end
      else
         itest = 350;
         jtest = 350;
      end

   case 3
      Model    = 'FR1a0.03'
      topo_dir = [M_INP 'FramStrait_Hyc2.2.12/' Model '/topo/'];
      idm      = 400;
      jdm      = 320;

      %%for plotting:
      clon  = 0.5;
      clat  = 78.5;
      rad   = 8.2;
      rad2  = .5;

      %%test one point
      if 0
         [itest,jtest]  = find(ipiv==1);
         if 0%~isempty(itest)
            itest = itest(1);
            jtest = jtest(1);
         else
            itest = round(idm*.3);
            jtest = round(jdm*.5);
         end
      else
         itest = 240;
         jtest = 80;
      end

end

%%get model grid
gridfile = [topo_dir 'regional.grid.a']
plon     = loada(gridfile,1,idm,jdm);
plat     = loada(gridfile,2,idm,jdm);

%%get pivots
%pivdir   = [dpath '/pivots/'];
pivdir   = './';
pivfile  = [pivdir Model '_pivots_' dname '.a']
ipiv     = loada(pivfile,1,idm,jdm);
jpiv     = loada(pivfile,2,idm,jdm);
dist     = loada(pivfile,3,idm,jdm);

tlon  = plon(itest,jtest);
tlat  = plat(itest,jtest);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1%%plot original variable;
   warning off
   figure(1);
   m_proj('stereographic','lon',clon,'lat',clat,'rad',rad,'rec','on','rot',0);
   fs = 24;%%font size;
   lw = 1.0;%% line width of grid lines;

   %%get lon/lat;
   lon   = NCget(ncfil,{d_lon_name,rank_lonlat});
   lat   = NCget(ncfil,{d_lat_name,rank_lonlat});
   if rank_lonlat==1
      [lon,lat]   = meshgrid(lon,lat);
   end
   nx = size(lon,1);
   ny = size(lon,2);

   %% get variable;
   if rank_var==3%%3d (time-dep) variable;
      [vbl0,sf,ao] = NCget_v2(ncfil2,{var_name,3});
      vbl   = squeeze(vbl0(irec,:,:));%%1st time;
      clear vbl0;
   else
      [vbl,sf,ao] = NCget_v2(ncfil2,{var_name,2});
   end

   j_missing      = find(vbl==d_missing);
   vbl(j_missing) = NaN;
   if ~isempty(sf)
      vbl   = vbl*sf;
   end
   if ~isempty(ao)
      vbl   = vbl+ao;
   end
   if 1
      %lon_rng  = [min(min(lon)) max(max(lon))]
      %lat_rng  = [min(min(lat)) max(max(lat))]
      var_rng  = [min(min(vbl)) max(max(vbl))]
   end

   P     = m_pcolor(lon,lat,vbl);

   if 1
      %%show where the test point is;
      hold on;
      [x0,y0]  = m_ll2xy(tlon,tlat);
      plot(x0,y0,'xb');
      plot(x0,y0,'ob');
      hold off;
   end

   set(P,'linestyle','none');
   shading flat;

   m_grid('fontname','times','fontsize',fs,'linewidth',lw);
   m_coast('patch',[ .2 .2 .2]);
   box off;

   colorbar
   colormap cisconc
end

%%plot interpolated variable;
figure(2);

vinterp  = zeros(idm,jdm);
for i=1:idm
   for j=1:jdm
      ii = ipiv(i,j);
      jj = jpiv(i,j);
      if ii*jj>0
         vv             = vbl(jj,ii);
         vinterp(i,j)   = vv;
   %     if vv>0
   %        disp(vv);
   %     end
      else
         vinterp(i,j)   = NaN;
      end
   end
end

m_proj('stereographic','lon',clon,'lat',clat,'rad',rad,'rec','on','rot',0);
P     = m_pcolor(plon,plat,vinterp);

if 1
   %%show where the test point is;
   hold on;
   [x0,y0]  = m_ll2xy(tlon,tlat);
   plot(x0,y0,'xb');
   plot(x0,y0,'ob');
   hold off;
end

set(P,'linestyle','none');
shading flat;

m_grid('fontname','times','fontsize',fs,'linewidth',lw);
m_coast('patch',[ .2 .2 .2]);
box off;

colorbar
colormap cisconc

vmax  = max(max(vinterp));
vmin  = min(min(vinterp));
caxis([vmin vmax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
caxis([vmin vmax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ip    = ipiv(itest,jtest);
jp    = jpiv(itest,jtest);
if ip==0
   disp('itest,jtest are out of the data grid');
   return;
end

figure(3);
glon  = lon(jp,ip);
glat  = lat(jp,ip);

%%
disp('Single point test:')
disp('(can also use this to check HYCOM is reading/interpolating data correctly)')
itest,jtest
piv_test = [ip jp]

ll_test  = [tlon tlat;glon glat]
dtest_km = [dist(itest,jtest) GEN_great_circle_dist(tlon,tlat,glon,glat)]/1e3
disp(['test <<' var_name '>> at itest/jtest:']);
vtest    = [vinterp(itest,jtest) vbl(jp,ip)]
%%
if test_xtra_vars
   for j=1:Nv
      %% get variable;
      if rank_var==3%%3d (time-dep) variable;
         [vbl0,sf,ao]   = NCget_v2(ncfil3{j},{var_name2{j},3});
         vbl            = squeeze(vbl0(irec,:,:));
         clear vbl0;
      else
         [vbl,sf,ao] = NCget_v2(ncfil3{j},{var_name2{j},2});
      end

      %%scale etc
      if ~isempty(sf)
         vbl   = vbl*sf;
      end
      if ~isempty(ao)
         vbl   = vbl+ao;
      end

      disp(['test <<' var_name2{j} '>> at itest/jtest:']);
      vtest    = vbl(jp,ip)
   end
end

disp('Plotting grids:')
m_proj('stereographic','lon',tlon,'lat',tlat,'rad',rad2,'rec','on','rot',0);

[X,Y] = m_ll2xy(plon,plat);
plot(X,Y,'or');
disp([Model,' = red']);
hold on;

[X,Y] = m_ll2xy(lon,lat);
plot(X,Y,'.g');

for j=1:length(Model)
   ss(j) = ' ';
end
ss(1:4)  = 'Data';
disp([ss ' = green']);

[X,Y] = m_ll2xy(tlon,tlat);
plot(X,Y,'^');
[X,Y] = m_ll2xy(glon,glat);
plot(X,Y,'v');

m_grid('fontname','times','fontsize',fs,'linewidth',lw);
m_coast('patch',[ .2 .2 .2]);
box off;
hold off;
