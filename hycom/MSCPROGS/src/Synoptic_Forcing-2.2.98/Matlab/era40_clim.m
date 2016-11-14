function era40_clim(var_file)

% Variable to process.
if (strcmp(var_file,'2T'))
   var_id='T2M_sfc'    % 
   var_id_controller  = var_id;
   var_tp='ans'
elseif (strcmp(var_file,'CI'))
   var_id='CI_sfc' %
   var_id_controller  = var_id;
   var_tp='ans'
elseif (strcmp(var_file,'relhum'))
   var_id='relhum'
elseif (strcmp(var_file,'TCC'))
   var_id='TCC_sfc' %
   var_id_controller  = var_id;
   var_tp='ans'
elseif (strcmp(var_file,'TP'))
   var_id='TP' %
   var_id_controller  = var_id;
   var_tp='fcs'
elseif (strcmp(var_file,'MSL'))
   var_id='MSL_sfc' %
   var_id_controller  = var_id;
   var_tp='ans'
elseif (strcmp(var_file,'10U'))
   var_id='U10M_sfc' %
   var_id_controller  = var_id;
   var_tp='ans'
elseif (strcmp(var_file,'10V'))
   var_id='V10M_sfc' %
   var_id_controller  = var_id;
   var_tp='ans'
elseif (strcmp(var_file,'taux'))
   var_id  ='taux';
elseif (strcmp(var_file,'tauy'))
   var_id  ='tauy';
elseif (strcmp(var_file,'wndspd'))
   var_id  ='wndspd';
else
   disp(['Unknown var_file ' var_file]);
   return
end

files4=[];
files3=[];
files2=[];
files=[];
if (strcmp(var_id,'taux') | strcmp(var_id,'tauy') | strcmp(var_id,'wndspd')  )
   var_tp1='ans'; var_file1='10U'; 
   var_tp2='ans'; var_file2='10V';
   files =dir([ var_tp1 '.6h.*.' var_file1 '.nc']);
   files2=dir([ var_tp2 '.6h.*.' var_file2 '.nc']);
   var_id1  ='U10M_sfc';
   var_id2  ='V10M_sfc';
   var_id_controller  =var_id1;
elseif (strcmp(var_id,'relhum'))
   var_tp1='ans'; var_file1='2T';  var_id1='T2M_sfc'
   var_tp2='ans'; var_file2='2D';  var_id2='D2M_sfc'
   var_tp3='ans'; var_file3='MSL'; var_id3='MSL_sfc'
   var_tp4='ans'; var_file4='CI'; var_id4='CI_sfc'
   files =dir([ var_tp1 '.6h.*.' var_file1 '.nc']);
   files2=dir([ var_tp2 '.6h.*.' var_file2 '.nc']);
   files3=dir([ var_tp3 '.6h.*.' var_file3 '.nc']);
   files4=dir([ var_tp4 '.6h.*.' var_file4 '.nc']);
   var_id_controller  =var_id1;
else
   files=dir([ var_tp '.6h.*.' var_file '.nc'])
end

% Go through files
meanfld=[];
maxyear=1000;
minyear=3000;
cd=0.0012;
airdns=1.2;
mcnt=zeros(12);
for i=1:size(files,1);
%for i=1:1 % Test
%for i=35:36 % Test

   disp(files(i).name);
   if (~isempty(files2))
      disp(files2(i).name);
   end
   if (~isempty(files3))
      disp(files3(i).name);
   end
   if (~isempty(files4))
      disp(files4(i).name);
   end

   % Open netcdf file
   if (strcmp(var_id,'taux') | strcmp(var_id,'tauy') | strcmp(var_id,'wndspd')  )
      nc=netcdf(files(i).name);
      nc2=netcdf(files2(i).name);
   elseif (strcmp(var_id,'relhum') )
      nc=netcdf(files(i).name);
      nc2=netcdf(files2(i).name);
      nc3=netcdf(files3(i).name);
      nc4=netcdf(files4(i).name);
   else
      nc=netcdf(files(i).name);
   end

   % Go through fields
   for j=1:size(nc{var_id_controller},1);

      % Read time info
      rtime=nc{'valtime'}(j);


      % Time is in hours since 1992-01-01
      %[Y,M,D,H] = datevec(rtime/24.+1.,1992);
      %Y=Y+1992;
      datepiv='01/01/1992';
      datenum0= datenum(datepiv);
      [Y,M,D,H] = datevec(rtime/24.+datenum0);
      




      if (mod(j,200)==0)
      disp([ num2str(rtime,'%9.2f') ' ' num2str(Y,'%4.4d') ' ' num2str(M,'%2.2d') ...
             ' ' num2str(D,'%2.2d') ' ' num2str(H,'%2.2d')])
       end

      if (strcmp(var_id,'taux') | strcmp(var_id,'tauy') | strcmp(var_id,'wndspd')  )
         U10 =nc{var_id1,1}(j,:,:);
         V10 =nc2{var_id2,1}(j,:,:);
         fld=U10;
      elseif (strcmp(var_id,'relhum') )
         T2  =nc{var_id1,1}(j,:,:);
         D2  =nc2{var_id2,1}(j,:,:);
         MSL =nc3{var_id3,1}(j,:,:);
         CI  =nc4{var_id4,1}(j,:,:);
         fld=T2;
      else
         fld=nc{var_id,1}(j,:,:);
      end

      if (isempty(meanfld)) % First pass
         meanfld=zeros(12,size(fld,1),size(fld,2));
         lon=nc{'lon'}(:);
         lat=nc{'lat'}(:);
      end 


      if (strcmp(var_id,'taux') | strcmp(var_id,'tauy') | strcmp(var_id,'wndspd')  )
         U10=reshape(U10,1,size(U10,1),size(U10,2));
         V10=reshape(V10,1,size(V10,1),size(V10,2));
         wspd=sqrt(U10.*U10+V10.*V10);
         wndfac=(1.+sign(wspd-11))*.5;
         cd_new=(0.49+0.065*wspd).*1e-3.*wndfac + cd.*(1-wndfac);

         if (strcmp(var_id,'taux'))
            meanfld(M,:,:)=meanfld(M,:,:)+wspd.*U10.*cd_new.*airdns;
         elseif (strcmp(var_id,'tauy'))
            meanfld(M,:,:)=meanfld(M,:,:)+wspd.*V10.*cd_new.*airdns;
         elseif (strcmp(var_id,'wndspd'))
            meanfld(M,:,:)=meanfld(M,:,:)+wspd;
         end
      elseif (strcmp(var_id,'relhum') )
         T2=reshape(T2,1,size(T2,1),size(T2,2));
         D2=reshape(D2,1,size(D2,1),size(D2,2));
         MSL=reshape(MSL,1,size(MSL,1),size(MSL,2));
         CI=reshape(CI,1,size(CI,1),size(CI,2));
         svpair=satvap(T2,CI);
         svpdew=satvap(D2,CI);
         relhum=relhumid(svpair,svpdew,MSL).*0.01;
         meanfld(M,:,:)=meanfld(M,:,:)+relhum;
      else
         meanfld(M,:,:)=meanfld(M,:,:)+reshape(fld,1,size(fld,1),size(fld,2));
      end
      mcnt(M)=mcnt(M)+1;

      maxyear=max(Y,maxyear);
      minyear=min(Y,minyear);


   end
   if (strcmp(var_id,'taux') | strcmp(var_id,'tauy') )
      close(nc);
      close(nc2);
   elseif (strcmp(var_id,'relhum') )
      close(nc);
      close(nc2);
      close(nc3);
      close(nc4);
   else
      close(nc)
   end
end

for M=1:12
   meanfld(M,:,:)=meanfld(M,:,:)/mcnt(M);
end
   

% Create a climatology file - use first file as basis
nc=netcdf(files(1).name);
nc2=netcdf(['era40_climatology_' var_id '.nc'],'clobber');

copy(nc('lon'),nc2);
copy(nc('lat'),nc2);
copy(nc('nav'),nc2);
nc2('month') = 0;

copy(nc{'lon'},nc2);
copy(nc{'lat'},nc2);
copy(nc{'Ni'},nc2);
copy(nc{'Nj'},nc2);
copy(nc{'La1'},nc2);
copy(nc{'Lo1'},nc2);
copy(nc{'La2'},nc2);
copy(nc{'Lo2'},nc2);
copy(nc{'Di'},nc2);
copy(nc{'Dj'},nc2);

nc2{'lon'}(:)=nc{'lon'}(:);
nc2{'lat'}(:)=nc{'lat'}(:);
nc2{'lat'}(:)=nc{'lat'}(:);
nc2{'Ni'}(:)=nc{'Ni'}(:);
nc2{'Nj'}(:)=nc{'Nj'}(:);
nc2{'Di'}(:)=nc{'Di'}(:);
nc2{'Dj'}(:)=nc{'Dj'}(:);
nc2{'Lo1'}(:)=nc{'Lo1'}(:);
nc2{'La1'}(:)=nc{'La1'}(:);
nc2{'Lo2'}(:)=nc{'Lo2'}(:);
nc2{'La2'}(:)=nc{'La2'}(:);

nc2{var_id} = {'month','lat','lon'};
nc2{'month'} = {'month'};
if (strcmp(var_id,'taux') | strcmp(var_id,'tauy') | strcmp(var_id,'wndspd')  )
   if (strcmp(var_id,'taux'))
      nc2{var_id}.units='N m**-2';
      nc2{var_id}.long_name='Eastward stress component';
   elseif (strcmp(var_id,'tauy'))
      nc2{var_id}.units='N m**-2';
      nc2{var_id}.long_name='Northward stress component';
   elseif (strcmp(var_id,'wndspd'))
      nc2{var_id}.units='m s**-1';
      nc2{var_id}.long_name='Wind Speed';
   end
elseif (strcmp(var_id,'relhum,'))
   nc2{var_id}.units='';
   nc2{var_id}.long_name='Relative humidity';
else
   copy(nc{var_id}.units,nc2{var_id});
   copy(nc{var_id}.long_name,nc2{var_id});
   copy(nc{var_id}.navigation,nc2{var_id});
   copy(nc{var_id}.FillValue_,nc2{var_id});
end
copy(nc.Conventions,nc2);
nc2.Title='ERA40 Monthly Climatology';
for M=1:12
   nc2{var_id}(M,:,:) = meanfld(M,:,:); 
   nc2{'month'}(M)=M;
end
close(nc)
close(nc2)



function relhumid=relhumid(sva,svd,msl);
%real function relhumid(sva,svd,msl)
%! This routine calculates the relative humidity by the 
%! dew point temperature and the mean sea level pressure.
%! Modified: Anita Jacob, June '97
%
%! Input: sva: saturatn vapour press at air temp [K]
%! svd: saturatn vapour press at dew pt temp [K]
%! msl: pressure at mean sea level [Pa]
%! Output: relhumid: Relative Humidity
%
%! We use the Tetens formula:
%! es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual
%!              es(Tdew)        p - es(Tair)
%! RH = 100 *  -----------   *  ------------
%!             p - es(tdew)       es(Tair)

      aaa=msl - svd;
      aaa = svd./aaa;
      bbb = (msl - sva)./sva;
      relhumid = 100 .* aaa .* bbb;







function satvap=satvap(t,ci)
%! This function calculates the saturation vapour pressure
%! [Pa] from the temperature [deg K].
%! Modified: Anita Jacob, June '97
%
%! Input: t: temperature [deg K]
%! Output: satvap: saturation vapour pressure at temp. t
%
%! es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual


      c1=610.78;
      t00=273.16;

      %if (t < t00) then
      %   c3 = 21.875
      %   c4 = 7.66
      %else
      %   c3 = 17.269
      %   c4 = 35.86
      %endif
      %iceflag=(1.-sign(t-t00))/2.
      iceflag=(1.-sign(.03-ci))/2.;

      %figure(1);
      %[C,H]=contourf(reshape(ci,size(ci,2),size(ci,3)),0:.01:1); set(H,'LineStyle','none');
      %hold on;
      %[C,H]=contour(reshape(iceflag,size(iceflag,2),size(iceflag,3)),[1 1]);
      %set(H,'EdgeColor','k') ; set(H,'LineWidth',3);
      %pause(1);

      c3=iceflag.*21.875 + (1.-iceflag).*17.269;
      c4=iceflag.* 7.66  + (1.-iceflag).*35.86 ;

      aa = c3 .* (t - t00);
      bb = t - c4;
      
      cc=aa./bb;
      %if (cc < -20.0) then
      %   satvap=0.0
      %else
      %   satvap = c1 * exp(aa/bb)
      %endif
      satvap = c1 .* exp(cc);




