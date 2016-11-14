function [dist,depthc,fld,lon,lat]=load_p2s(filename,varname,record)
%[dist,depthc,fld[,lon,lat]]=load_p2s(filename,varname,record)
%
%Routine reads variables lon, lat, distance along sections, depth interfaces and a 
%variable corresponding to the sone specified by "varname" Input files are thse created by 
%utility "pak2sec".
%
%Input args:
%  -filename: The name of the pak2sec-generated file to read (typically section001.nc
%             section002.nc ... )
%  -varname : Name of variable to read (ex 'SAL' )
%  -record  : The record number to read (if you use pak2sec on several files)
%
%Output vars:
%  -dist    : distance along section
%  -depthc  : depth interfaces
%  -fld     : variable corresponding to varname
%  -lon     : longitude (optional)
%  -lat     : latitude  (optional)
%
%All output variables are formatted for contour plotting (lon, lat, dist, depthc a
%and fld will all have the same dimension... ). Use pcolor (pcolor(dist,depthc,fld)) 
%for most realistic plotting.
%
%Example:
%Read variable salinity from section001.nc, use record number 1
%[dist,depthc,fld,lon,lat]=load_p2s('section001.nc','SAL',1);
%pcolor(dist,depthc,fld) % - To plot ...
%
%NB: If you give a variable name not present in the section-file,
%    you will be presented with a list of possible variable names.

dist=[];
depthc=[];
fld=[];
lon=[];
lat=[];

% Open file
nc=netcdf(filename);
if (isempty(nc))
   disp (['No such file ' filename ]);
   return
end

% check for variable corr to "varname"
v=nc{varname,1};
if (isempty(v))
   disp (['No such variable ' varname ' in file ' filename ]);
   disp('Available variables are:')
   tmp=var(nc);
   for i=1:size(tmp,2)
      disp(char(ncnames(tmp(i))))
   end
   return
end

% Check that record dimension < record
nrec=size(v,1);
if (record>nrec | record<1 )
   disp(['Record to small or large, nrec = ' num2str(nrec)])
   return
end

% Get vars 
lon=nc{'longitude'}(:);
lat=nc{'latitude'}(:);
dist=nc{'distance'}(:);
depthc=nc{'depthc'}(record,:,:);
fld=v(record,:,:);
fillv=fillval(v);
I=find(fld==fillv);
fld(I)=nan;

% Find Last mass-containing layer 
depthdim=size(fld,1) ;
distdim=size(fld,2) ;
for i=1:distdim
   I=find(~isnan(fld(1:depthdim,i)));
   if (~isempty(I))
      I=max(max(I));
      %fld(I:depthdim,i)=fld(I,i);
   end
end

%for i=2:distdim-1
%for j=2:depthdim-1
%   if (depthc(j,i+1)-depthc(j+1,i+1)<12) 
%      fld(i,j)=nan;
%   end 
%end
%end


%size(lon)
%size(lat)
%size(depthc)
%size(fld)

% Replicate lon lat dist along depth dimension
%depthdim=prod(size(fld))/prod(size(lon)) ; % Hee hee
lon=repmat(lon,1,depthdim)';
lat=repmat(lat,1,depthdim)';
dist=repmat(dist,1,depthdim)';
dist=dist/1000 ; % To km



disp([ 'Ok - use pcolor for most realistic display']);
