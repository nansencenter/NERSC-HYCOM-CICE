function [lon,lat,depths,varargout]=loadlevel2(file,type,level,varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [lon,lat,depths,fld1,fld2,...]=loadlevel2(file,type,level,varargin,vname1,vname2...);
%
%
% Wrapper for various file-reading tools. Reads daily,restart and weekly average files, then 
% calculates the value at a given depth level. For now there is now vertical interpolation
% so the results might appear "choppy"..
%
% Output arguments: lon, lat, depths, and a variable number of output fields (eg saln,temp etc etc)
% Input  arguments: 
%                   filename -- a weekly, average or restart file
%                   type     -- specify either 'nersc_daily', 'nersc_weekly','archv'
%                               or 'restart' - depending on what file you have
%                   level    -- at what depth level to extract data
%                   vname1,vname2 .... -- variable names to extract (eg 'saln', 'temp')
%
% Examples:
% --To extract salinity and temperature from a weekly average file at 100 meter depth:
%   [lon,lat,depths,saln,temp]=loadlevel2('TP4AVE_1996_07_1.a','weekly',100,'saln','temp');
%
% --To extract salinity from a daily average file at 100 meter depth:
%   [lon,lat,depths,saln]=loadlevel('C1ADAILY_2006_119_2006_120.a','daily',100,'saln');
%
% --To extract salinity u and v from a daily average file at 100 meter depth:
%   [lon,lat,depths,saln,utot,vtot]=loadlevel('C1ADAILY_2006_119_2006_120.a','daily',100,'saln','utot','vtot');
%
% NB: Number of output arguments fld1, fld2, ... must match number of input variable names vname1,vname2...
%     Also, fld1 will contain fields of variable specified in vname1.
%     To see what variable names are in a given file, look at the ".b" extension, or look at routines
%     loaddaily, loadrestart and loadweekly with the 'list' option.
%
%     Finally, this routine works on 3d fields, it will probably complain if you try to read 2d fields (ex barotropic
%     velocities)
%
% PS: This routine differs from "loadlevel" in that it reads several variables at the same time. It is therefore
%     more efficient than "loadlevel" in that pressure interfaces (or layer thickness) will only be read once. Also,
%     loadlevel2 has logic which will skip reading layers when we dont need them.
%
% Knut Liseter, 16.05.2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

onem=9806;

lon=[];
lat=[];
depths=[];

nvarargout=nargout-3;
nvarargin=nargin-3;
if (nvarargout==0) 
   disp(['Number of variable output args must be > 0)'])
   return
elseif (nvarargout~=nvarargin)
   disp(['Number of variable output args must must match variable input args'])
   return
end

obj=abfile('regional.grid.a','regional_grid');
lon=obj.getfield('plon',[],[]);
lat=obj.getfield('plat',[],[]);
obj=abfile('regional.depth.a','raw');
depths=obj.getfield([],1,[]);

% File to read
obj=abfile(file,type);
kdm=max(max(obj.getlevels()));

cont=1;
k=1;
while (k<=kdm & cont)


   % Always read pressures
   if(strcmp(type,'restart'))
      dp=obj.getfield('dp',k,1);
      dp=dp/onem;
   elseif(strcmp(type,'archv'))
      dp=obj.getfield('thknss',k,1);
      dp=dp/onem;
   elseif(strcmp(type,'nersc_weekly'))
      dp=obj.getfield('pres',k,1);
   elseif(strcmp(type,'nersc_daily'))
      dp=obj.getfield('pres',k,1);
      if (k>1) 
         dp=obj.getfield('pres',k-1,1);
      else
         intf=zeros(size(dp));
      end
      dp=dp-intf;
      dp=dp/onem;
   end

   if (k==1) 
      I=find(dp>.1);
      disp([ 'Thickness of 1st layer: ' num2str(min(min(dp(I)))) ]);
      disp([ '(If this is very high p is in pressure coords = modify level!)']);
   end


   % p is lower interface of layer
   if (k==1)
      p=dp;
      for ivar=1:nvarargout
         varargout{ivar}=zeros(size(p));
      end 
      nofill=ones(size(p)); % Where not to fill with values ...
   else
      p=p+dp;
   end

   I=find(nofill==1); % Points which havent been filled yet
   J=find(p>level);   % Points where lower interface is below level
   K=intersect(I,J);  % Intersection is where to fill with values
   nofill(K)=0;       % Makes sure we wont overwrite already valid values later on

   %Fill in values if size(K)>0
   if (size(K)>0) 
      disp(['Reading layer ' num2str(k) ' of ' num2str(kdm)])
      for ivar=1:nvarargout
         fldin=obj.getfield(varargin{ivar},k,1);
         varargout{ivar}(K)=fldin(K);
      end 
   else
      disp(['Skipping layer ' num2str(k)])
      fldin=zeros(size(p));;
   end 


   % Check if we need to read next level
   L=find(abs(p-depths)<1.);
   M=union(J,L);
   if (prod(size(M)) == prod(size(fldin)))
      cont=0;
      disp('Skipping remaining layers')
   end

   k=k+1;
end


K=find(nofill==1);
for ivar=1:nvarargout
   varargout{ivar}(K)=nan;
end

% Blank out when level>topo
K=find(level>depths);
for ivar=1:nvarargout
   varargout{ivar}(K)=nan;
end
