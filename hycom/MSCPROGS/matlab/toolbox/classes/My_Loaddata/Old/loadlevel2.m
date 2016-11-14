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
%                   type     -- specify either 'daily', 'weekly' or 'restart' - depending on what file you have
%                   level    -- at what depth level to extract data
%                   vname1,vname2 .... -- variable names to extract (eg 'saln', 'temp')
%
% Examples:
% --To extract salinity and temperature from a weekly average file at 100 meter depth:
%   [lon,lat,depths,saln,temp]=loadlevel('TP4AVE_1996_07_1.a','weekly',100,'saln','temp');
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

% Open blkdat.input to retrieve dimensions
fid=fopen([ 'blkdat.input']);
if (fid~=-1)
   for k=1:6
      A=fgetl(fid);
   end
   A=fgetl(fid);
   cidm=A(1:7);
   idm=str2num(cidm);

   A=fgetl(fid);
   cjdm=A(1:7);
   jdm=str2num(cjdm);

   A=fgetl(fid);
   A=fgetl(fid);
   A=fgetl(fid);
   ckdm=A(1:7);
   kdm=str2num(ckdm);

else
   disp(['Error - could not read blkdat.input  - I quit'])
   return
end
fclose(fid);

nvarargout=nargout-3;
nvarargin=nargin-3;
if (nvarargout==0) 
   disp(['Number of variable output args must be > 0)'])
   return
elseif (nvarargout~=nvarargin)
   disp(['Number of variable output args must must match variable input args'])
   return
end

% To get lon, lat, depths
if(strcmp(type,'restart'))
   [fldin, lon , lat, depths]=loadrestart(file,varargin(1),1);
elseif(strcmp(type,'weekly'))
   [fldin, lon , lat, depths]=loadweekly(file,varargin(1),1);
elseif(strcmp(type,'daily'))
   [fldin, lon , lat, depths]=loaddaily(file,varargin(1),1);
else
   disp(['Uknown filetype ' type ]);
   return
end

cont=1;
k=1;
while (k<=kdm & cont)
%for k=1:kdm


   % Always read pressures
   if(strcmp(type,'restart'))
      dp   =loadrestart(file,'dp',k);
      dp=dp/onem;
   elseif(strcmp(type,'weekly'))
      dp   =loadweekly(file,'pres',k);
   elseif(strcmp(type,'daily'))
      dp   =loaddaily(file,'pres',k);
      if (k>1) 
         [intf , lon , lat, depths]=loaddaily(file,'pres',k-1);
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
         % Always read pressures
         if(strcmp(type,'restart'))
            fldin=loadrestart(file,varargin(ivar),k);
         elseif(strcmp(type,'weekly'))
            fldin=loadweekly(file,varargin(ivar),k);
         elseif(strcmp(type,'daily'))
            fldin=loaddaily(file,varargin(ivar),k);
         end
         varargout{ivar}(K)=fldin(K);
      end 
   else
      disp(['Skipping layer ' num2str(k)])
   end 


   % Check if we need to read next level
   L=find(abs(p-depths)<1.);
   M=union(J,L);
   %prod(size(M))
   %prod(size(fldin))
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
