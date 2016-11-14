function [slon,slat,sdist,sdepth,sname,varargout]=loadsection(file,type,secnum,varargin);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [slon,slat,sdist,sdepth,sname,vout1,vout2...]=loadsection(filename,type,secnum,vname1,vname2...);
%
%
% Wrapper for various file-reading tools. Reads daily,restart and weekly average files, then 
% uses data files from the "section_intersect" routine to get points along section.
%
% Output arguments: lon, lat, depths, name and a variable number of output fields (eg vout1,vout2 etc etc)
% along the section.
% Input  arguments: 
%                   filename -- a weekly, average or restart file
%                   type     -- specify either 'daily', 'weekly' or 'restart' - depending on what file you have
%                   secnum   -- a number specifying which section to use  - see output from section_intersect
%                   vname1,vname2 .... -- variable names to extract (eg 'saln', 'temp')
%
% Examples:
% --To extract salinity and temperature from a weekly average file for section 1:
%   [slon,slat,sdist,sdepths,sname,ssaln,stemp]=loadsection('TP4AVE_1996_07_1.a','weekly',1,'saln','temp');
%
% --To extract salinity from a daily average file in section 1:
%   [slon,slat,sdist,sdepths,sname,ssaln]=loadsection('C1ADAILY_2006_119_2006_120.a','daily',1,'saln');
%
% --To extract salinity u and v from a daily average file in section 1:
%   [slon,slat,sdist,sdepths,ssaln,sname,utot,vtot]=loadsection('C1ADAILY_2006_119_2006_120.a','daily',1,'saln','utot','vtot');
%
% NB: Number of output arguments vout1, vout2, ... must match number of input variable names vname1,vname2...
%     Also, vout1 will contain fields of variable specified in vname1.
%     To see what variable names are in a given file, look at the ".b" extension, or look at routines
%     loaddaily, loadrestart and loadweekly with the 'list' option.
%
%     Finally, this routine works on 3d fields, it will probably complain if you try to read 2d fields (ex barotropic
%     velocities).
%
% NB2: Remember to run section_intersect first (or pak2sec which uses this rutine as well). This routine reads
%      from data files "section001.dat", "section002.dat" (depending on the number of the sectiom, secnum), which
%      are created by section_intersect
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

nvarargout=nargout-5;
nvarargin=nargin-3;
if (nvarargout==0) 
   disp(['Number of variable output args must be > 0)'])
   return
elseif (nvarargout~=nvarargin)
   disp(['Number of variable output args must must match variable input args'])
   return
end


%Read section intersect file
secname=['section' num2str(secnum,'%3.3i') '.dat'];
fid=fopen(secname);
if (fid<0)
   disp(['Can not open section file ' secname ' You should run section_intersect first'])
   return
end

% First line has number of values
tline = fgetl(fid);
numentry=str2num(tline(2:6));

% Remaining lines contain section details
for i=1:numentry
   tline = fgetl(fid);
   [ipiv(i) jpiv(i) slon1(i) slat1(i) sdist1(i) tmp1 tmp2 sname] = strread(tline,'%d%d%f%f%f%f%f%s',1);
   %A = sscanf(tline,'%f',5);
   %ipiv(i) = A(1);
   %jpiv(i) = A(2);
   %slon(i) = A(3);
   %slat(i) = A(4);
   %cdist(i)= A(5);
end
%sname

%index into 2D matrix
PIV=sub2ind([idm jdm],ipiv,jpiv);
sdepth=[];
sdist=[];
slon=[];
slat=[];
for ivar=1:nvarargout
   varargout{ivar}=[];
end
for k=1:kdm


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
      sdepth=[ sdepth ; zeros(size(dp(PIV))) ];
      p=dp;
      sdist=[ sdist ; sdist1];
      slon=[ slon ; slon1];
      slat=[ slat ; slat1];
   else
      p=p+dp;
   end


   sdepth=[ sdepth ; p(PIV)];
   sdist=[ sdist ; sdist1];
   slon=[ slon ; slon1];
   slat=[ slat ; slat1];




   for ivar=1:nvarargout
      % Always read pressures
      if(strcmp(type,'restart'))
         fldin=loadrestart(file,varargin{ivar},k);
      elseif(strcmp(type,'weekly'))
         fldin=loadweekly(file,varargin{ivar},k);
      elseif(strcmp(type,'daily'))
         fldin=loaddaily(file,varargin{ivar},k);
      end
      if (k==1)
         varargout{ivar}=[ varargout{ivar} ; fldin(PIV) ];
      end
      varargout{ivar}=[ varargout{ivar} ; fldin(PIV) ];
   end 
end
