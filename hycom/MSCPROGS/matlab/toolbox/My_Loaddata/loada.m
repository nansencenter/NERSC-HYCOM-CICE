function [fld,lon,lat,depths]=loada(afile,rec,idm,jdm);
%function [fld,lon,lat,depths]=loada(afile,rec)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Routine loada:
%%%%%%%%%%%%%%%%%
%This routine reads data from 4-byte big-endian daily average files (.a - files)
%The filename and the record number is supplied .
%
%The routine also tries to read the depths from the model depths-file
%and lon/lat from the file newpos.uf. If these files are absent, the
%fields lon,lat and/or depths will be empty.
%
%NB: This is a dumb routine.. If you need to read restart/weekly/daily files,
%use loadweekly/loaddaily/loadrestart
%
%Examples: 
% 
%To read 3rd record from a file
%[fld,lon , lat, depths]=loada('N32DAILY_1958_000_1958_000.b',3,idm,jdm);
%
%Knut Liseter, 30.10.2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lon=[];
lat=[];
depths=[];
fld=[];
onem=9806;
listmode=0;

if (nargin~=4  )
   disp(['loada needs 4 input  arguments']);
   return
end

%disp(['Opening data file ' afile ]);
% A quirk of the mod_za module, data is dumped containing a whole multiple of 4096 values.
% This is used to skip to the correct record, but we will only read idm*jdm values...
n2drec=floor((idm*jdm+4095)/4096)*4096;
bytes_per_float=4;
fid=fopen(afile,'r','ieee-be'); % Big-endian

%layerind_file

% Skip to indices in layerind_file
fld=zeros(idm,jdm);
stat=fseek(fid,n2drec*bytes_per_float*(rec-1),'bof'); % Skip to corr record
if (stat==0)
   [fld,count]=fread(fid,[idm jdm],'single'); 
else
   disp(['Couldnt skip to the indicated record in ' afile ])
   disp( 'You probably tried to read past the End-Of-File')
   fld=[];
end
fclose(fid);



if (nargout>2)
   % Try to retrieve lon/lat from newpos.uf
   fid=fopen('newpos.uf','r','ieee-be');
   if (fid~=-1)
      stat=fseek(fid,4,'bof'); % Skip fortran 4-byte header 
      lat=fread(fid,[idm jdm],'double');
      lon=fread(fid,[idm jdm],'double');
      %contourf(lon)
      fclose(fid);
   else
      try
         [lon]=loada('regional.grid.a',1,idm,jdm);
         [lat]=loada('regional.grid.a',2,idm,jdm);
      catch
         disp(['newpos.uf or regional.grid.a not found -- lon lat will be empty']);
         lon=[];
         lat=[];
      end
   end
end

if (nargout>1)
   fdepths=['depths' num2str(idm,'%3.3i') 'x' num2str(jdm,'%3.3i') '.uf'];
   fid=fopen(fdepths, 'r','ieee-be');
   if (fid~=-1)
	  disp(['depths file is ' fdepths]);
      stat=fseek(fid,4,'bof'); % Skip fortran 4-byte header 
      depths=fread(fid,[idm jdm],'double');
      fclose(fid);
   else
      disp([fdepths ' not found -- depths will be empty']);
      depths=[];
   end
end

