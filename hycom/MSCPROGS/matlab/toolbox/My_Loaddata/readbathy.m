function [ depths ] = readbathy(idm,jdm,old);
%function writebathy(depths);
%
% Routine to read a bathymetry to  matlab 


% This opens the unformatted file type "depthsXXXxXXX.uf" - sequential access
fname=['depths' num2str(idm,'%3.3d') 'x'  num2str(jdm,'%3.3d') '.uf' ];
fid=fopen(fname,'r','ieee-be');
if (fid~=-1)
   nent=fread(fid,1,'integer*4');  % Header - Total read to be written
   depthsuf=fread(fid,[ idm jdm],'double');
   nent2=fread(fid,3,'integer*4'); % "Trailer" - for what its worth
   fclose(fid);
else
   depthsuf=[];
end

% This opens .b and .a (direct access) - the latter is easier to dump to
fida=fopen('regional.depth.a','r','ieee-be');
if (fida~=-1) 
   n2drec=floor((idm*jdm+4095)/4096)*4096;
   npad=n2drec-idm*jdm;
   depthsa=fread(fida,[idm jdm],'float');
   fclose(fida);
else
   disp('Can not open regional.depth file ')
   depthsa=[];
end 

if (~isempty(depthsuf) & ~isempty(depthsa))
   if (abs(max(max(depthsuf-depthsa)))>1e-4) 
      disp('Warning depths file mismatch between regional and depths...uf')
   end
end

%default - return depthsa
depths=depthsa;

% Unless uf is requested
if (nargin==3 & old==1 )
   depths=depthsuf ; % could be empty
end
