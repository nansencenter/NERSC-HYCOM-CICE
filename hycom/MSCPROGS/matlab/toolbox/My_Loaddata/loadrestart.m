function [fld,lon,lat,depths]=loadrestart(rstfile,varname,layer1,layer2);
%function [fld,lon,lat,depths]=loadrestart(rstfile,varname,layer1,layer2);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Routine loadrestart:
%%%%%%%%%%%%%%%%%
%This routine reads data from 4-byte big-endian restart files (.ab - files)
%The filename and the name of the variable is supplied
%as arguments. Also, the user must supply a layer argument to the routine. 
%For 3D vars, the layer is the model layer. For 2D vars, the layer must
%be equal to 0. You can also specify a range of layers by specifying 
%the first and last layer to be read.
%
%The routine also tries to read the depths from the model depths-file
%and lon/lat from the file newpos.uf. If these files are absent, the
%fields lon,lat and/or depths will be empty.
%
%Examples: 
% 
%To read Salinity for layer 1 from the file  N32DAILY_1958_000_1958_000.[ab], three approaches:
%[fld,lon , lat, depths]=loadrestart('N32DAILY_1958_000_1958_000.b','saln',1);   
%-- or --
%[fld,lon , lat, depths]=loadrestart('N32DAILY_1958_000_1958_000.b','saln',1,1);
%fld=loadrestart('N32DAILY_1958_000_1958_000.b','saln',1,1); % Skips lon/lat/depths
%
%To read temperature for layer 1 to 3 from the file  N32DAILY_1958_000_1958_000.[ab]:
%[fld,lon , lat, depths]=loadrestart('N32DAILY_1958_000_1958_000.b','temp',1,3);
% 
%To read SSH (2-Dimensional variable) from the file  N32DAILY_1958_000_1958_000.[ab]:
%[fld,lon , lat, depths]=loadrestart('N32DAILY_1958_000_1958_000.b','ssh',1);
%
% There is also a list mode, which shows the variables in the file, example:
%[fld,lon , lat, depths]=loadrestart('N32DAILY_1958_000_1958_000.b','list');
%
%
%Knut Liseter, 17.08.2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lon=[];
lat=[];
depths=[];
fld=[];
onem=9806;
listmode=0;

if (nargin==2 & strcmp(varname,'list'))
   listmode=1;
   layer1=1;
   layer2=1;
elseif (nargin==3)
   layer2=layer1;
elseif (nargin~=4  )
   disp(['loadpak needs 3 or 4 input  arguments']);
   %disp(['loadpak needs      4 output arguments']);
   return
end

% layer index
layerind=layer1:1:layer2;
layermatch=zeros(prod(size(layerind)),1);
layerind_file=zeros(prod(size(layerind)),1);



% Convert name if necessary
i=findstr(rstfile,'.a')-1;
if (~isempty(i))
   rstfile=rstfile(1:i);
end
i=findstr(rstfile,'.b')-1;
if (~isempty(i))
   rstfile=rstfile(1:i);
end

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




%disp(['Opening header file ' rstfile '.b']);
fid=fopen([rstfile '.b']);
if (fid~=-1)


   for k=1:2
      A=fgetl(fid);
   end

   aindex=1;
   match=0;
   %Read until match
   while (~match | listmode)
      

      A=fgetl(fid);
      if (A==-1) 
         if (~listmode)
            disp([ 'Could not find wanted variable ' varname ])
            disp([ 'Use List mode to get list of variables'  ])
         end
         return
      end

      fldid=A(1:8);
      fldidtmp=fldid;
      fldid=strtrim(fldid);
      i=findstr(A,'=')+1;
      tst=sscanf(A(i:prod(size(A))),'%f');
      layer=tst(1);
      timeind=tst(2);
      bmin=tst(3);
      bmax=tst(4);
      %disp([ fldid ' ' num2str(layer) ' ' num2str(match) ]);

      if (listmode)
         S=[ 'Variable name :' fldidtmp ' layer: ' num2str(layer) ];
         disp(S);
      end

      % Get field(s)
      %disp([ fldid ' ' varname ]);
      if (strcmp(fldid,varname))
         ind=layer-layer1+1;
         if (ind>0 & ind <= layer2-layer1+1)
            %disp(['var:' num2str(ind) ' ' num2str(layer) ' ' fldid]);
            layerind_file(ind)=aindex;
            layermatch(ind)=1;
         end
         if (sum(layermatch)==prod(size(layermatch)))
            match=1;
         end 
      end

      aindex=aindex+1;
   end

   fclose(fid);
else
   disp(['Error - could not read header file - I quit'])
   return
end

%disp(['Opening data file ' rstfile '.a' ]);
% A quirk of the mod_za module, data is dumped containing a whole multiple of 4096 values.
% This is used to skip to the correct record, but we will only read idm*jdm values...
n2drec=floor((idm*jdm+4095)/4096)*4096;
bytes_per_float=4;
fid=fopen([rstfile '.a'],'r','ieee-be'); % Big-endian

%layerind_file

% Skip to indices in layerind_file
fld=zeros(prod(size(layerind_file)),idm,jdm);
for i=1:prod(size(layerind_file))
   stat=fseek(fid,n2drec*bytes_per_float*(layerind_file(i)-1),'bof'); % Skip to corr record
   fldtmp=fread(fid,[idm jdm],'single');
   %size(fldtmp)
   %imagesc(fldtmp./avecount)

   fld(i,:,:)=fldtmp;
end
fclose(fid);


if (size(fld,1)==1) 
   fld=reshape(fld,size(fld,2),size(fld,3));
end




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
      disp(['newpos.uf not found -- lon lat will be empty']);
      lon=[];
      lat=[];
   end
end

if (nargout>1)
   fdepths=['depths' num2str(idm,'%3.3i') 'x' num2str(idm,'%3.3i') '.uf'];
   fid=fopen(['depths' num2str(idm,'%3.3i') 'x' num2str(jdm,'%3.3i') '.uf'], ...
             'r','ieee-be');
   if (fid~=-1)
      stat=fseek(fid,4,'bof'); % Skip fortran 4-byte header 
      depths=fread(fid,[idm jdm],'double');
      fclose(fid);
   else
      disp([fdepths ' not found -- depths will be empty']);
      depths=[];
   end
end

