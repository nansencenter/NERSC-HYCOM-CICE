function [lon,lat,depths,fld]=loadpak(pakfile,varname,layer1,layer2);
%function [lon,lat,depths,fld]=loadpak(pakfile,varname,layer1,layer2);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Routine loadpak:
%%%%%%%%%%%%%%%%%
%This routine reads pak-files and converts from 6-bit ascii-encoded
%numbers to floats. The filename and the name of the variable is supplied
%as arguments. Also, the user must supply a layer argument to the routine. 
%For 3D vars, the layer is the model layer. For 2D vars, the layer must
%be equal to 1. You can also specify a range of layers by specifying 
%the first and last layer to be read.
%
%The routine also tries to read the depths from the model depths-file
%and lon/lat from the file newpos.uf. If these files are absent, the
%fields lon,lat and/or depths will be empty.
%
%Examples: 
% 
%To read Salinity for layer 1 from the file  FORy1994d300h00, two approaches:
%[lon , lat, depths, fld]=loadpak('FORy1994d300h00','SAL',1);   
%-- or --
%[lon , lat, depths, fld]=loadpak('FORy1994d300h00','SAL',1,1);
%
%To read temperature for layer 1 to 3 from the file  FORy1994d300h00:
%[lon , lat, depths, fld]=loadpak('FORy1994d300h00','TEM',1,3);
% 
%To read SSH (2-Dimensional variable) from the file  FORy1994d300h00:
%[lon , lat, depths, fld]=loadpak('FORy1994d300h00','SSH',1);
%
%
%Note that there are some sanity checks on the arguments supplied to 
%this routine. Also, if you have a header-file accompanying the pak-file, 
%its contents will be checked to see wether the field you specify (varname) 
%corresponds to a field in the pakfile. If not, a list of possible 2D and 
%3D fields will be displayed and the routine will exit.
%
%In case you are missing a header file, loadpak will do its best to read 
%the fields specified. This means all sanity checks on the arguments
%to this routine is up to you, the user ...
%
%Knut Liseter, 24.09.2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lon=[];
lat=[];
depths=[];
fld=[];

if (nargin<4)
   layer2=layer1;
elseif (nargin<3 | nargin>4 | nargout ~=4  )
   disp(['loadpak needs 3 or 4 input  arguments']);
   disp(['loadpak needs      4 output arguments']);
   return
end


disp(['Opening header file ' pakfile '.hdr']);
fid=fopen([pakfile '.hdr']);
if (fid~=-1)

   % Some potentially useful stuff here - for now its skipped
   A=fscanf(fid,'%f%d');
   ver=A(1);
   nlines=A(2);
   rungen=fscanf(fid,'%s%s',2);
   fname=fscanf(fid,'%s',1);
   tmp=fscanf(fid,'%s',2);

   %2D flds
   n2dfld=fscanf(fid,'%d',1);
   is2d=0;
   for i=1:n2dfld
      tmp=fscanf(fid,'%s',1);
      flds2D(i)=cellstr(tmp);
      is2d=is2d | strcmp(strtok(tmp),varname);
   end

   %3D flds
   tmp=fscanf(fid,'%s',2);
   n3dfld=fscanf(fid,'%d',1);
   is3d=0;
   for i=1:n3dfld
      tmp=fscanf(fid,'%s',1);
      flds3D(i)=cellstr(tmp);
      is3d=is3d | strcmp(strtok(tmp),varname);
   end

   %Domain size
   tmp=fscanf(fid,'%31c',1);
   %tmp(13:31)
   idm=sscanf(tmp(13:17),'%d');
   %tmp(20:31)
   jdm=sscanf(tmp(20:23),'%d');
   %tmp(27:31)
   kdm=sscanf(tmp(27:31),'%d');
   
   %return
   %tmp=fscanf(fid,'%s',2);
   %idm=fscanf(fid,'%d',1);
   %tmp=fscanf(fid,'%s',1);
   %jdm=fscanf(fid,'%d',1);
   %tmp=fscanf(fid,'%s',1);
   %kdm=fscanf(fid,'%d',1);

   % Finished reading header (remander is skipped)
   fclose(fid);

   if (~is2d & ~is3d)
      disp(['Variable ' varname ' is neither a 2D or 3D var' ]);
      disp('Possible variable names are:')
      disp('2D:');
      for i=1:n2dfld
         disp(char(flds2D(i)));
      end
      disp(' ');
      disp('3D:');
      for i=1:n3dfld
         disp(char(flds3D(i)));
      end
      return
   end

   % If a 3d Var, make sure that layers make sense
   if (is3d) 
      if (layer1>layer2)
         disp(['layer1 > layer2 ... switching ']);
         tmp=layer1;
         layer1=layer2;
         layer2=tmp;
      end 

      if (layer1 < 1 | layer2 < 1 )
         disp([ 'layers < 1 , setting layers to  max(1,layer)']);
         layer1=max(layer1,1);
         layer2=max(layer2,1);
      end 

      if (layer1 > kdm | layer2 > kdm )
         disp([ 'layers > kdm , setting layers to  min(kdm,layer)']);
         layer1=min(layer1,kdm);
         layer2=min(layer2,kdm);
      end 
   end

   if (is2d)
      if (layer1~=layer2 | layer2~=1)
         disp([ 'For 2D plots layers should be equal to 1 .. fixing this ']);
         layer1=1;
         layer2=1;
      end
   end 
else
   disp(['Error - could not read header file - Guessing from now on...'])
end

disp(['Opening data file ' pakfile ]);

%Open pakfile
fid=fopen(pakfile);
if (fid==-1)
   disp(['Could not open file ' pakfile]);
   return
end

ilayer=layer1;
foundvar=0;
fld=[];
while (~foundvar)

   % Read field header
   fieldhead=fscanf(fid,'%80c',1);
   [vname,readlayer,idm,jdm]=invfldh(fieldhead);

   % Read data
   foundvar=strcmp(varname,strtok(vname)) & (readlayer==ilayer) ;

   if (foundvar)
      %disp(['hey'])


      % Read data
      util=fscanf(fid,'%2c',(idm*jdm+14));

      % Convert to floats
      [base,scale,vals2]=unpakk(util,jdm,idm);


      % Put into field "fld"
      if(~isempty(fld))
         fld=cat(3,fld,vals2);
      else
         fld=vals2;
      end
      %size(vals2)
      %disp([num2str(ilayer)]);
      %size(fld)


      disp([ 'got variable ' varname ' for layer ' num2str(ilayer)]);
      ilayer=ilayer+1; % Next layer to seek
      foundvar=ilayer>layer2; % Unset if no other layers are needed
         
   else
      %Skip this number of bytes if field (from fieldheader) not used
      status=fseek(fid,(idm*jdm+14)*2,0); 
   end 
end



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



lon=lon';
lat=lat';
depths=depths';
if (~isempty(depths))
   I=find(depths<1.);
   for i=1:size(fld,3);
      tmp=fld(:,:,i);
      tmp(I)=nan;
      fld(:,:,i)=tmp;
   end
end





function [vname,readlayer,idm,jdm]=invfldh(fieldhead)
   vname=sscanf(fieldhead(4:8),'%5c');
   readlayer=sscanf(fieldhead(9:10),'%d');
   ix=sscanf(fieldhead(11:12),'%d');
   iy=sscanf(fieldhead(13:14),'%d');
   ia=sscanf(fieldhead(15:18),'%d');
   ib=sscanf(fieldhead(19:22),'%d');
   ic=sscanf(fieldhead(23:26),'%d');
   length=sscanf(fieldhead(27:80),'%d');
   idm=ib;
   jdm=ic;


function [base,scale,vals2]=unpakk(util,jdm,idm);
   base=sscanf(util(1:14),'%f');
   scale=sscanf(util(15:28),'%f');
   values=util(29:2*(idm*jdm+14));

   % Rearrange 
   values2=reshape(values,2,jdm,idm);
   
   vals=values2-0; % Trick to turn chars into ascii codes
   I=find(vals>96);
   vals(I)=vals(I)-6;
   I=find(vals>64);
   vals(I)=vals(I)-7;
   vals=vals-14;

   vals2=scale*(64*(vals(1,:,:)-32)+(vals(2,:,:)-32))+base;
   vals2=reshape(vals2,jdm,idm);
   %size(vals2)


