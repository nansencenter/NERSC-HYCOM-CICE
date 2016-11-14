function [fld,lon,lat,depths]=loadweekly(weeklyfile,varname,layer1,layer2);
%function [fld,lon,lat,depths]=loadweekly(weeklyfile,varname,layer1,layer2);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Routine loadweekly:
%%%%%%%%%%%%%%%%%
%This routine reads data from 4-byte big-endian weekly average files (.ab - files)
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
%There is also a list mode, called in this way:
%function [fld,lon,lat,depths]=loadweekly(weeklyfile,'list');
%
%Examples: 
% 
%To read Salinity for layer 1 from the file  N31AVE_1958_08_2.[ab], two approaches:
%[fld, lon , lat, depths]=loadweekly('N31AVE_1958_08_2.b','saln',1);   
%-- or --
%[fld, lon , lat, depths]=loadweekly('N31AVE_1958_08_2.b','saln',1,1);
%
%To read temperature for layer 1 to 3 from the file  N31AVE_1958_08_2.[ab]:
%[fld, lon , lat, depths]=loadweekly('N31AVE_1958_08_2.a','temp',1,3);
% 
%To read SSH (2-Dimensional variable) from the file  N31AVE_1958_08_2.[ab]:
%[fld, lon , lat, depths]=loadweekly('N31AVE_1958_08_2.a','ssh',0);
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
   layer1=1; layer2=1;
elseif (nargin<4)
   layer2=layer1;
elseif (nargin<3 | nargin>4  )
   disp(['loadpak needs 1 3 or 4 input  arguments']);
   %disp(['loadpak needs      4 output arguments']);
   return
end

% layer index
layerind=layer1:1:layer2;
layermatch=zeros(prod(size(layerind)),1);
layerind_file=zeros(prod(size(layerind)),1);
presmatch=zeros(prod(size(layerind)),1);
presind_file=zeros(prod(size(layerind)),1);



% Convert name if necessary
i=findstr(weeklyfile,'.a')-1;
if (~isempty(i))
   weeklyfile=weeklyfile(1:i);
end
i=findstr(weeklyfile,'.b')-1;
if (~isempty(i))
   weeklyfile=weeklyfile(1:i);
end

%disp(['Opening header file ' weeklyfile '.b']);
fid=fopen([weeklyfile '.b']);
if (fid~=-1)


   A=fgetl(fid);
   A=fgetl(fid);
   A=fgetl(fid);
   A=fgetl(fid);
   A=fgetl(fid);
   A=fgetl(fid);

   %Yearflag
   A=fgetl(fid);
   yrflag=sscanf(A,'%d');

   %idm
   A=fgetl(fid);
   idm=sscanf(A,'%d');

   %jdm
   A=fgetl(fid);
   jdm=sscanf(A,'%d');

   %kdm
   A=fgetl(fid);
   kdm=sscanf(A,'%d');

   A=fgetl(fid);
   A=fgetl(fid);

   %AVE counter
   A=fgetl(fid);
   avecount=sscanf(A,'%d');

   A=fgetl(fid);
   
   aindex=1;
   match=0;

   if (strcmp(varname,'saln') | strcmp(varname,'temp')    | ...  
       strcmp(varname,'utot') | strcmp(varname,'vtot')    | ... 
       strcmp(varname,'pres') | strcmp(varname,'kinetic') | ...
       strcmp(varname,'nit')  | strcmp(varname,'pho' )    | ...
       strcmp(varname,'sil')  | strcmp(varname,'dia' )    | ...
       strcmp(varname,'fla')  | strcmp(varname,'oxy' ) )

      dpmatch=0;
      var3d=1;
   else
      dpmatch=1;
      var3d=0;
   end 


   %Read until match
   while (~match | ~dpmatch | listmode )
      

      A=fgetl(fid);
      if (A==-1) 
         if (~listmode)
            disp([ 'Could not find wanted variable ' varname ])
            disp([ 'Use List mode to get list of variables'  ])
            disp([ 'Also note that layer thickness = "pres"' ])
         end
         return
      end

      fldid=A(1:8);
      fldidtmp=fldid;
      fldid=strtrim(fldid);
      A=A(11:prod(size(A)));
      tst=sscanf(A,'%f');
      nstep=tst(1);
      dtime=tst(2);
      layer=tst(3);
      dens=tst(4);
      bmin=tst(5);
      bmax=tst(6);
      %disp([ fldid ' ' num2str(layer) ' ' num2str(match) ' ' num2str(dpmatch)]);

      if (listmode)
         S=[ 'Variable name :' fldidtmp ' layer: ' num2str(layer) ];
         disp(S);
      end

      % Get field(s)
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

      % Get dp fields as well (used for ocean vars)
      % NB - a mod-average blooper here, pres is actually "dp"
      if (strcmp(fldid,'pres'))
         ind=layer-layer1+1;
         if (ind>0 & ind <= layer2-layer1+1)
            %disp([ num2str(ind) ' ' num2str(layer) ]);
            presind_file(ind)=aindex;
            presmatch(ind)=1;
         end
         if (sum(presmatch)==prod(size(presmatch)))
            dpmatch=1;
         end 
      end


      aindex=aindex+1;
   end

   fclose(fid);
else
   fid
   disp(['Error - could not read header file ' weeklyfile '.b' ' - I quit'])
   return
end


%disp(['Opening data file ' weeklyfile '.a' ]);
% A quirk of the mod_za module, data is dumped containing a whole multiple of 4096 values.
% This is used to skip to the correct record, but we will only read idm*jdm values...
n2drec=floor((idm*jdm+4095)/4096)*4096;
bytes_per_float=4;
fid=fopen([weeklyfile '.a'],'r','ieee-be'); % Big-endian
if (fid==-1)
   disp(['Error - could not open data file - I quit'])
   return
end

%layerind_file
%presind_file

% Skip to indices in layerind_file
fld=zeros(prod(size(layerind_file)),idm,jdm);
for i=1:prod(size(layerind_file))
   stat=fseek(fid,n2drec*bytes_per_float*(layerind_file(i)-1),'bof'); % Skip to corr record
   fldtmp=fread(fid,[idm jdm],'single');
   %imagesc(fldtmp./avecount)

   % Seek layer thickness
   if (var3d)
      stat=fseek(fid,n2drec*bytes_per_float*(presind_file(i)-1),'bof'); % Skip to corr record
      dp=fread(fid,[idm jdm],'single');

      %imagesc(fldtmp./(dp +1))
      if (~strcmp(varname,'pres'))
         %fld(i,:,:)=fldtmp./(dp+1); ! NBNB ! usually in meters
         fld(i,:,:)=fldtmp./(dp+1e-4);
      else
         fld(i,:,:)=fldtmp./(avecount);
      end
   else
      fld(i,:,:)=fldtmp./avecount;
   end 

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

   % try to retrieve from regional.grid.a
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
   fdepths=['depths' num2str(idm,'%3.3i') 'x' num2str(idm,'%3.3i') '.uf'];
   fid=fopen(['depths' num2str(idm,'%3.3i') 'x' num2str(jdm,'%3.3i') '.uf'], ...
             'r','ieee-be');
   if (fid~=-1)
      stat=fseek(fid,4,'bof'); % Skip fortran 4-byte header 
      depths=fread(fid,[idm jdm],'double');
      fclose(fid);
   else
      try
         [depths]=loada('regional.depth.a',1,idm,jdm);
      catch
         disp([fdepths ' or regional.depth.a not found -- depths will be empty']);
         depths=[];
      end 
   end
end

