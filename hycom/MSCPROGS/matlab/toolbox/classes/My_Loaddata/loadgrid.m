function [lon,lat,depths]=loadgrid();
%function [lon,lat,depths]=loadgrid();
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Routine loadgrid:
%%%%%%%%%%%%%%%%%
% Reads data from typical grid files (newpos.uf, depths....uf regional...)
% and returns longitude latitude (in p-points) and depths
%
%
%Knut Liseter, 17.08.2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lon=[];
lat=[];
depths=[];

files=dir('depths*.uf');
if isempty(files) 
   % Get grid dimensions from regional
   try
      fid=fopen('regional.grid.b');
      A=fgetl(fid); idm=sscanf(A,'%f');
      A=fgetl(fid); jdm=sscanf(A,'%f');
      fclose(fid)
   catch
      files=dir('depths*.uf')
      disp('can not get grid size');
      return
   end
else
   if (prod(size(files(1).name))==16) 
      idm=str2num(files(1).name(7:9));
      jdm=str2num(files(1).name(11:13));
   elseif (prod(size(files(1).name))==18) 
      idm=str2num(files(1).name(7:10));
      jdm=str2num(files(1).name(12:15));
   else
      disp('can not get grid size');
      return
   end
end




% Try to retrieve lon/lat from newpos.uf
fid=fopen('newpos.uf','r','ieee-be');
if (fid~=-1)
   try
      stat=fseek(fid,4,'bof'); % Skip fortran 4-byte header 
      lat=fread(fid,[idm jdm],'double');
      lon=fread(fid,[idm jdm],'double');
      %contourf(lon)
      fclose(fid);
   catch
      disp(['could not read newpos.uf -- lon lat will be empty']);
      lon=[];
      lat=[];
   end

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

fdepths=['depths' num2str(idm,'%3.3i') 'x' num2str(idm,'%3.3i') '.uf'];
fid=fopen(['depths' num2str(idm,'%3.3i') 'x' num2str(jdm,'%3.3i') '.uf'], ...
          'r','ieee-be');
if (fid~=-1)
   try
      stat=fseek(fid,4,'bof'); % Skip fortran 4-byte header 
      depths=fread(fid,[idm jdm],'double');
      fclose(fid);
   catch
      disp([fdepths ' could not read depths file -- depths will be empty']);
      depths=[];
   end
else
   try
      [depths]=loada('regional.depth.a',1,idm,jdm);
   catch
      disp([fdepths ' or regional.depth.a not found -- depths will be empty']);
      depths=[];
   end 
end
end
