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

