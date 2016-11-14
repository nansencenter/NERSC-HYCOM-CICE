function fld=read_restart_ice(filename,varname,idm,jdm,record)
%function fld=read_restart_ice(filename,varname,idm,jdm)
% Read ice fields from hycom ICE restart files

% Size of one record entry (5 fields, double precision)
recsize=5*idm*jdm*8;

fid = fopen(filename,'r','ieee-be') ;
status=fseek(fid,(record-1)*recsize,-1);
ficem=fread(fid,[idm jdm],'double');
hicem=fread(fid,[idm jdm],'double');
hsnwm=fread(fid,[idm jdm],'double');
ticem=fread(fid,[idm jdm],'double');
tsrfm=fread(fid,[idm jdm],'double');
fclose(fid)

if(strcmp(varname,'ficem'))
   fld=ficem;
elseif(strcmp(varname,'hicem'))
   fld=hicem;
elseif(strcmp(varname,'hsnwm'))
   fld=hsnwm;
elseif(strcmp(varname,'ticem'))
   fld=ticem;
elseif(strcmp(varname,'tsrfm'))
   fld=tsrfm;
else
   disp([ 'Unknown field ' varname ])
end

