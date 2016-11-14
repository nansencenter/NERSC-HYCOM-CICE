function writebathy(depths,prefix);
%function writebathy(depths,prefix);
%
% Routine to write a bathymetry modified by matlab to a new file. This routine 
% writes  three files depthsXXXxXXX.uf.new , regional.depth.a.new and regional.depth.b.new
%
% If optional argument pre is given, files will be prefixed with that. Otherwise they will
% be prefixed with ``new''
idm=size(depths,1);
jdm=size(depths,2);

if (nargin==2)
   pre=prefix;
else
   pre='new';
end



% This opens the unformatted file type "depthsXXXxXXX.uf" - sequential access
if (idm>999 | jdm>999) 
   fnameuf=[pre '.depths' num2str(idm,'%4.4d') 'x'  num2str(jdm,'%4.4d') '.uf' ];
else
   fnameuf=[pre '.depths' num2str(idm,'%3.3d') 'x'  num2str(jdm,'%3.3d') '.uf' ];
end
fnamebase=[ pre '.regional.depth.a'];

fid=fopen(fnameuf,'w','ieee-be');
fwrite(fid,8*idm*jdm,'integer*4'); % Header - Total bytes to be written
fwrite(fid,depths,'double');
fwrite(fid,8*idm*jdm,'integer*4'); % "Trailer" - what were these Fortran guys thinking?
fclose(fid);
disp(['depths dumped to ' fnameuf ]);


% This opens .b and .a (direct access) - the latter is easier to dump to
fida=fopen([ fnamebase '.a'],'w','ieee-be');
n2drec=floor((idm*jdm+4095)/4096)*4096;
npad=n2drec-idm*jdm;
fwrite(fida,depths,'float');
if (npad>0)
   fwrite(fida,zeros(npad,1),'float');
end
fclose(fida);

fidb=fopen([fnamebase '.b'],'w');
fprintf(fidb,'%s\n\n\n\n\n','Bathymetry dumped by matlab');
fprintf(fidb,'min,max depth =     %8.3f %8.3f',min(min(depths)) , max(max(depths)));
fclose(fidb);

disp(['depths dumped to ' fnamebase '[.ab]']);


