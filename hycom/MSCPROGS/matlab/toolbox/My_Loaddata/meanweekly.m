function [avefld,lon,lat,depths]=meanweekly(files,varname,layer);
nfld=0;
for i=1:prod(size(files))
   disp( [ num2str(i)  ' '  num2str(prod(size(files))) ] );
   if (i==1) 
      [fld, lon , lat, depths]=loadweekly(files(i).name,varname,layer);
      avefld=fld;
   else
      fld=loadweekly(files(i).name,varname,layer);
      avefld=avefld+fld;
   end
   nfld=nfld+1;
end

avefld=avefld/nfld;



