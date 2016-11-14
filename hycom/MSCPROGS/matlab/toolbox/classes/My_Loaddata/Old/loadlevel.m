function [fld,lon,lat,depths]=loadlevel(file,type,varname,level);
%function [fld,lon,lat,depths]=loadlevel(file,type,varname,level);
% Use loadlevel2 in stead

onem=9806;

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

% To get lon, lat, depths
if(strcmp(type,'restart'))
   [fldin, lon , lat, depths]=loadrestart(file,varname,1);
elseif(strcmp(type,'weekly'))
   [fldin, lon , lat, depths]=loadweekly(file,varname,1);
elseif(strcmp(type,'daily'))
   [fldin, lon , lat, depths]=loaddaily(file,varname,1);
else
   disp(['Uknown filetype ' type ]);
   return
end

for k=1:kdm

   disp(['Reading layer ' num2str(k) ' of ' num2str(kdm)])

   if(strcmp(type,'restart'))
      fldin=loadrestart(file,varname,k);
      dp   =loadrestart(file,'dp',k);
      dp=dp/onem;
   elseif(strcmp(type,'weekly'))
      fldin=loadweekly(file,varname,k);
      dp   =loadweekly(file,'pres',k);
   elseif(strcmp(type,'daily'))
      fldin=loaddaily(file,varname,k);
      dp   =loaddaily(file,'pres',k);
      if (k>1) 
         [intf , lon , lat, depths]=loaddaily(file,'pres',k-1);
      else
         intf=zeros(size(dp));
      end
      dp=dp-intf;
      dp=dp/onem;
   end

   if (k==1) 
      I=find(dp>.1);
      disp([ 'Thickness of 1st layer: ' num2str(min(min(dp(I)))) ]);
      disp([ '(If this is very high p is in pressure coords = modify level!)']);
   end



   if (k==1)
      p=dp;
      fld=zeros(size(p));
      nofill=ones(size(p)); % Where not to fill with values ...
   else
      p=p+dp;
   end


   I=find(nofill==1); % Points which havent been filled yet
   J=find(p>level);   % Points where lower interface is below level
   K=intersect(I,J);  % Intersection is where to fill with values

   fld   (K)=fldin(K); % These are valid values...
   nofill(K)=0;        % Makes sure we wont overwrite already valid values later on
end

K=find(nofill==1);
fld(K)=nan;
