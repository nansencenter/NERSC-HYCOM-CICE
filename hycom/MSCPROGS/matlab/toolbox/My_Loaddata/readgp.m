%function [z,data,time] =readgp(mon,var,region,steps);
%function [z,data,time] =readgp(mon,var,region,steps);



% Read record file
file='T06gprecord.asc';
fid=fopen(file)
irec=0;
line=fgetl(fid); % skip first line
while isstr(line)
   line=fgetl(fid);
   if isstr(line)
      irec=irec+1;
      ind=min(min(strfind(line,' ')));
      tmp=strtrim(line(1:ind));
      vname(irec)=cellstr(tmp);
      %tmp=sscanf(line(ind:prod(size(line))),'%f');
      vdim(irec)=sscanf(line(ind:prod(size(line))),'%f');
   end
end
nrec=irec;
fclose(fid);

listvars=1==1;
if (listvars) 
   disp('Variables available:')
   disp('dimension name      ')
   for irec=1:nrec
      disp([ num2str(vdim(irec),'%2.2d')  '  ----  ' char(vname(irec))]);
   end 
end

% Three first recs are integers, the rest are real*4
nentries=sum(sum(vdim(1:nrec)));
nentriesint =3;
nentriesreal=nentries-nentriesint;;

file='gp_1968_09_+000.064x+070.610.uf';
fid=fopen(file,'r','ieee-be'); % big-endian
[tmp1]=fread(fid,nentriesint ,'integer*4');
[tmp2]=fread(fid,nentriesreal,'real*4');
fclose(fid);

% Shape the resulting data into arrays with names from vname
% These are fairly fixed
itime=tmp1(1);
ipiv=tmp1(2);
jpiv=tmp1(3);

cnt=1;
for irec=3:nrec % Starting from first float
   startrec(irec)=cnt
   cnt=cnt+vdim(irec);
   endrec  (irec)=cnt-1;
end







%mont=mod((mon-1),12)+1;
%y=floor((mon-1)/12);
%
%rungen='NW4GP';
%workd='/work/ceciha/NWSNO/';
%lon=[2.001];
%lat=[66.013];
%reglo=lon(region);      
%regla=lat(region);      
%
%if(reglo<0)
%   reg=['_-00',num2str(abs(reglo),'%3.3f'),'x+0',num2str(regla,'%3.3f'),'.uf'];
%elseif(reglo>0)
%   reg=['_+00',num2str(reglo,'%3.3f'),'x+0',num2str(regla,'%3.3f'),'.uf'];
%end
%
%dirc = [workd rungen num2str(y,'%4.4i') '_'];
%file=[dirc num2str(mont,'%2.2i') ...
%      '/gp_' num2str(y,'%4.4i') '_' ...
%      num2str(mont,'%2.2i') reg];
%
%%eval(['ls -l ' file])
%fid=fopen(file,'r','ieee-be'); 
%iostat=feof(fid);
%
%fstep=0;
%cnt=0;
%NUL=23; %number of layers
%
%while (iostat==0)
%      cnt=cnt+1;
%      time(cnt)=fread(fid,1,'int32');
%      i=fread(fid,2,'int32');
%      %a=fread(fid,23*5,'float32');
%      Z(cnt,:)=fread(fid,23,'float32');
%      U(cnt,:)=fread(fid,23,'float32');
%      V(cnt,:)=fread(fid,23,'float32');
%      TEMP(cnt,:)=fread(fid,23,'float32');
%      SALN(cnt,:)=fread(fid,23,'float32');
%      a=fread(fid,3,'float32');
%      NIT(cnt,:)=fread(fid,23,'float32');
%      PHO(cnt,:)=fread(fid,23,'float32');
%      SIL(cnt,:)=fread(fid,23,'float32');
%      DET(cnt,:)=fread(fid,23,'float32');
%      SIS(cnt,:)=fread(fid,23,'float32');
%      FLA(cnt,:)=fread(fid,23,'float32');
%      DIA(cnt,:)=fread(fid,23,'float32');
%      OXY(cnt,:)=fread(fid,23,'float32');
%      SED(cnt,:)=fread(fid,23,'float32');
%      YEL(cnt,:)=fread(fid,23,'float32');
%      iostat=feof(fid);
%      if(cnt>=steps)
%         iostat=1;
%      end
%end
%fclose(fid);
%
%
%   z=Z';
%if(strcmp(var,'uvel'))   
%   data=U';
%elseif(strcmp(var,'vvel'))
%   data=V';
%elseif(strcmp(var,'temp'))
%   data=TEMP';
%elseif(strcmp(var,'saln'))
%   data=SALN';
%elseif(strcmp(var,'dens'))
%   [alph,rho]=eos80(0,TEMP',SALN');
%   data=rho;
%elseif(strcmp(var,'nit'))
%   data=NIT';
%elseif(strcmp(var,'pho'))
%   data=PHO';
%elseif(strcmp(var,'sil'))
%   data=SIL';
%elseif(strcmp(var,'det'))
%   data=DET';
%elseif(strcmp(var,'sis'))
%   data=SIS';
%elseif(strcmp(var,'fla'))
%   data=FLA';
%elseif(strcmp(var,'dia'))
%   data=DIA';
%elseif(strcmp(var,'oxy'))
%   data=OXY';
%elseif(strcmp(var,'sed'))
%   data=SED';
%elseif(strcmp(var,'yel'))
%   data=YEL';
%end
%%
