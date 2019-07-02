function showconfmap(gridroutine,zfac)
%function showconfmap(gridroutine,zfac)
%
%Routine for running confgrid routine  and displaying its location
%on a stereographic map (using m_map).
%
% First argument is the full path to the to the conformal mapping routine
% Second argument is a scale factor. < 1 to zoom in, > 1 to zoom out


% Run grid routine
if (nargin>=1)
   eval(['!' gridroutine ]);
end

if (nargin==2)
   zoomfac=zfac;
else
   zoomfac=1.2;
end

global MAP_PROJECTION;
warning off MATLAB:FINITE:obsoleteFunction;

% Get grid size - 1st line of latlon.fat
A=textread('latlon.dat','%d',2);
idm=A(1);
jdm=A(2);
disp(['Grid size ' num2str(idm) ' ' num2str(jdm)])



% This opens the unformatted file type "newpos.uf" - sequential access
fname='newpos.uf';
fid=fopen(fname,'r','ieee-be');
if (fid~=-1)
   nent=fread(fid,1,'integer*4');  % Header - Total read to be written (here we skip it)
   plat=fread(fid,[ idm jdm],'double');
   plon=fread(fid,[ idm jdm],'double');
   fclose(fid);
else
   disp('Could not open newpos.uf')
end

centerlon=plon(round(idm/2),round(jdm/2));
centerlat=plat(round(idm/2),round(jdm/2));
disp(['Center ' num2str(centerlon) ' ' num2str(centerlat)])

% Calculate approximate max degree between center and boundary
[s] = m_idist(centerlon,centerlat,plon,plat);
maxdist=max(max(s));
maxdeg=maxdist/100000;
maxdeg=max(maxdeg,20);
disp(['Projection radius in degrees ' num2str(maxdeg)])


figure(1) ; clf
% This would make user-defined map_projections possible - turned off for now
%MAP_PROJECTION
%if isempty(MAP_PROJECTION) ; 
   m_proj('stereographic','lon',centerlon,'lat',centerlat ,'rad',maxdeg*zoomfac,'rec','on')
%end
m_coast('LineWidth',2,'Color','k');
hold on;

for i=idm/10:idm/10:idm*9/10
   m_plot(plon(floor(i),:),plat(floor(i),:),'LineWidth',1,'Color','r');
end
for j=jdm/10:jdm/10:jdm*9/10
   m_plot(plon(:,floor(j)),plat(:,floor(j)),'LineWidth',1,'Color','r');
end


m_plot(plon(1,:),plat(1,:),'LineWidth',2,'Color','r');
m_plot(plon(:,1),plat(:,1),'LineWidth',2,'Color','r');
m_plot(plon(:,jdm),plat(:,jdm),'LineWidth',2,'Color','r');
m_plot(plon(idm,:),plat(idm,:),'LineWidth',2,'Color','r');

m_grid
m_elev();
