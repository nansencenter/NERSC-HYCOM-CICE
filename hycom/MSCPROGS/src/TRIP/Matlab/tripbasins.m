load trips.tec;
load trips_rivers.tec;

% Put trips data into a matrix
nx=max(trips(:,1));
ny=max(trips(:,2));
lons=zeros(nx,ny);
lats=zeros(nx,ny);
basins=zeros(nx,ny);
dirnum=zeros(nx,ny);
for I=1:size(trips,1)
   i=trips(I,1);
   j=trips(I,2);
   lons  (i,j)=trips(I,3);
   lats  (i,j)=trips(I,4);
   basins(i,j)=trips(I,9);
   dirnum(i,j)=trips(I,8);
end

dx=lons(2,1)-lons(1,1)
dy=lats(1,1)-lats(1,2)



%% Find land points
%I=find(trips(:,8)~=0);
%
%% Plot Basins  - river outlets
figure(1) ; clf
pcolor(lons,lats,basins); shading flat; hold on;
tmparea=trips_rivers(:,7);
tmparea=max(15,tmparea/max(tmparea)*500);
S2=scatter(trips_rivers(:,3),trips_rivers(:,4),tmparea, ...
           trips_rivers(:,5),'filled')
set(S2,'MarkerEdgeColor','k')
colormap(colorcube)

figure(2); clf
pcolor(lons,lats,dirnum); shading flat; hold on;
I=find(trips(:,8)~=0);
quiver(trips(I,3)+dx/2,trips(I,4)+dy/2,trips(I,10),trips(I,11));
