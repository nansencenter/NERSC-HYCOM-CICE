%Program that create the lon lat value from the IBCAO projection input
% Francois Counillon 07/05/2013

nx=11617    %model dimmension in x
ny=11617    %model dimension in y
e=0.0818191908426215  %eccentricity
t_scale=75  %true scale projection
c_mer=0   %central meridien
radius=6378137; %Radius of the earth
x_range=[-2904000 2904000] %Range in X of the stereographic projection
y_range=[-2904000 2904000] %Range in Y of the stereographic projection
dx=(x_range(2)-x_range(1))/(nx-1);
dy=(y_range(2)-y_range(1))/(ny-1);
x=x_range(1):dx:x_range(2);
%y=y_range(2):-dy:y_range(1);
y=y_range(1):dy:y_range(2);
fid=fopen('Lat_ibcao.uf','w','ieee-be');
fid2=fopen('Lon_ibcao.uf','w','ieee-be');
fwrite(fid,8*nx*ny,'integer*4'); % Header - Total bytes to be written
fwrite(fid2,8*nx*ny,'integer*4'); % Header - Total bytes to be written
for j=1:ny
 for i=1:nx
      [lat(i),lon(i)]=polarstereo_inv(x(i),y(j),radius,e,t_scale,c_mer);
 end
 fwrite(fid2,lon,'double');
 fwrite(fid ,lat,'double');
 j
end
fwrite(fid,8*nx*ny,'integer*4'); % "Trailer" - what were these Fortran guys thinking?
fwrite(fid2,8*nx*ny,'integer*4'); % "Trailer" - what were these Fortran guys thinking?
fclose(fid);
fclose(fid2);

