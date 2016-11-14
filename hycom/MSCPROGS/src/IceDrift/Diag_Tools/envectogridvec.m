function [ iu, jv, X, Y ] = envectogridvec(lon,lat,eu,nv)
%
%[ iu, jv ] = envectogridvec(lon,lat,eu,ev)
%
% evectogridvec=East/North VECtor TO GRID VECtor
%
% Routine to rotate velocities in east/north directions
% to the directions of the grid (increasing grid indices
% are positive)
%   lon  - longitude position of grid points
%   lat  - latitude  position of grid points
%   eu   - easterly  vector component
%   nu   - northerly vector component
%
%   iu   - grid vector component -- increasing 1st dimension (i)
%   iu   - grid vector component -- increasing 2nd dimension (j)



% Rotate velocities to easterly / westerly directions
nx=size(lon,1);
ny=size(lon,2);


radian=pi/180.;

% Longitude increment for one increment of 1st dimension
dlon_ip1=(lon(2:nx,:)-lon(1:nx-1,:));
%max(max(dlon_ip1))
%min(min(dlon_ip1))
I=find(dlon_ip1>180); dlon_ip1(I)=dlon_ip1(I)-360;
I=find(dlon_ip1<-180); dlon_ip1(I)=dlon_ip1(I)+360;
%figure(1); clf; pcolor(dlon_ip1) ; shading interp; colorbar
dlon_ip1=dlon_ip1.*cos( radian*.5*(lat(2:nx,:)+lat(1:nx-1,:)) ); %Rough approx
dlon_ip1(nx,:)=dlon_ip1(nx-1,:);

% Latitude increment for one increment of 1st dimension
dlat_ip1=lat(2:nx,:)-lat(1:nx-1,:);
dlat_ip1(nx,:)=dlat_ip1(nx-1,:);





% Longitude increment for one increment of 2nd dimension
dlon_jp1=(lon(:,2:ny)-lon(:,1:ny-1));
%max(max(dlon_jp1))
%min(min(dlon_jp1))
I=find(dlon_jp1> 180); dlon_jp1(I)=dlon_jp1(I)-360;
I=find(dlon_jp1<-180); dlon_jp1(I)=dlon_jp1(I)+360;
dlon_jp1=dlon_jp1.*cos( radian*.5*(lat(:,2:ny)+lat(:,1:ny-1)) ); %Rough approx
%figure(1); clf; pcolor(dlon_ip1) ; shading interp; colorbar
dlon_jp1(:,ny)=dlon_jp1(:,ny-1);

% Latitude increment for one increment of 1st dimension
dlat_jp1=lat(:,2:ny)-lat(:,1:ny-1);
dlat_jp1(:,ny)=dlat_jp1(:,ny-1);

% Right-handedness test
for i=1:nx
for j=1:ny
   %i=nx/2;
   %%disp([ num2str(theta_ip1(nx/2,j)) ' ' num2str(theta_jp1(nx/2,j))  ]);
   %disp(num2str(cross( [dlon_ip1(i,j) dlat_ip1(i,j) 0],[dlon_jp1(i,j) dlat_jp1(i,j) 0])))
   tmp=cross( [dlon_ip1(i,j) dlat_ip1(i,j) 0],[dlon_jp1(i,j) dlat_jp1(i,j) 0]);
   dirc(i,j)=tmp(3);
end
end

I=find(dirc<0.);
J=find(dirc>=0.);
disp([' # crprd>0 ' num2str(prod(size(J)))]);
disp([' # crprd<0 ' num2str(prod(size(I)))]);

% Matlab democratic vote 
if (prod(size(J))>prod(size(I))) 
   righth=1;
else
   righth=0;
end

righth

if (righth==1) 
   theta_ip1=-atan2(dlat_ip1,dlon_ip1); % Angle one dx displacement makes with constant latitude lines
   theta_jp1=-atan2(dlat_jp1,dlon_jp1); % Angle one dy displacement makes with constant latitude lines
   nv=-nv; % Hmm
else
   theta_ip1=atan2(dlat_ip1,dlon_ip1); % Angle one dx displacement makes with constant latitude lines
   theta_jp1=atan2(dlat_jp1,dlon_jp1); % Angle one dy displacement makes with constant latitude lines
   %theta_ip1=-atan2(dlon_ip1,dlat_ip1); % Angle one dx displacement makes with constant latitude lines
   %theta_jp1=-atan2(dlon_jp1,dlat_jp1); % Angle displacement makes with constant latitude lines
end 

%for j=1:ny
%disp([ num2str(theta_jp1(nx/2,j)*90) ' '  num2str(theta_ip1(nx/2,j)*90) ])
%end





%%Easterly
iu=eu.*cos(theta_ip1) + nv.*cos(theta_jp1);
jv=eu.*sin(theta_ip1) + nv.*sin(theta_jp1);

X=repmat([1:nx]',1,ny);
Y=repmat(1:ny,nx,1);

%figure(2); clf; pcolor(lon') ; shading interp; colorbar
%figure(3); clf; pcolor(theta_ip1') ; shading interp; colorbar;
%  cmap=colormap; cmap(30:35,:)=1; colormap(cmap)
%figure(4); clf; pcolor(theta_jp1') ; shading interp; colorbar;
%  cmap=colormap; cmap(30:35,:)=1; colormap(cmap)
%figure(5); quiver(X,Y,iu,jv)
