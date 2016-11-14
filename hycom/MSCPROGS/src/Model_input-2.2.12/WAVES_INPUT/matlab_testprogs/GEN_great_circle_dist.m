function s  = GEN_great_circle_dist(lon1,lat1,lon2,lat2);
%%CALL: s  = GEN_great_circle_dist(lon1,lat1,lon2,lat2);
%% s = great circle distance [m];
%% lon1,lat1 are scalars
%% lon2,lat2 are matrices of any rank - s is the same size as them;

phi1  = pi/180*lat1;
phi2  = pi/180*lat2;
%%
lam1  = pi/180*lon1;
lam2  = pi/180*lon2;
%%
dphi  = pi/180*(lat2-lat1);
dlam  = pi/180*(lon2-lon1);
%%
r_earth  = 1e3*GEN_radius_earth();%m
s        = r_earth*acos(sin(phi1).*sin(phi2)+cos(phi1).*cos(phi2).*cos(dlam));

if s<=5e3
   s  = r_earth*2*asin(sqrt(sin(dphi/2).^2+cos(phi1).*cos(phi2).*sin(dlam/2).^2));
end
