function Plot_3D(filename1,filename2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                 %%%
%%%    Program to make 3D plots                                     %%%
%%%                                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Read the netCDF file
nc=netcdf(filename1);
%ncdump(nc);
mlon=nc{'longitude',1}(:,:);
mlat=nc{'latitude',1}(:,:);
mdepth=nc{'depth',1}(:,:);
mdepth2=mdepth;


nx=size(mdepth,1);
ny=size(mdepth,2);

mlon   = mlon   (4:nx-3,4:ny-3);
mlat   = mlat   (4:nx-3,4:ny-3);
mdepth = mdepth (4:nx-3,4:ny-3);
mdepth2= mdepth2(4:nx-3,4:ny-3);


max(max(mdepth))
min(min(mdepth))
I=find(mdepth<0);
mdepth (I)=0.;
mdepth2(I)=0.;
I=find(mdepth>0.1);
mdepth2(I)=nan;

S=surf(mlon,mlat,-mdepth);
set(S,'FaceColor',[.5 .5 .5]);
set(S,'FaceLighting','gouraud');
set(S,'EdgeColor','none');
%set(S,'EdgeColor',get(S,'FaceColor'));
%set(S,'EdgeLighting',get(S,'FaceLighting'));


hold on;
S2=surf(mlon,mlat,-mdepth2+.1);
set(S2,'FaceColor',[.5 .8 .5]);
set(S2,'EdgeColor',get(S2,'FaceColor'));
set(S2,'FaceLighting','none');
%material metal;
%lighting gouraud;
%light;
%lightangle(210,30);
%view(210,30);

nc2=netcdf(filename2);
%ncdump(nc2);

slon=nc2{'longitude',1}(:);
slat=nc2{'latitude',1}(:);
sdepth=nc2{'depthc',1}(1,:,:);
stemp=nc2{'TEM',1}(1,:,:);

slon2=repmat(slon,1,size(sdepth,1));
slat2=repmat(slat,1,size(sdepth,1));

I=find(stemp<-1000);
stemp(I)=NaN;


S=surf(slon2',slat2',sdepth,stemp);
set(S,'FaceColor')
set(S,'FaceLighting')
set(S,'FaceColor','interp')
set(S,'FaceLighting','none')
set(S,'EdgeColor','none')
axis off;

material metal;
lighting gouraud;
light;
lightangle(210,30);
view(210,30);
caxis([5 8])


