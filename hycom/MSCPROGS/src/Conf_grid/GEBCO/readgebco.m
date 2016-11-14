

% The GEBCO dataset is a perfect example of how NOT to use netcdf....
nc=netcdf('gridone.grd');
dims=nc{'dimension',1}(:)
xrange=nc{'x_range',1}(:)
yrange=nc{'y_range',1}(:)
spacing=nc{'spacing',1}(:)



minlon=-30;
maxlon=30;
lonskip=1;
minlat=  40;
maxlat=  65;
latskip=1;

minloni=min(dims(1),max(round(minlon-xrange(1))/spacing(1),1))
maxloni=min(dims(1),max(round(maxlon-xrange(1))/spacing(1),1))
minlati=min(dims(2),max(round(minlat-yrange(1))/spacing(2),1))
maxlati=min(dims(2),max(round(maxlat-yrange(1))/spacing(2),1))

tmp=maxlati
%minlati=1
%maxlati=5400
maxlati=dims(2)-minlati
minlati=dims(2)-tmp


% Read the data, one row at a time
fld=[];
%for irow=minloni:lonskip:maxloni
%   tmp1 = nc{'z',1}( dims(1)*(irow-1)+minlati:latskip:dims(1)*(irow-1)+maxlati )  ;
%   fld=[fld  tmp1];
%   %disp([ num2str(irow) ' ' num2str(maxloni) ]);
%end

for irow=minlati:latskip:maxlati
   tmp1 = nc{'z',1}( dims(1)*(irow-1)+minloni:lonskip:dims(1)*(irow-1)+maxloni )  ;
   fld=[fld  tmp1];
   disp([ num2str(irow) ' ' num2str(maxlati) ]);
end
lats=max(yrange)-[minlati:latskip:maxlati]*spacing(2);
lons=[minloni:lonskip:maxloni]*spacing(1)+min(xrange);


%skip=10;
%xdm=floor((dims(1)-1)/skip)+1;
%ydm=floor((dims(2)-1)/skip)+1;

%test=zeros(xdm,ydm);
%for j = 1:ydm
%   j2=j*skip;
%   tmp1 = nc{'z',1}( dims(1)*(j2-1)+1:dims(1)*j2 )  ;
%   test(:,j) = tmp1(1:skip:dims(1));
%end

%test=min(test,0);
%fld=min(fld,200);
%fld=max(fld,-1200);
%contourf(-fld)
pcolor(rot90(fld,1)) ; shading flat
colorbar


%test=zeros(ydm,xdm);
%for j = 1:xdm
%   j2=j*skip;
%   tmp1 = nc{'z',1}( dims(2)*(j2-1)+1:dims(2)*j2 )  ;
%   test(:,j) = tmp1(1:skip:dims(2));
%end

