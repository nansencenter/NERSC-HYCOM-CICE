function trip_riverbycatch(irec);

nc=netcdf('test.nc');

%lon, lat

lon=nc{'riverbycatch_lon'}(irec);
lat=nc{'riverbycatch_lat'}(irec);
lon


figure(1);clf ; hold on
m_proj('stereographic','lon',lon,'lat',lat,'rad',35,'rec','on');
m_gshhs_i('patch',[.3 .3 .3]); %,'EdgeColor',[0.3 0.3 0.3])
m_grid();
hold on;
M=m_plot(lon,lat,'r+');
set(M,'MarkerSize',40)
set(M,'MarkerFaceColor','r')
set(M,'LineWidth',2)
get(M)
pos=get(gca,'Position');
pos(4)=pos(4)-0.2;
pos(2)=pos(2)+0.2;
set(gca,'Position',pos);


pause(0.5)

tseries=nc{'riverbycatch'}(:,irec);
t=nc{'record'}(:,irec);
axes('Position',[0.13 0.05 0.775 0.18])
plot(t,tseries*1e-3,'LineWidth',2);
grid on
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
title('Discharge [1000 m^3 s^-1]')

%set(gca,'XColor','r')
%set(gca,'YColor','r')
