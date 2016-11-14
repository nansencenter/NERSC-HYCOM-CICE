% NetCDF file contains prior/after fields for assimilation
nc=netcdf('tst.nc') ; 
m_proj('stereographic','lon',-45,'lat',90,'rad',30,'rec','on');

lon=nc{'longitude'}(:,:);
lat=nc{'latitude'}(:,:);

temp=nc{'temp01'}(:,:);
saln=nc{'saln01'}(:,:);
fice=nc{'ficem00'}(:,:);
Hice=nc{'hicem00'}(:,:);

fillval=nc{'temp01'}.FillValue_(:);


[x,y]=m_ll2xy(lon,lat);
I=find(temp==fillval); temp(I)=nan;
I=find(saln==fillval); saln(I)=nan;
I=find(fice==fillval); fice(I)=nan;
I=find(hice==fillval); hice(I)=nan;

figure(1); clf
P=pcolor(x,y,reshape(fice(1,:,:),size(x)));
set(gca,'FontSize',16) ; set(gca,'FontWeight','bold');
set(P,'LineStyle','none')
m_coast('patch','k'); m_grid;
title('c_i 7th Feb 2007')
print('-depsc2','cice20070207.eps')
print('-dpng','cice20070207.png')

figure(2); clf
P=pcolor(x,y,reshape(fice(2,:,:)-fice(1,:,:),size(x)));
set(gca,'FontSize',16) ; set(gca,'FontWeight','bold');
set(P,'LineStyle','none')
caxis([-.5 .5]);
cmap=colormap; cmap(30:35,:)=1; colormap(cmap);
colorbar
m_coast('patch','k'); m_grid;
title(' c_i Update 7th Feb 2007')
print('-depsc2','deltacice20070207.eps')
print('-dpng','deltacice20070207.png')

figure(3); clf
tmp=reshape(hice(1,:,:),size(x));
I=find(tmp<.05); tmp(I)=nan;
P=pcolor(x,y,tmp);
set(gca,'FontSize',16) ; set(gca,'FontWeight','bold');
set(P,'LineStyle','none')
caxis([0 4]);
colorbar
m_coast('patch','k'); m_grid;
title('h_i 7th Feb 2007')
print('-depsc2','hice20070207.eps')
print('-dpng','hice20070207.png')

figure(4); clf
P=pcolor(x,y,reshape(hice(2,:,:)-hice(1,:,:),size(x)));
set(gca,'FontSize',16) ; set(gca,'FontWeight','bold');
set(P,'LineStyle','none')
caxis([-2 2]);
cmap=colormap; cmap(30:35,:)=1; colormap(cmap);
colorbar
m_coast('patch','k'); m_grid;
title('h_i Update 7th Feb 2007')
print('-depsc2','deltahice20070207.eps')
print('-dpng','deltahice20070207.png')



figure(5); clf
P=pcolor(x,y,reshape(temp(1,:,:),size(x)));
set(gca,'FontSize',16) ; set(gca,'FontWeight','bold');
set(P,'LineStyle','none')
m_coast('patch','k'); m_grid;
title('sst 7th Feb 2007')
print('-depsc2','sst20070207.eps')
print('-dpng','sst20070207.png')

figure(6); clf
P=pcolor(x,y,reshape(temp(2,:,:)-temp(1,:,:),size(x)));
set(gca,'FontSize',16) ; set(gca,'FontWeight','bold');
set(P,'LineStyle','none')
caxis([-1 1]);
cmap=colormap; cmap(30:35,:)=1; colormap(cmap);
colorbar
m_coast('patch','k'); m_grid;
title(' sst Update 7th Feb 2007')
print('-depsc2','deltasst20070207.eps')
print('-dpng','deltasst20070207.png')



