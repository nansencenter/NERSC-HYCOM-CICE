% Load nestbat diag file 
load nestbat.mat;

m_proj('stereographic','lon',40,'lat',75,'rad',25,'rec','on')

[lx,ly]=m_ll2xy(llon,llat);
[gx,gy]=m_ll2xy(glon,glat);


figure(1)
P=pcolor(gx,gy,gdepths) ; caxis( [0 4000] ); 
%shading interp; % -- May be useful
set(P,'EdgeAlpha',.2) % Transparent grid lines
hold on
nxl=size(lx,1);
nyl=size(lx,2);
P=pcolor(lx(2:nxl,2:nyl),ly(2:nxl,2:nyl),depths_final(2:nxl,2:nyl)) ; caxis( [0 4000] );
set(P,'EdgeAlpha',.2) % Transparent grid lines
m_coast('patch','k')
m_grid
title('Final local grid on top of global')
hold off


figure(2)
P=pcolor(gx,gy,gdepths) ; caxis( [0 4000] ); 
%shading interp; % -- May be useful
set(P,'EdgeAlpha',.2) % Transparent grid lines
hold on
nxl=size(lx,1);
nyl=size(lx,2);
P=pcolor(lx(2:nxl,2:nyl),ly(2:nxl,2:nyl),ldepths(2:nxl,2:nyl)) ; caxis( [0 4000] );
set(P,'EdgeAlpha',.2) % Transparent grid lines
m_coast('patch','k')
m_grid
title('Initial local grid on top of global')
hold off

