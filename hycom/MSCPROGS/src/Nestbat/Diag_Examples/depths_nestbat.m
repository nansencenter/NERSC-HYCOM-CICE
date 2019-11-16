
% Load nestbat diag file 
load nestbat.mat;

% Original depths
figure(1);
I=find(ldepths<.1); ldepths(I)=nan;    
P=pcolor(ldepths'); set(P,'LineStyle','none') ; caxis([0 4000])     
title('Original depth')

% From global
figure(2);
I=find(depths_global_to_local<.1); depths_global_to_local(I)=nan;    
P=pcolor(depths_global_to_local'); set(P,'LineStyle','none') ; caxis([0 4000])     
title('Global depths on local grid')

% Final depth
figure(3)
I=find(depths_final<.1); depths_final(I)=nan;    
P=pcolor(depths_final'); set(P,'LineStyle','none') ; caxis([0 4000])     
title('Final depths')

% Final - Original
figure(4)
I=find(isnan(ldepths)); ldepths(I)=0.;
I=find(isnan(depths_final)); depths_final(I)=0.;
tmp=depths_final-ldepths;
I=find(ldepths==0. & depths_final==0.);
tmp(I)=nan;
P=pcolor(tmp'); set(P,'LineStyle','none') ; colorbar;
title('Final depths - Original depths')

