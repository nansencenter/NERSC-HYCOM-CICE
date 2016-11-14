function cmp_nest(fname1,fname2,indlw,indup,cvar)
%%% Simple function to compare nesting conditions of a given boundary file
%function cmp_nest(fname1,fname2,indlw,indup,cvar)
%function show_nest(fname)
%
% fname = name of nesting file
%   



% Sample diag plot - fixed time, varying boundary

m_proj('stereographic','lon',30,'lat',75,'rad',25,'rec','on')

for i2=indlw:indup

   cind=num2str(i2,'%4.4i');

   load(fname1);
   eval(['fld1 = ' cvar cind ';']);
   eval(['int1 = interface' cind ';']);
   eval(['lon1 = longitude' cind ';']);
   eval(['lat1 = latitude' cind ';']);
   eval(['dst1 = grid_distance;']);


   load(fname2);
   eval(['fld2 = ' cvar cind ';']);
   eval(['int2 = interface' cind ';']);
   eval(['lon2 = longitude' cind ';']);
   eval(['lat2 = latitude' cind ';']);
   eval(['dst2 = grid_distance;']);



   figure(1);
   clf;
   P=pcolor(dst1,-int1/9806,fld1); set(P,'LineStyle','none')
   %caxis(cax);
   hold on;
   for j=1:size(dst1,2)
      plot(dst1(:,j),-int1(:,j)/9806,'k')
   end
   colorbar;

   figure(2);
   clf;
   P=pcolor(dst2,-int2/9806,fld2); set(P,'LineStyle','none')
   %caxis(cax);
   hold on;
   for j=1:size(dst2,2)
      plot(dst2(:,j),-int2(:,j)/9806,'k')
   end
   colorbar;

   figure(4);
   clf;
   hold on;
   for j=1:size(dst2,2)
      plot(dst2(:,j),-int2(:,j)/9806,'k')
   end
   for j=1:size(dst1,2)
      plot(dst1(:,j),-int1(:,j)/9806,'k--')
   end

   if (i2==indlw)
      figure(3); clf;
      m_gshhs_i('patch','k')
      m_grid;
      hold on;
   end
   figure(3); 
   m_plot(lon1(:,j),lat1(:,j),'LineWidth',2,'Color','r');

   %mind=min(min(dst2))
   %maxd=max(max(dst2))
   %mini=min(min(int2))
   %maxi=max(max(int2))
   %[XI,YI]=meshgrid(mind:(maxd-mind)/100:maxd, ...
   %                 mini:(maxi-mini)/100:maxi);
   %ZI1=interp2(dst1,int1,fld1,XI,YI);
   %figure(5) ; clf; pcolor(XI,-YI/9806,ZI1); shading flat;

   disp('press a key')
   pause
end
