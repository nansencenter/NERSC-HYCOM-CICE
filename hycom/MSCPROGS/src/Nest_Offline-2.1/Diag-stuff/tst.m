% Sample diag plot - fixed boundary, varying time
cbnd='jj'
files=dir(['nest_9999_*' cbnd '_test.mat']);
ind=206;

cind=num2str(ind,'%4.4i')
for i=1:prod(size(files))


   disp(files(i).name);
   load(files(i).name);
   eval(['saln = saln' cind ';']);
   eval(['interface = interface' cind ';']);
   eval(['longitude = longitude' cind ';']);
   eval(['latitude = latitude' cind ';']);

   figure(1);
   P=pcolor(grid_distance,-interface/9806,saln); set(P,'LineStyle','none')
   caxis([33 35]);
   colorbar;

   if (i==1) 
      % Adjust  ...
      m_proj('stereographic','lon',30,'lat',75,'rad',25,'rec','on')
      figure(2);
      m_plot(longitude(:,1),latitude(:,1),'LineWidth',2,'Color','r');
      m_coast('patch','k')
      m_grid;
   end

   disp('press a key')
   pause

end
