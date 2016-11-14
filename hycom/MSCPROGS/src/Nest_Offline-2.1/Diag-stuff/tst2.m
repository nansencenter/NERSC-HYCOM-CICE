% Sample diag plot - fixed time, varying boundary
cbnd='i1' % Boundary
files=dir(['nest_9999_*' cbnd '_test.mat']);
indlw=1;
indup=20;
m_proj('stereographic','lon',30,'lat',75,'rad',25,'rec','on')
cvar='vtot'
for i=1:prod(size(files))

   disp(files(i).name);
   load(files(i).name);
   for i2=indlw:indup

      cind=num2str(i2,'%4.4i');
      disp(cind);

      if (strcmp(cvar,'saln'))
         eval(['fld = saln' cind ';']);
         cax=[33 35];
      elseif (strcmp(cvar,'temp'))
         eval(['fld = temp' cind ';']);
         cax=[-2 10];
      elseif (strcmp(cvar,'utot'))
         eval(['fld = utot' cind ';']);
         cax=[-.2 .2];
      elseif (strcmp(cvar,'vtot'))
         eval(['fld = vtot' cind ';']);
         cax=[-.2 .2];
      elseif (strcmp(cvar,'ubaclin'))
         eval(['fld = ubaclin' cind ';']);
         cax=[-.2 .2];
      elseif (strcmp(cvar,'vbaclin'))
         eval(['fld = vbaclin' cind ';']);
         cax=[-.2 .2];
      else
         disp('unknown field')
      end
         
      eval(['interface = interface' cind ';']);
      eval(['longitude = longitude' cind ';']);
      eval(['latitude = latitude' cind ';']);

      figure(1);
      P=pcolor(grid_distance,-interface/9806,fld); set(P,'LineStyle','none')
      caxis(cax);
      colorbar;

         % Adjust  ...
      if (i2==indlw)
         figure(2); clf
      else
         figure(2); hold on; 
      end

      m_plot(longitude(:,1),latitude(:,1),'LineWidth',2,'Color','r');

      if (i2==indup)
         m_coast('patch','k')
         m_grid;
      end

      disp('press a key')
      pause

   end
end
