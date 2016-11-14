function show_nest(fname,indlw,indup,cvar,loncenter,latcenter,radius)
%%% Simple function to plot nesting conditions of a given boundary file
%function show_nest(fname,indlw,indup,cvar,loncenter,latcenter,radius)
%
% fname = name of nesting file
% indlw = lower index
% indup = lower index
% cvar  = name of variable to plot
% loncenter  = longitude center of projection
% latcenter  = latitude  center of projection
% radius     = "radius" of projection in degrees
%
% fname is the name of a file produced by checknest.sh
%
% indlw, indup are indexes which depend on the boundary - for "i1" and "j1", 
% it ranges from 1 to 20. For ii and jj it depends on grid size. If idm is 400, 
% then indlw and indup must be in the range 381 to 400. More generally, for ii,
% indlw, indup must be in the range idm-20,idm. For jj it must be in the range 
% jdm-20, jdm.
%
% The boundaries are illustrated as lines on the model grid - on a stereographic projection. 
% loncenter, latcenter and radious are inputs to this projection using m_map (see m_map
% documentation)



% Sample diag plot - fixed time, varying boundary
files=dir(fname);
m_proj('stereographic','lon',loncenter,'lat',latcenter,'rad',radius,'rec','on')
for i=1:prod(size(files))

   disp(files(i).name);
   load(files(i).name);
   Plist=[];
   for i2=indlw:indup

      cind=num2str(i2,'%4.4i');
      disp(cind);

      fldtype='section';
      if (strcmp(cvar,'saln'))
         eval(['fld = saln' cind ';']);
         cax=[33 35];
      elseif (strcmp(cvar,'temp'))
         eval(['fld = temp' cind ';']);
         cax=[-2 6];
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
      elseif (strcmp(cvar,'fice'))
         eval(['fld = fice' cind ';']);
         cax=[0 1];
         fldtype='track';
      elseif (strcmp(cvar,'hice'))
         eval(['fld = hice' cind ';']);
         cax=[0 1];
         fldtype='track';
      else
         disp('unknown field')
      end
         
      eval(['interface = interface' cind ';']);
      eval(['longitude = longitude' cind ';']);
      eval(['latitude = latitude' cind ';']);

      if (strcmp(fldtype,'section'))
         figure(1);
         hold off;
         P=pcolor(grid_distance,-interface/9806,fld); set(P,'LineStyle','none')
         caxis(cax);
         hold on;
         for j=1:size(grid_distance,2)
            %j
            plot(grid_distance(:,j),-interface(:,j)/9806,'k')
         end
         colorbar;
      elseif (strcmp(fldtype,'track'))
         figure(1);
         if (i2==indlw)
            clf;
            hold on;
         else
            for q=1:prod(size(Plist))
               q
               set(Plist(q),'Color',[.8 .8 .8]);
               set(Plist(q),'LineWidth',1);
            end
         end 
         P=plot(grid_distance(:,1),fld,'r','LineWidth',2)
         Plist=[Plist P]

      end 


      if (i2==indlw)
         figure(2); clf;
         m_coast('patch','k')
         m_grid;
         hold on;
      end
      figure(2); 
      m_plot(longitude(:,1),latitude(:,1),'LineWidth',2,'Color','r');

      disp('press a key')
      pause

   end
end
