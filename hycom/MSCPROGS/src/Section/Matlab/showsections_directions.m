function fig=showsections_directions(lonrange,latrange)
% A function which plots time record <record_number>
% from section <section_number> 
if (nargin==0)
   lonrng=[-180 180];
   latrng=[-85 85];
elseif (nargin==1)
   lonrng=lonrange;
   latrng=[-85 85];
elseif (nargin==2)
   lonrng=lonrange;
   latrng=latrange;
else
   disp('Max 2 arguments');
   return
end




% Check for section file
nsec=0;
for sec=1:100
   csec=num2str(sec,'%3.3d');
   fname=['section' csec '.dat'];
   fid=fopen(fname);
   if (fid>0)
      nsec=sec;
      % First line is header
      line=fgetl(fid);
      line=fgetl(fid);
      count=1;
      [data, ncols, errmsg, nxtindex] = sscanf(line, '%f');
      while (~isempty(data) & isstr(line))

         line=fgetl(fid);
         if (isstr(line))
            [data, ncols, errmsg, nxtindex] = sscanf(line, '%f');
            ipiv(sec,count) = data(1);
            jpiv(sec,count) = data(2);
            lon(sec,count) = data(3);
            lat(sec,count) = data(4);
            xang(sec,count)=data(6);
            yang(sec,count)=data(7);
            ncount(sec)=count;
            count=count+1;
            %disp([ num2str(sec) ' ' num2str(ncount(sec))]);
         end 
      end 
      fclose(fid);
   end 
end
for sec=1:nsec
   for is=1:ncount(sec)
      test1(ipiv(sec,is),jpiv(sec,is)) = xang(sec,is)/pi;
      test2(ipiv(sec,is),jpiv(sec,is)) = yang(sec,is)/pi;

      x1(ipiv(sec,is),jpiv(sec,is)) = cos(xang(sec,is));
      y1(ipiv(sec,is),jpiv(sec,is)) = sin(xang(sec,is));
      x2(ipiv(sec,is),jpiv(sec,is)) = cos(yang(sec,is));
      y2(ipiv(sec,is),jpiv(sec,is)) = sin(yang(sec,is));

      iind(ipiv(sec,is),jpiv(sec,is)) = ipiv(sec,is);
      jind(ipiv(sec,is),jpiv(sec,is)) = jpiv(sec,is);

   end 
end
I=find(abs(test1)<0.001); test1(I)=nan; 
I=find(abs(test2)<0.001); test2(I)=nan;


figure(1); clf; P=pcolor(test1');
set(P,'LineStyle','none')
caxis([-2 2]);
title(' X - angle');
xlabel(' X - direction');
ylabel(' Y - direction');
C=colorbar;
hold on;
for sec=1:nsec
   text(ipiv(sec,1),jpiv(sec,1),[num2str(sec) ' start' ]);
end 
I=find(x1.^2+y1.^2 > 0.001);
%x1(I)=nan;
%y1(I)=nan;
quiver(iind(I),jind(I),x1(I),y1(I));
hold off;
axes(C);
xlabel('ang [\pi]');

figure(2);clf; P=pcolor(test2');
set(P,'LineStyle','none')
caxis([-2 2]);
title(' Y - angle');
xlabel(' X - direction');
ylabel(' Y - direction');
C=colorbar;
hold on;
for sec=1:nsec
   text(ipiv(sec,1),jpiv(sec,1),[num2str(sec) ' start' ]);
end 
I=find(x2.^2+y2.^2 > 0.001);
%x1(I)=nan;
%y1(I)=nan;
quiver(iind(I),jind(I),x2(I),y2(I));
hold off;
axes(C);
xlabel('ang [\pi]');




%figure(3);clf;
%I=find(x1.^2+y1.^2 > 0.001);
%%x1(I)=nan;
%%y1(I)=nan;
%quiver(iind(I),jind(I),x1(I),y1(I));



%% Initialize projection
%gcf; clf;
%m_proj('Miller','lon',lonrng,'lat',latrng);
%hold on;
%
%%Mask
%m_coast('patch','k');
%m_grid
%
%% 
%for sec=1:nsec
%   % Lines
%   m_track(lon(sec,1:ncount(sec)),lat(sec,1:ncount(sec)),'color','r', ...
%           'linew',2);
%
%   % For text
%   [x,y]=m_ll2xy(lon(sec,1:ncount(sec)),lat(sec,1:ncount(sec)));
%
%
%   % Trial and error, the rsesulting box around the number
%   % should be ok with "zoomfac" logic below
%   % -- There probably is an easier way to do this
%   XLim=get(gca,'XLim');
%   YLim=get(gca,'YLim');
%   zoomfac=max((XLim(2)-XLim(1))/6,(YLim(2)-YLim(1))*1.4/6);
%
%
%   % Backdrop for text 
%   P=patch(x(1)+.08*zoomfac*[-1 1 1 -1],y(1)+.08*zoomfac*[-1 -1 1 1],[.8 .8 .8]);
%
%   % Text -- centered on patch above
%   T=text(x(1),y(1),num2str(sec));
%   set(T,'FontWeight','bold');
%   set(T,'FontSize',20);
%   %set(T,'FontUnits','normalized');
%   set(T,'FontUnits','points');
%   set(T,'HorizontalAlignment','center');
%   set(T,'VerticalAlignment','middle');
%
%end
%%set(T,'FontUnits')
%%set(P)
%
if (nsec==0) 
   disp('No sections found...');
end
