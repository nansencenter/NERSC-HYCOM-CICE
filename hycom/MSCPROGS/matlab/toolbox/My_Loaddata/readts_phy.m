function[]=readts_phy(var,region,year);
if(nargin~=3)
   disp('specify var(temp,saln,nitr etc.), region and year')
   return
end
%close all;
syear=year;
eyear=year;
start_mon=1;
end_mon=11;
smont=syear*12+start_mon;
emont=eyear*12+end_mon;

months=emont-smont+1;
if (smont/12 < 1000)
   steps=[719 719 719 719 719 719 ...
          719 719 719 719 719 719];
elseif((year/4)-floor(year/4)==0) %skuddår
   steps=[743 695 743 719 743 719 ...
          743 743 719 743 719 743];
else
   steps=[743 671 743 719 743 719 ...
          743 743 719 743 719 743];
end
msteps=max(steps);

NOL=23;

tot_steps=sum(steps(start_mon:end_mon))

z=zeros(NOL,tot_steps);
data=zeros(NOL,tot_steps);
time=zeros(1,tot_steps);
SP=1;EP=steps(mod(smont-1,12)+1);
for mon=smont:emont
   mon
step=steps(mod(mon-1,12)+1);
     %[A , B]=readts_mon_phy(mon,var,region,step);
   [z(:,SP:EP),data(:,SP:EP),time(SP:EP)]=readts_mon_phy(mon,var,region,step);
     SP=EP+1;
     EP=EP+steps(mod(mon,12)+1);
end

%step=12;
%R=1:step:tot_steps;

myear=floor(time/10e5);
month=floor((time-myear*10^6)/10e3);
day=floor((time-myear*10^6-month*10^4)/10^2);
hour=floor(time-myear*10^6-month*10^4-day*10^2);

% --- finding daily values

intday=datenum(myear,month,day);
intday=intday+1; % From hycom date to normal date
min_day=min(intday)
max_day=max(intday)
nr_days=max_day-min_day;

h=0;
for i=min_day:max_day
   h=h+1;
   I=find(intday==i);
   for k=1:NOL
      test(k,h)=sum(data(k,I))/length(I);
      depth(k,h)=sum(z(k,I))/length(I);
   end
end

R=min_day:1:max_day;

[x,y]=meshgrid(R,1:NOL); 
start_day=datenum(year,1,1);
x_new=x-start_day;
%pcolor(x./24,-depth,test);shading flat;hold on;
pcolor(x_new,-depth,test);shading flat;hold on;
ylim([-600 0]);

%plot(R/24,data(1,R)-data(5,R))
%period=datenum(syear:eyear,01,01);
%if (strcmp(var,'SEDI')|strcmp(var,'YELL')) 
%  pcolor(x./24,-squeeze(z(1:NOL,R)),log(data(1:NOL,R)));shading flat;hold on;
%  ylabel(var)
% % log_colorbar;
%else
%  %pcolor(x./24,-squeeze(z(1:NOL,R)),data(1:NOL,R)); shading flat;
%  pcolor(x,-depth(1:NOL,R),test(1:NOL,R)); shading flat;
%  if(strcmp(var,'saln'))
%     caxis([35.05 35.25]);
%     axis([min_day max_day -400 0])
%  elseif(strcmp(var,'temp'))
%     caxis([4 12]);
%     axis([min_day max_day -200 0])
%  elseif(strcmp(var,'nit'))
%     caxis([10 170]);
%     axis([min_day max_day -200 0])
%  elseif(strcmp(var,'dia'))   
%     caxis([0 50]);
%     axis([min_day max_day -120 0])
%  elseif(strcmp(var,'fla'))   
%     caxis([0 40]);
%     axis([min_day max_day -120 0])
%  elseif(strcmp(var,'pho'))   
%     caxis([0 35]);
%     axis([min_day max_day -120 0])
%  elseif(strcmp(var,'sil'))   
%     caxis([0 250]);
%     axis([min_day max_day -120 0])
%  end
%  title([var ' in NW4, year ' num2str(year) ' in gp ', num2str(region)])
%  ylabel('depth') 
%%  set(gca,'XTick',period,'XTicklabel',datestr(period,'yyyy'))
%%  axis([datenum(syear,01,01) datenum(eyear,01,01+31) -1200 0])
%  datetick('x',3);
%  xlabel([num2str(syear) '-' num2str(eyear)])
%  colorbar;
%  hold off;
%end 
%%print('-djpeg',['Figs_1995/',model,'_ts_',var,'_gp_',num2str(region),'.jpg']);
%%print('-dpng',['Figs_1995/',model,'_ts_',var,'_gp_',num2str(region),'.png']);
