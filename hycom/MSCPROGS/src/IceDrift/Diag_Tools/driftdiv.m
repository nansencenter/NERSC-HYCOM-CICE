% illustration of ice drift from model vs obs (CERSAT)

function driftdiv(imem);


cmem=num2str(imem,'%3.3i');
load([ 'drift' cmem '.mat' ])

K=find(o_zonal==max(max(o_zonal)));
o_zonal(K)=0;
K=find(o_merid==max(max(o_merid)));
o_merid(K)=0;

%Subsets
step=4;
I=[];
for i=1:step:size(o_longitude,1);
for j=1:step:size(o_longitude,2);

   disp([num2str(i) ' ' num2str(j) ' ' num2str(o_longitude(i,j)) ]);
   Itmp=sub2ind(size(o_longitude),i,j);
   I=[Itmp I];
end
end




slon=o_longitude(I);
slat=o_latitude(I);
somer=o_merid(I);
sozon=o_zonal(I);
smmer=m_merid(I);
smzon=m_zonal(I);

%slon=longitude(I);
%slat=latitude(I);
%somer=ones(1,prod(size(I)))*1000;
%sozon=ones(1,prod(size(I)))*1000;
%smmer=ones(1,prod(size(I)))*1000;
%smzon=ones(1,prod(size(I)))*1000;


sodist=sqrt(somer.^2+sozon.^2);
smdist=sqrt(smmer.^2+smzon.^2);

%m_quiver(longitude,latitude,m_zonal,m_merid,10)
%m_vec(0.02,longitude(:,300),latitude(:,300),m_zonal(:,300),m_merid(:,300))
figure(1); clf; set(gca,'FontSize',16); set(gca,'FontWeight','bold');
orient tall
m_proj('stereographic','lon',-45,'lat',90,'rad',30,'rec','on')

K=find(smdist>4000);
%Q=m_quiver(slon,slat,smzon*sfac,smmer*sfac,0,'Color','r','LineWidth',2);
m_vec(200000,slon(K),slat(K),smzon(K),smmer(K), ...
'shaftwidth',.6 , 'headwidth', 3, ...
'FaceColor','b','EdgeColor','b')


hold on
K=find(sodist>4000);
%Q=m_quiver(slon(K),slat(K),sozon(K)*sfac,somer(K)*sfac,0,'Color','b','LineWidth',2);
m_vec(200000,slon(K),slat(K),sozon(K),somer(K), ...
'shaftwidth',.6 , 'headwidth', 3, ...
'FaceColor','r','EdgeColor','r')


% Reference arrow
ang=-pi*45/180;
r=100000
[rv,rt] =m_vec(200000,0,68,r*cos(ang),r*sin(ang), ...
'shaftwidth',.6 , 'headwidth', 3, 'key', ' 100 km / 3 d', ...
'FaceColor','k','EdgeColor','k')

set(rt,'FontWeight','bold')


m_coast('patch','k')
m_grid

%legend('Member 1','Obs')

title(['Ice Drift -- Red: Obs, Blue:Ensemble Member ' cmem ])

print('-depsc2','drft.eps')
print('-dpng','drft.png')
