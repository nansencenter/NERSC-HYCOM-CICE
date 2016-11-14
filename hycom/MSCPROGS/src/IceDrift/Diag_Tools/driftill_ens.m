% illustration of ice drift from model vs obs (CERSAT)
figure(1); clf;
set(gca,'FontSize',16); set(gca,'FontWeight','bold');
orient tall
hold on

m_proj('stereographic','lon',15,'lat',80,'rad',10,'rec','on')

colors=['r' 'b' 'g' 'k' 'c' 'y' 'm']
Qset=[];

for iens=1:3

   cmem=num2str(iens,'%3.3i')

   cindex=mod(iens-1,7)+1;
   col=colors(cindex:cindex)

   load([ 'drift' cmem '.mat'])

   K=find(o_zonal==max(max(o_zonal)));
   o_zonal(K)=0;
   K=find(o_merid==max(max(o_merid)));
   o_merid(K)=0;





   if (iens==1) 
      %Subsets
      step=4;
      I=[];
      for i=1:step:size(o_longitude,1);
      for j=1:step:size(o_longitude,2);
         %disp([num2str(i) ' ' num2str(j) ' ' num2str(o_longitude(i,j)) ]);
         Itmp=sub2ind(size(o_longitude),i,j);
         I=[Itmp I];
      end
      end
   end




   slon=o_longitude(I);
   slat=o_latitude(I);
   somer=o_merid(I);
   sozon=o_zonal(I);
   smmer=m_merid(I);
   smzon=m_zonal(I);

   sodist=sqrt(somer.^2+sozon.^2);
   smdist=sqrt(smmer.^2+smzon.^2);


   sfac=5e-5;
   K=find(smdist>4000);
   Q=m_quiver(slon(K),slat(K),smzon(K)*sfac,smmer(K)*sfac,0,'Color',col,'LineWidth',1);

   Qset=[Qset Q];
   %m_vec(200000,slon(K),slat(K),smzon(K),smmer(K), ...
   %'shaftwidth',.2 , 'headwidth', 1, ...
   %'FaceColor',col,'EdgeColor',col)


end



m_gshhs_i('patch','k')
m_grid

legend(Qset, 'mem. 1','mem. 2','mem. 3')

%legend('Member 1','Obs')

%title('Ice Drift 20070204-20070207 -- Red: Obs, Blue:Ensemble Member 1')
print('-depsc2','drft-20070204-20070207-mem1-3.eps')
print('-dpng','drft-20070204-20070207-mem1-3.png')


