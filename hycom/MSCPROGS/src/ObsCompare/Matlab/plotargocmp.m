function plotargocmp(indir)
%List files
files=dir([ indir 'profile_date*'])

for i=1:prod(size(files))
   fname=files(i).name;
   dat=load([ indir fname] );


   % Get date from profile
   tmp=fname;
   tmp=strrep(tmp,'profile_date','');
   idate=regexprep(tmp,'_.*','');
   pos=regexprep(tmp,'.*pos','');
   epos=regexprep(pos,'Ex.*','');
   npos=regexprep(pos,'.*Ex',''); npos=regexprep(npos,'N.*','');
   disp([ fname ' ' epos ' ' npos]);

   %plot salinity
   I=find(dat(:,2)~=-999);
   J=find(dat(:,3)~=-999);
   if (prod(size(I))>1 & prod(size(J))>1)
      figure(1); clf ; hold on ;
      H1=[]; H2=[];
      legtxt1=[]; legtxt2=[];
      if (prod(size(I))>0)
         H1=plot(dat(I,2),-dat(I,1),'LineWidth',2);
         legtxt1='Data ';
      end
      if (prod(size(J))>0)
         H2=plot(dat(J,3),-dat(J,1),'LineWidth',2,'Color','r');
         legtxt2='Model';
      end
      set(gca,'FontSize',14);
      set(gca,'FontWeight','bold');
      grid on;
      H=[H1 H2];
      legtxt=[legtxt1 ; legtxt2];
      legend(H,legtxt); %,'Location','Best');
      xlabel('Salinity[psu]');
      ylabel('Depth[m]');
      print('-dpng','-r150',[ 'saltprofile_' pos '.png'])
      %print('-depsc2','-r300',[ 'saltprofile_' pos '.eps'])
      close(1);
   end

   %plot temperature
   I=find(dat(:,4)~=-999);
   J=find(dat(:,5)~=-999);
   H1=[]; H2=[];
   legtxt1=[]; legtxt2=[];
   if (prod(size(I))>1 & prod(size(J))>1)
      figure(2); clf ; hold on ;
      if (prod(size(I))>0)
         H1=plot(dat(I,4),-dat(I,1),'LineWidth',2);
         legtxt1='Data ';
      end
      if (prod(size(J))>0)
         H2=plot(dat(J,5),-dat(J,1),'LineWidth',2,'Color','r');
         legtxt2='Model';
      end
      set(gca,'FontSize',14);
      set(gca,'FontWeight','bold');
      grid on;
      H=[H1 H2];
      legtxt=[legtxt1 ; legtxt2];
      legend(H,legtxt); %,'Location','Best');
      xlabel('Temperature[C]');
      ylabel('Depth[m]');
      print('-dpng','-r150',[ 'tempprofile_' pos '.png'])
      close(2);
   end
end
