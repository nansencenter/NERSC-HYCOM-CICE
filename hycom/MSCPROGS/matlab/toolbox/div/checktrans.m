function checktrans(isec,idm,jdm);
%function checktrans(isec,idm,jdm);
% Function checks transport setup from the "section_intersect" routine
% (This routine is called whenever you run pak2sec, pak2trans and pak2trans2)
%
%Remember to run that routine before running this matlab routine...
%
%input:
% isec - number of section to check
% idm  - 1st model grid dimension
% jdm  - 2nd model grid dimension
%
% The arrows indicate positive direction of transport
% Grey squares indicate gridcells involved in transport calculation

css=num2str(isec,'%3.3i');

fdepths=['depths' num2str(idm,'%3.3i') 'x' num2str(jdm,'%3.3i') '.uf'];
fid=fopen(['depths' num2str(idm,'%3.3i') 'x' num2str(jdm,'%3.3i') '.uf'], ...
          'r','ieee-be');
if (fid~=-1)
   stat=fseek(fid,4,'bof'); % Skip fortran 4-byte header 
   depths=fread(fid,[idm jdm],'double');
   fclose(fid);
else
   disp([fdepths ' not found ']);
   return
end
depths=depths';

% 1) Load transportxxx.dat
fid=fopen(['transport' css '.dat']);
if fid==-1 
   disp(['Could not open section'  css '.dat'])
   return
end
A=' ';
jpwest=[];
jpsouth=[];
cellx=[];
celly=[];
while (isstr(A))
   A=fgetl(fid);
   if (isstr(A))
      B=sscanf(A,'%f%f%f%f');
      cellx=[cellx B(1)];
      celly=[celly B(2)];
      jpwest =[jpwest  B(3)];
      jpsouth=[jpsouth B(4)];
   end
end
fclose(fid);

clf
I=find(depths<.5);
depths(I)=nan;
X=[1:size(depths,2)]+.5;
Y=[1:size(depths,1)]+.5;
P=pcolor(X,Y,depths); set(P,'LineStyle','none')
hold on
% Mark all cells involved
for i=1:prod(size(celly))
   P=patch([cellx(i)-.5 cellx(i)+.5 cellx(i)+.5 cellx(i)-.5], ...
           [celly(i)-.5 celly(i)-.5 celly(i)+.5 celly(i)+.5],[.6 .6 .6]);
   set(P,'FaceAlpha',.7)
end

% Mark vectors from u/v locations
Q1=quiver(cellx-.5,celly   ,jpwest,zeros(size(jpsouth)),.5,'Color','r','LineWidth',2);
Q2=quiver(cellx   ,celly-.5,zeros(size(jpwest)),jpsouth,.5,'Color','g','LineWidth',2);
