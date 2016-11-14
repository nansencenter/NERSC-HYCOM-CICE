function proftest(num)
%function proftest(num)
%
% Uses diag output from vremap_mercator


cnum=num2str(num,'%5.5i');

O=load([ 'oldprof' cnum '.dat']);
N=load([ 'newprof' cnum '.dat']);

whos

figure(1) ; clf; 
stairs(O(:,2),-O(:,1)) ; hold on
stairs(N(:,2),-N(:,1),'Color','r') ; 
set(gca,'XLim',[ min(O(:,2))-.1 max(O(:,2))+.1])

N(:,1)
O(:,1)

figure(2) ; clf; 
plot(O(:,2),-O(:,1)) ; hold on
plot(N(:,2),-N(:,1),'Color','r') ; 
set(gca,'XLim',[ min(O(:,2))-.1 max(O(:,2))+.1])
