function tstint_fld


XI=1:80;
YI=1:120;
[XI2,YI2]=meshgrid(XI,YI);
XI2=XI2'; % This way XI2(m+1,n) > XI2(m,n)
YI2=YI2'; % This way YI2(m,n+1) > YI2(m,n)

% GLobal vars for nested funcs
opt=2;
k=6; l=8;
x0=31 ; y0=31; r0=1;



old   =f(XI2,YI2);
olddx =fx(XI2,YI2);
olddy =fy(XI2,YI2);
olddxy=fxy(XI2,YI2);
%olddx=zeros(size(XI2));
%olddy=zeros(size(XI2));
%olddxy=zeros(size(XI2));

figure(1); clf; pcolor(XI2,YI2,old); shading flat


% Sampled set
SXI=29:.20:35;
SYI=25:.20:35;
[SXI2,SYI2]=meshgrid(SXI,SYI);
SXI2=SXI2'; % This way XI2(m+1,n) > XI2(m,n)
SYI2=SYI2'; % This way YI2(m,n+1) > YI2(m,n)

for i=1:prod(size(SXI))
for j=1:prod(size(SYI))

   t=floor(SXI(i));
   u=floor(SYI(j));

   g=  [   old(t,u)    old(t+1,u)    old(t,u+1)    old(t+1,u+1)];
   gx= [ olddx(t,u)  olddx(t+1,u)  olddx(t,u+1)  olddx(t+1,u+1)];
   gy= [ olddy(t,u)  olddy(t+1,u)  olddy(t,u+1)  olddy(t+1,u+1)];
   gxy=[olddxy(t,u) olddxy(t+1,u) olddxy(t,u+1) olddxy(t+1,u+1)];

   t=SXI(i)-floor(SXI(i));
   u=SYI(j)-floor(SYI(j));

   solbl(i,j)=bilin(g,t,u);
   c=bicubiccoeffs(g,gx,gy,gxy);
   solbc(i,j)=bicubic(c,t,u);
end
end

figure(1); hold on ; 
P=patch([SXI(1) SXI(end) SXI(end) SXI(1)] , [SYI(1) SYI(1) SYI(end) SYI(end)],'m');
set(P,'FaceColor','none');
set(P,'EDGEColor','m');
set(P,'LineWidth',2);

figure(2); 
old2=f(SXI2,SYI2);
subplot(1,2,1) ; pcolor(SXI2,SYI2,solbl); shading flat
subplot(1,2,2) ; pcolor(SXI2,SYI2,old2-solbl) ; shading flat ; colorbar
figure(3); 
subplot(1,2,1) ; pcolor(SXI2,SYI2,solbc); shading flat
subplot(1,2,2) ; pcolor(SXI2,SYI2,old2-solbc ); shading flat ; colorbar



%%%%%%%%% Nested funcs %%%%%%%%%%%%%%%
   function f=f(x,y)
      if(opt==1) 
         f=  sin(x*2*pi/k) .* sin(y*2*pi/l);
      elseif (opt==2)
         f=exp(- (  (x-x0).^2 + (y-y0).^2 ) /  r0^2);
      end
   end

   function fx=fx(x,y)
      if(opt==1) 
         fx=cos(x*2*pi/k) .* sin(y*2*pi/l)*2*pi/k ;
      elseif (opt==2)
         fx=-exp(- (  (x-x0).^2 + (y-y0).^2 ) /  r0^2) .* ...
           2.*(x-x0)/r0^2;
      end
   end

   function fy=fy(x,y)
      if(opt==1) 
         fy=sin(x*2*pi/k) .* cos(y*2*pi/l)*2*pi/l ;
      elseif (opt==2)
         fy=-exp(- (  (x-x0).^2 + (y-y0).^2 ) /  r0^2) .* ...
           2.*(y-y0)/r0^2;
      end
   end

   function fxy=fxy(x,y)
      if(opt==1) 
         fxy=cos(x*2*pi/k).* cos(y*2*pi/l)*(2*pi)^2/(k*l) ;
      elseif (opt==2)
         fxy=exp(- (  (x-x0).^2 + (y-y0).^2 ) /  r0^2) .* ...
           4.*(y-y0).*(x-x0)/r0^4;
      end
   end


end
