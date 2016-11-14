function coeffs=bicubiccoeffs(g,gx,gy,gxy);
%f  (x,y) = sum_i=1^4 sum_j=1,4 c_ij t^(i-1) t^(j-1)
%fx (x,y) = sum_i=2^4 sum_j=1,4 c_ij (i-1) * t^(i-2) t^(j-1)
%fy (x,y) = sum_i=1^4 sum_j=2,4 c_ij (j-1) * t^(i-1) t^(j-2)
%fxy(x,y) = sum_i=2^4 sum_j=2,4 c_ij (j-1) * (j-1) * t^(i-2) t^(j-2)

% Get matrix
c=zeros(16,16);
c(1,:) =f(0,0);   rhs(1) =g  (1);
c(2,:) =f(1,0);   rhs(2) =g  (2);
c(3,:) =f(0,1);   rhs(3) =g  (3);
c(4,:) =f(1,1);   rhs(4) =g  (4);
c(5,:) =fx(0,0);  rhs(5) =gx (1);
c(6,:) =fx(1,0);  rhs(6) =gx (2);
c(7,:) =fx(0,1);  rhs(7) =gx (3);
c(8,:) =fx(1,1);  rhs(8) =gx (4);
c(9,:) =fy(0,0);  rhs(9) =gy (1);
c(10,:)=fy(1,0);  rhs(10)=gy (2);
c(11,:)=fy(0,1);  rhs(11)=gy (3);
c(12,:)=fy(1,1);  rhs(12)=gy (4);
c(13,:)=fxy(0,0); rhs(13)=gxy(1);
c(14,:)=fxy(1,0); rhs(14)=gxy(2);
c(15,:)=fxy(0,1); rhs(15)=gxy(3);
c(16,:)=fxy(1,1); rhs(16)=gxy(4);

d=inv(c);
%d(1:8)
%d(1,1:8)
%d(1:8,1)

coeffs=d*rhs';

%for l=1:16
%   cst ='';
%   for k=1:16
%      cst= [ cst sprintf( '%2d',d(k,l)) ',' ];
%   end
%   disp(cst);
%end
%ME = MException('MATLAB:TST', 'threw exit');
%ME.throw

%NB - the coeffs matrix may be hardcoded

% Return coeffs for f(x,y)
function c=f(t,u)
   for i=1:4
   for j=1:4
      c(i,j) = t^(i-1) * u^(j-1);
   end 
   end
   c=reshape(c,1,16);


% Return coeffs for fx(x,y)
function cx=fx(t,u)
   cx=zeros(4,4);
   for i=2:4
   for j=1:4
      cx(i,j) = (i-1) * t^(i-2) * u^(j-1);
   end 
   end
   cx=reshape(cx,1,16);


% Return coeffs for f(x,y)
function cy=fy(t,u)
   cy=zeros(4,4);
   for i=1:4
   for j=2:4
      cy(i,j) = (j-1) * t^(i-1) * u^(j-2);
   end 
   end
   cy=reshape(cy,1,16);
   


% Return coeffs for f(x,y)
function cxy=fxy(t,u)
   cxy=zeros(4,4);
   for i=2:4
   for j=2:4
      cxy(i,j) = (j-1) * (i-1) * t^(i-2) * u^(j-2);
   end 
   end
   cxy=reshape(cxy,1,16);


   


   

