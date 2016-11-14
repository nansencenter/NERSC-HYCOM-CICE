function J=gaudy(nbase)
% GAUDY color table
%
%  Similar to the color table used by Rainer Bleck 
%  Can be found in most of the NCAR-based tools of hycom
%

if nargin < 1
   nbase = size(get(gcf,'colormap'),1);
end

num=8;
red=[ 0 0 1 1 0 0 1 1];
grn=[ 0 1 0 1 0 1 0 1];
blu=[ 0 0 0 0 1 1 1 1];

% --- 'expo' controls color interpolation. The larger the value, the
% --- more abrupt the transition between the prescribed color values
expo=1.0;

J=[];

for i=1:nbase
   x=1.+(num-1)*(i-1)/(nbase-1);
   %m=min(num-1,int(x));
   m=min(num-1,floor(x));
   x=x-m;
   if (x<0.5) 
    xx=.5*(2.*x)^expo;
   else
    xx=1.-.5*(2.*(1.-x))^expo;
   end
   r=red(m)*(1.-xx)+red(m+1)*xx;
   g=grn(m)*(1.-xx)+grn(m+1)*xx;
   b=blu(m)*(1.-xx)+blu(m+1)*xx;
   J=[J ; r g b];
end

