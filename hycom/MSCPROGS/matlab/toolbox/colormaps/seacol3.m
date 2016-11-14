function col=anom(m)
%JET    Variant of HSV.
%   BCM(M), a variant of HSV(M), is the colormap used in BCM
%   BCM, by itself, is the same length as the current colormap.
%   Use COLORMAP(BCM).
%
%   See also HSV, HOT, JET, PINK, FLAG, COLORMAP, RGBPLOT.

%   Tore Furevik 4/9/2001 

if nargin < 1,
  m=size(get(gcf,'colormap'),1);
end

%c=[  0   0 255;
%   100 100 255;
%   100 255 255;
%   150 255 150;
%   255 255 255;
%   255 255 127;
%   255 220   0;
%   255 127   0;
%   150   0   0]/255;
c=[  %0   0 200;
     0 100 255;
     0 255 255;
   100 255 100;
   250 255 250]/255;

for i=1:3;
 col(:,i)=interp1(linspace(0,1,length(c)),c(:,i)',linspace(0,1,m))';
end;
