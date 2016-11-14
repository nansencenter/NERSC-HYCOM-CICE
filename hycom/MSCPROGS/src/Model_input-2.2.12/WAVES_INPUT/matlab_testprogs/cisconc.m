function col=cisconc(m)
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

c=[205 205 255;
    80  80 255;
     0 180   0;
   255 240   0;
   255 120   0;
   180   0   0;
   255 180 180]/255;
   
%c=[  0 160 255;
%   250 255 255;
%   187   0   0;
%   255 120   0;
%   255 255  50;
%     0 120   0;
%   228 255 220]/255;
   
for i=1:3;
 col(:,i)=interp1([0 0.05 0.15 0.45 0.80 0.92 1],c(:,i)',linspace(0,1,100))';
end;
