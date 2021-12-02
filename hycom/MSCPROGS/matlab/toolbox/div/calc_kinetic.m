function [tke,mke,eke]=calc_kinetic(uin,vin)
% function [tke,mke,eke]=calc_kinetic(uin,vin)
% 
% routine to calculate total kinetic energy (tke), kinetic energy of the mean flow (mke),
% and eddy kinetic energy (eke).
%
% tke can be split into mke and eke as follows:
% tke = mke + eke
%
% mke reflects the energy from the mean flow and eke reflects the energy from the
% mesoscale variability (eddies).
%
% applying Reynold's averaging: U = Ubar + U'
% where Ubar is the time mean of U, and U' the mesoscale deviations (U' = U - Ubar)
% the same applies for V
% 
% mke = (Ubar*Ubar + Vbar*Vbar)/2
%
% eke = (U'U'bar + V'V'bar)/2
% because U'U'bar (the time mean of the velocity correlations) is difficult to calculate
% first calculate Ubar (time average U) and UUbar (time average U*U) and then
% U'U'bar = UUbar - Ubar*Ubar (same for V)
% 

% calculate Ubar / Vbar
ubar=nanmean(uin,3);
vbar=nanmean(vin,3);

% calculate mke = (Ubar*Ubar + Vbar*Vbar)/2
mke=(ubar.*ubar+vbar.*vbar)./2;

% calculate UUbar and VVbar
uubar=nanmean(uin.*uin,3);
vvbar=nanmean(vin.*vin,3);

% calculate eke =  (U'U'bar + V'V'bar)/2
% where U'U'bar = UUbar - Ubar*Ubar
% and V'V'bar = VVbar - Vbar*Vbar
uprimeuprimebar=uubar-ubar.*ubar;
vprimevprimebar=vvbar-vbar.*vbar;
eke=(uprimeuprimebar+vprimevprimebar)./2;

% calculate tke
tke=mke+eke;
