
function satvap_formulas()


t=273.16-40:1:273.16+40;

satvappw1=satvappw(t);
satvappw2=satvappw_gill(t);
satvappw3=satvappw_jacob(t);

plot(t,satvappw1);

hold on;

plot(t,satvappw2,'Color','red');
plot(t,satvappw3,'Color','green');



% #######################################################################
% Saturation vapor pressure over water
% #######################################################################
   function satvappw=satvappw(tx)
   satvappw=611.*10.^(7.5.*(tx-273.16)./(tx-7.50)) ; % tx in K



% #######################################################################
% Saturation vapor pressure over water
% #######################################################################
   function satvappw_gill=satvappw_gill(tx)
   tx=tx-273.16;
   satvappw_gill=10.^( (0.7859 + 0.03477.*tx)./(1 + 0.00412.*tx) );
   satvappw_gill= satvappw_gill * 100.;



  function satvappw_jacob=satvappw_jacob(t)
% This function calculates the saturation vapour pressure
% [Pa] from the temperature [deg K].
% Modified: Anita Jacob, June '97
%
% Input: t: temperature [deg K]
% Output: satvap: saturation vapour pressure at temp. t
%
% es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual


%      implicit none
%      real, intent(in):: t
%      real :: aa,bb,c1,c3,c4,t00,cc
%     data c1/610.78/,t00/273.16/
      c1=610.78;
      t00=273.16;

      %if (t < t00) then
      %   c3 = 21.875
      %   c4 = 7.66
      %else
         c3 = 17.269
         c4 = 35.86
      %endif
      aa = c3 .* (t - t00);
      bb = t - c4;
      
      cc=aa./bb;
      %if (cc < -20.0) then
      %   satvap=0.0
      %else
      %   satvap = c1 * exp(aa/bb)
      %endif
      satvappw_jacob = c1 .* exp(aa./bb);
