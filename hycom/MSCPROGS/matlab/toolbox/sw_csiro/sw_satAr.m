%$$$ 
%$$$ #undef __PR
%$$$ #include "VARIANT.h"

function c = sw_satAr(S,T)

% SW_SATAr   Satuaration of Ar in sea water
%=========================================================================
% sw_satAr $Revision: 1.1 $  $Date: 1998/04/22 02:15:56 $
%          Copyright (C) CSIRO, Phil Morgan 1998.
%
% USAGE:  satAr = sw_satAr(S,T,P)
%
% DESCRIPTION:
%    Solubility (satuaration) of Argon (Ar) in sea water
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (IPTS-68)]
%
% OUTPUT:
%   satAr = solubility of Ar  [ml/l] 
% 
% AUTHOR:  Phil Morgan 97-11-05  (morgan@ml.csiro.au)
%
%$$$ #include "disclaimer_in_code.inc"
%
% REFERENCES:
%    Weiss, R. F. 1970
%    "The solubility of nitrogen, oxygen and argon in water and seawater."
%    Deap-Sea Research., 1970, Vol 17, pp721-735.
%=========================================================================

% CALLER: general purpose
% CALLEE: 

%$$$ #ifdef VARIANT_PRIVATE
%$$$ %***********************************************************
%$$$ %$Id: sw_satAr.M,v 1.1 1998/04/22 02:15:56 morgan Exp $
%$$$ %
%$$$ %$Log: sw_satAr.M,v $

%$$$ %
%$$$ %***********************************************************
%$$$ #endif

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=2
   error('sw_satAr.m: Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   error('sw_satAr: S & T must have same dimensions')
end %if

% IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
Transpose = 0;
if ms == 1  % row vector
   T       =  T(:);
   S       =  S(:);   
   Transpose = 1;
end %if

%------
% BEGIN
%------

% convert T to Kelvin
T = 273.15 + T; 

% constants for Eqn (4) of Weiss 1970
a1 = -173.5146;
a2 =  245.4510;
a3 =  141.8222;
a4 =  -21.8020;
b1 =   -0.034474;
b2 =    0.014934;
b3 =   -0.0017729;

% Eqn (4) of Weiss 1970
lnC = a1 + a2.*(100./T) + a3.*log(T./100) + a4.*(T./100) + ...
      S.*( b1 + b2.*(T./100) + b3.*((T./100).^2) );

c = exp(lnC);

return

