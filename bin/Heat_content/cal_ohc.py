"""
Calculation OHC
================

This module contains several python tools used for evaluation of
Oean Heat Content from 2D-array inputs:[depth,grid]

Requires the gsw (Gibbs SeaWater, TEOS-10) package:
    pip install --user gsw
"""

import numpy as np
import gsw

def insect_levs(Nk,Dlev):
   """
   Nk   : dimension in the verticl levels needs to be more than 2
   Dlev : levels on vertical [Nk]

   Return  the insect levels located at the middles.
         Plev[Nk-1]
   """
   # data control
   if  np.size(Dlev) != Nk:
      print("Wrong input the dimesion on vertical!")
      print("broken in cal_ohc")
      quit()
   if Nk>1:
      #Plev=np.zeros(Nk-1,dtype=float)
      Plev=[]
      #print(Dlev)
      for ii in range(Nk-1):
         tmp=(Dlev[ii]+Dlev[ii+1])/2
         Plev.append(tmp)
   else:
      Plev=Dlev[:]
   return Plev


def ohc(Nk,Nx,T,S,Dlev,Dmin,Dmax,Marea,Mdep):
   """
   Nk,Nx: dimensions of T/S in the verticl and the horizontal directions
   T, S : temperature and salinity
   Dlev : grid levels on vertical [Nk]
   Marea: model grid area size by unit m2 [Nx]
   Mdep : grid topography as one limit to the deepest depth [Nx]
   Dmin,Dmax:  pointed the depth limits for upper and bottom

   Return: a vector of ocean heat content on the grid
           Sohc with size of [Nx]
   """
   Cp=3991.87  # Specific heat capacity for ocean
               # J kg-1 K-1

   # data control
   if (Nx != T.shape[1] or Nx != S.shape[1] or
      Nx != Marea.shape[0] or Nx != Mdep.shape[0]):
      print("Check the array dimensions!")
      print("broken in ohc calculation")
      quit()
   elif Nk != Dlev.shape[0]:
      print("Check the depth level dimensions!")
      print("broken in ohc calculation")
      quit()
   if np.sum((Dlev>=Dmin)&(Dlev<=Dmax)) < 1:
      print("Invalid level limits for OHC")
      print("broken in ohc calculation")
      quit()

   print("Calculation the OHC in %.0f~%.0f m "%(Dmin,Dmax))

   # Per-grid-point effective bottom depth
   maxLev=np.minimum(Dmax,Mdep)   # [Nx]

   # Build cell interfaces: midpoints between adjacent levels,
   # with Dmin at top and Dmax at bottom  [Nk+1]
   interfaces=np.empty(Nk+1)
   interfaces[0]=Dmin
   interfaces[1:Nk]=(Dlev[:-1]+Dlev[1:])/2.0
   interfaces[Nk]=Dmax

   # Clip each interface to [Dmin, maxLev] per grid point  [Nk+1, Nx]
   iface=np.clip(interfaces[:,np.newaxis], Dmin, maxLev[np.newaxis,:])

   # Layer thicknesses (negative values clamped to zero)  [Nk, Nx]
   delz=np.maximum(np.diff(iface,axis=0), 0.0)

   # Validity mask: level within depth range and T/S finite  [Nk, Nx]
   valid=((Dlev[:,np.newaxis]>=Dmin) &
          (Dlev[:,np.newaxis]<=maxLev[np.newaxis,:]) &
          np.isfinite(T) & np.isfinite(S))

   # Replace NaN with 0 before density calculation to avoid propagation
   T_safe=np.where(valid, T, 0.0)
   S_safe=np.where(valid, S, 0.0)

   # Compute density for all levels and grid points at once  [Nk, Nx]
   # Pressure (dbar) ≈ depth (m); lon/lat set to 0 as SA correction is small
   p=Dlev[:,np.newaxis]*np.ones_like(T_safe)
   SA=gsw.SA_from_SP(S_safe, p, 0, 0)
   rho=gsw.rho_t_exact(SA, T_safe, p)

   # Integrate rho*T*dz; invalid cells contribute zero
   integrand=np.where(valid, rho*T_safe*delz, 0.0)
   Sohc=np.sum(integrand, axis=0)*Marea*Cp   # [Nx]

   return Sohc
