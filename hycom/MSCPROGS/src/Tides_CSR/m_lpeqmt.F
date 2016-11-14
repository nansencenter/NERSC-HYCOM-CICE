      module m_lpeqmt
      contains
      subroutine lpeqmt( ts, dlat, tlp )
*
*  function - computes the long-period equilibrium ocean tides.
*
*  arguments -
*     name      type  i/o               description
*     ----      ----  ---               -----------
*     ts         d     i   modified julian day, in seconds, denoting
*                          time at which tide is to be computed.
*
*     dlat       d     i   latitude in degrees (positive north) for the
*                          position at which tide is computed.
*
*     tlp        d     o   computed long-period tide, in centimeters.
*
*
*  processing logic -
*     fifteen tidal spectral lines from the cartwright-tayler-edden
*     tables are summed over to compute the long-period tide.
*
*  technical references -
*     cartwright & tayler, geophys. j. r.a.s., 23, 45, 1971.
*     cartwright & edden, geophys. j. r.a.s., 33, 253, 1973.
*
*  history -
*    version   date       programmer         change description
*    -------   ----       ----------         ------------------
*      1.0   11/27/90     d. cartwright    documented by r. ray
*
*
cGE      implicit double precision (a-h,o-z)
      dimension  phc(4),dpd(4),shpn(4)
      save
*
      parameter (rad= 0.017453292519943 )
      parameter (psol=283.*rad)       ! solar perigee held constant
      data phc/290.21,280.12,274.35,343.51/,
     *     dpd/13.1763965,0.9856473,0.1114041,0.0529539/,
     *     tc/40431744.d2/
*
*     compute 4 principal mean longitudes in radians at time td
*     ---------------------------------------------------------
      td = (ts - tc)/864.d2             ! days past 1 jan 87 00:00
      do 1 n=1,4
         ph = phc(n) + td*dpd(n)
         shpn(n) = rad*mod(ph,360.0)
    1 continue
*
*     assemble long-period tide potential from 15 cte terms > 1 mm.
*     nodal term is included but not the constant term.
*     -------------------------------------------------
      zlp = 2.79 * cos( shpn(4) )
     *     -0.49 * cos( shpn(2) - psol )
     *     -3.10 * cos( 2.*shpn(2) )
      ph = shpn(1)
      zlp = zlp - 0.67 * cos( ph - 2.*shpn(2) + shpn(3) )
     *          -(3.52 - 0.46*cos(shpn(4))) * cos( ph - shpn(3) )
      ph = ph + shpn(1)
      zlp = zlp - 6.66 * cos( ph )
     *          - 2.76 * cos( ph + shpn(4) )
     *          - 0.26 * cos( ph + 2.*shpn(4) )
     *          - 0.58 * cos( ph - 2.*shpn(2) )
     *          - 0.29 * cos( ph - 2.*shpn(3) )
      ph = ph + shpn(1)
      zlp = zlp - 1.27 * cos( ph - shpn(3) )
     *          - 0.53 * cos( ph - shpn(3) + shpn(4) )
     *          - 0.24 * cos( ph - 2.*shpn(2) + shpn(3) )
*
*     multiply by gamma_2 * sqrt(5/4 pi) * p20(lat)
*     ---------------------------------------------
      s = sin(dlat*rad)
      tlp = 0.437*zlp*(1.5*s*s - 0.5)
      end subroutine
      end module
