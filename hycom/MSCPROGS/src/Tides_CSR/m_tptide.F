      module m_tptide
      contains
      subroutine tptide( dlat,dlon,time,tide,isdata,u,v,pseudo,
     .                   radiation )
      use m_tptid1
*
* this code is based upon the gstide program used to compute the
*      cartwright&ray geosat tide solution.
*
*
*  name - tptide              name derivation - topex/poseidon derived tide
*
*  function -  to compute the ocean tidal height at a given time
*              and location.
*
*  language - fortran 77
*
*  arguments -
*     name      type  i/o               description
*     ----      ----  ---               -----------
*     dlat       d     i    north latitude (in degrees, -90 to  90)
*
*     dlon       d     i    east longitude (in degrees,   0 to 360).
*
*     time       d     i    desired time, in seconds of mjd.  e.g.,
*                            jan 1, 1986 00:00 is 4011638400.
*
*     tide       d     o    computed tidal height, in cm.
*
*     isdata     l     o    logical denoting whether tide data exist at
*                           desired location.  if false, then tide is
*                           not modified.
*
*     u,v        d     o    the orthoweights interpolated to dlat,dlon
*                           dimensioned (3,2)
*
*     pseudo     l     o    true for c&r style "orthoweights"
*                           this can be passed to admit.f for proper
*                           interpretation of the orthoweights.
*                           set per integer on orthoweight file
*
*     radiation  l     o    true if the correction for radiation
*                           potential is activated.
*                           set per integer on orthoweight file
*
*
*  usage notes -
*     using the user supplied orthoweight file, this routine computes
*     the height of the (ocean + load) tide at the desired location.
*     this quantity is appropriate for use in satellite altimetry.
*     the computed tide is composed of the 30 largest spectral lines
*     within both the diurnal and semidiurnal bands, sufficient for
*     representing all major constituents including nodal modulations.
*
*  processing logic -
*     the tidal height is computed at four grid points surrounding
*     the desired location, with the final result being a bilinear
*     interpolation.
*
*  file references -
*     the orthoweight data are read on initial call to routine:
*     file iunit - global array of orthoweights defining tidal
*                  constants at all valid locations. (see parameter statement)
*
*  important local variables -
*     array uv is loaded with orthoweights as follows:
*         uv(i,j,l,m,n), where
*           i = 1 for u; 2 for v
*           j = 1,2 or 3  (3 complex coeffs per species.)
*           l = tidal species (1 or 2)
*           m = 1 to nx, longitude index number
*           n = 1 to ny, latitude index
*
*     for machines like the cray, the compiler option to ignore double
*         precision should be selected.
*
*  error processing - none.
*
*  technical references -
*     d. cartwright and r. ray, oceanic tides from geosat altimetry,
*        journal of geophysical research, 95, 3069-3090, 1990.
*        (particularly appendix a).
*
*  history -
*   version   date       programmer        change description
*   -------   ----       ----------        ------------------
*     1.0    8/17/90   ray/cartwright  initial version.
*     1.1   12/18/90   ray/cartwright  enhanced documentation;
*                                      use ascii input for portability.
*     1.2   12/31/91   r. ray        changed min w to 0.5 (old value of
*                                    0.2 allowed too much extrapolation)
*     2.0    2/25/94   r. eanes      changed to use t/p derived orthoweights
*     2.1    3/02/94   r. eanes      generalized grid size up to 1 x 1 degree
*                                    and added data statements for the harmonic
*                                    amplitudes and phases (a la the jpl mods)
*                                    interpolate orthoweights instead of tides
*     2.2    3/03/94   r. eanes      return orthoweights (u,v) for other use
*                                    added itype to switch to c&r style orthoweights
*                                    change min w (wmin) back to 0.2.
*     2.3    4/18/94   r. eanes      take out t/p kludge for latitudes > 66.2
*     2.4    5/20/94   r. eanes      change treatment of extrapolating past the southern
*                                    and northernmost grid points.  set wmin to 0.1.
*     2.5                            added r.ray's radiation tide correction
*                                    wmin changed to 0.5 (again)
*     3.0    1/24/95   r. eanes      parameter statements for dimensions and unit
*                                    number.  compact format option for orthoweight file.
*                                    increased dimensions to allow 0.5 x 0.5 deg grid.
*                                    changed wmin to 0.2 to make rsc94b work at gauges.
*

cGE      implicit double precision (a-h,o-z)
      parameter (nxmax=720,nymax=361)
      parameter (iunit=51)
      character*80     title1,title2
      real             uv(2,3,2,nxmax,nymax), undef
      double precision latmin,latmax,lonmin,lonmax
      dimension        u(3,2),v(3,2),iuv(12)
      save
      logical    init,isdata,pseudo,radiation
      data       init/.true./, undef/99999./, eps /1.d-10/, wmin /0.2/

*     on first call, read t/p   -derived ortho-weights
*     ------------------------------------------------
      if (init) then
         init = .false.
         do 4 j=1,nymax
         do 4 i=1,nxmax
            uv(1,1,1,i,j) = undef
    4    continue
         read (iunit,14) title1,title2
         write (6,14) title1,title2
         read (iunit,11) nx,ny,latmin,latmax,lonmin,lonmax,idiff,
     .                   itype,irad,ifmt
c
c...orthoweight file format options
c   idiff - if .eq. 1 then zero tide returned for points outside of domain
c                     and isdata set to true
c   itype - if .eq. 1 then cartwright&ray 1991 pseudo-orthoweights are assumed
c   irad  - if .eq. 1 then the cartwright&ray radiation potential modification is used
c   ifmt  - if .eq. 1 then compact integer format used for orthoweights
c
         if (itype.eq.1) then
            pseudo = .true.  ! c&r style orthoweights
         else
            pseudo = .false. ! corrected orthoweights
         endif
         if (irad.eq.1) then
            radiation = .true.  ! corrected for radiation potential
         else
            radiation = .false. ! not corrected for radiation potential
         endif
         if (nx.gt.nxmax .or. ny.gt.nymax) stop 1   ! must change dimensions
         dx = (lonmax-lonmin)/(nx-1)
         dy = (latmax-latmin)/(ny-1)
         xlats = latmin - dy                    ! southern limit of extrapolation
         xlatn = latmax + dy                    ! northern limit of extrapolation
    8    continue
         if (ifmt.eq.0) then                    ! original c&r format
            read(iunit,12,end=10) i,j,u1,v1,u2,v2,u3,v3,
     .                                u4,v4,u5,v5,u6,v6
            uv(1,1,1,i,j) = u1
            uv(2,1,1,i,j) = v1
            uv(1,2,1,i,j) = u2
            uv(2,2,1,i,j) = v2
            uv(1,3,1,i,j) = u3
            uv(2,3,1,i,j) = v3
            uv(1,1,2,i,j) = u4
            uv(2,1,2,i,j) = v4
            uv(1,2,2,i,j) = u5
            uv(2,2,2,i,j) = v5
            uv(1,3,2,i,j) = u6
            uv(2,3,2,i,j) = v6
         endif
         if (ifmt.eq.1) then                    ! compact integer format
            read(iunit,*,end=10) i,j,iuv
            uv(1,1,1,i,j) =  iuv(1)/1000.d0
            uv(2,1,1,i,j) =  iuv(2)/1000.d0
            uv(1,2,1,i,j) =  iuv(3)/1000.d0
            uv(2,2,1,i,j) =  iuv(4)/1000.d0
            uv(1,3,1,i,j) =  iuv(5)/1000.d0
            uv(2,3,1,i,j) =  iuv(6)/1000.d0
            uv(1,1,2,i,j) =  iuv(7)/1000.d0
            uv(2,1,2,i,j) =  iuv(8)/1000.d0
            uv(1,2,2,i,j) =  iuv(9)/1000.d0
            uv(2,2,2,i,j) = iuv(10)/1000.d0
            uv(1,3,2,i,j) = iuv(11)/1000.d0
            uv(2,3,2,i,j) = iuv(12)/1000.d0
         endif
         go to 8
   10    continue
   11    format(2i4,4f8.3,4i2)
   12    format(2i4,16x,6f8.3,/,24x,6f8.3)
   14    format(a80,/,a80)
      endif
      isdata = .true.
      tide = 0.0d0
      do i=1,3
         do j=1,2
            u(i,j) = 0.d0
            v(i,j) = 0.d0
         enddo
      enddo

*     compute indices for desired position
*     ------------------------------------
      xlat = dlat
      jlat1 = int((xlat - latmin)/dy) + 1
      if (jlat1.eq.0  .and. xlat.gt.xlats) jlat1 = 1        ! allow extrapolation to xlats
      if (jlat1.eq.ny .and. xlat.lt.xlatn) jlat1 = ny - 1   ! allow extrapolation to xlatn
      jlat2 = jlat1 + 1
      if (jlat1.lt.1 .or. jlat2.gt.ny) then
         isdata = .false.

*     if the orthoweight file repesents differences return zero
*     instead of saying there is no data
*     ---------------------------------------------------------
         if (idiff.eq.1) isdata = .true.
         return
      endif
      xlon = dlon
      if (xlon.lt.lonmin) xlon = xlon + 360.d0 - eps   ! eps to avoid problem when xlon=lonmin
      ilon1 = int((xlon - lonmin)/dx) + 1
      if (ilon1.gt.nx) ilon1 = ilon1 - 1               ! avoid problem for dlon=360.
      ilon2 = ilon1 + 1
      if (ilon2.gt.nx) ilon2 = 1                       ! longitude grid must span globe
      if (ilon1.lt.1 .or. ilon1.gt.nx) then
         write (6,*) 'error in tptide: incorrect longitude input'
         write (6,*) 'mjd,dlon,xlon,ilon1,ilon2',
     .                time/86400.d0,dlon,xlon,ilon1,ilon2
         stop 301
      endif

      if (uv(1,1,1,ilon1,jlat1).eq.undef .and.
     *    uv(1,1,1,ilon2,jlat1).eq.undef .and.
     *    uv(1,1,1,ilon1,jlat2).eq.undef .and.
     *    uv(1,1,1,ilon2,jlat2).eq.undef)  then
         isdata = .false.
         if (idiff.eq.1) isdata = .true.
         return
      endif

*  distances from corners in units of cell size
*  --------------------------------------------
      wx1 = ( dx - (xlon - real(ilon1-1)*dx - lonmin) )/dx
      wx2 = 1.d0 - wx1
      wy1 = ( dy - (xlat - real(jlat1-1)*dy - latmin) )/dy
      wy2 = 1.d0 - wy1
      if (wy2.lt.0.) then           ! north corner weights zero in southern extrapolation
         wy1 = wy1 + wy2
         wy2 = 0.
      endif
      if (wy1.lt.0.) then           ! south corner weights zero in northern extrapolation
         wy2 = wy1 + wy2
         wy1 = 0.
      endif

*  interpolation weights:
*  w1,w2,w3,w4 are for northwest,northeast,southeast,southwest corners.
*  --------------------------------------------------------------------
      w1 = wx1*wy2
      w2 = wx2*wy2
      w3 = wx2*wy1
      w4 = wx1*wy1

*  initialize sum of weights and orthoweights
*  ------------------------------------------
      w = 0.d0
      do 50 i=1,3
      do 50 j=1,2
      u(i,j) = 0.d0
   50 v(i,j) = 0.d0

*     interpolate the orthoweights
*     ----------------------------
      if (uv(1,1,1,ilon1,jlat1).ne.undef) then
         do 100 i=1,3
            u(i,1) = uv(1,i,1,ilon1,jlat1)*w4
            v(i,1) = uv(2,i,1,ilon1,jlat1)*w4
            u(i,2) = uv(1,i,2,ilon1,jlat1)*w4
            v(i,2) = uv(2,i,2,ilon1,jlat1)*w4
  100    continue
         w = w4
      endif
      if (uv(1,1,1,ilon1,jlat2).ne.undef) then
         do 200 i=1,3
            u(i,1) = uv(1,i,1,ilon1,jlat2)*w1 + u(i,1)
            v(i,1) = uv(2,i,1,ilon1,jlat2)*w1 + v(i,1)
            u(i,2) = uv(1,i,2,ilon1,jlat2)*w1 + u(i,2)
            v(i,2) = uv(2,i,2,ilon1,jlat2)*w1 + v(i,2)
  200    continue
         w = w + w1
      endif
      if (uv(1,1,1,ilon2,jlat2).ne.undef) then
         do 300 i=1,3
            u(i,1) = uv(1,i,1,ilon2,jlat2)*w2 + u(i,1)
            v(i,1) = uv(2,i,1,ilon2,jlat2)*w2 + v(i,1)
            u(i,2) = uv(1,i,2,ilon2,jlat2)*w2 + u(i,2)
            v(i,2) = uv(2,i,2,ilon2,jlat2)*w2 + v(i,2)
  300    continue
         w = w + w2
      endif
      if (uv(1,1,1,ilon2,jlat1).ne.undef) then
         do 400 i=1,3
            u(i,1) = uv(1,i,1,ilon2,jlat1)*w3 + u(i,1)
            v(i,1) = uv(2,i,1,ilon2,jlat1)*w3 + v(i,1)
            u(i,2) = uv(1,i,2,ilon2,jlat1)*w3 + u(i,2)
            v(i,2) = uv(2,i,2,ilon2,jlat1)*w3 + v(i,2)
  400    continue
         w = w + w3
      endif
      if (w.gt.wmin) then
         do 500 i=1,3
         do 500 j=1,2
         u(i,j) = u(i,j)/w
  500    v(i,j) = v(i,j)/w
         call tptid1( u,v,time,pseudo,radiation,tide )
      else
         isdata = .false.
         if (idiff.eq.1) then
            isdata = .true.
            tide = 0.0d0
            do i=1,3
               do j=1,2
                  u(i,j) = 0.d0
                  v(i,j) = 0.d0
               enddo
            enddo
         endif
      endif
      end subroutine
      end module
