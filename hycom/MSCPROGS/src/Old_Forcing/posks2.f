      subroutine posks2(llat,llong,offset,ri,rj,
     +                 gridn,ypivn,xpivn,ypivo)
      real llat,llong,ri,rj,ypivo,xpivn,ypivn,gridn
     +     radian,pi

      integer icp
c
c --- Author: Knud Simonsen, NERSC, nov 1993.
c --- This routine is based on the routine oldnew.f in the Bleck code,
c --- whice is reversed.
c
c --- This routine find the i,j - model coordinates (ri,rj) in a
c --- rotatet Mercator model grid. 
c
c --- Input: lat: latitude in degrees.
c ---        long: longitude in degrees.
c ---        offset =  0.0:  mass      points
c ---        offset = -0.5:  vorticity points
c
c --- output : ri: i-model coordinate
c ---          rj: j-model coordinate. 

      radian=57.2957795
      pi=3.1415926536

      gridnr=gridn/radian                       !gridres. in rad       
      
      xlongo = llong-ypivo           !ypivo is the meridiane, which
                                     !applied as 'model equator'
      if (xlongo.le.-180.) then      !Ensure that -180<xlongo<180
           xlongo = xlongo +360.    
      endif
      if (xlongo.gt.180.) then
           xlongo = xlongo -360.
      endif 
    
                                     !If crossing the meridians normal
                                     !to the ME then icp=2, else icp=0.
      icp = AINT(1.+SIGN(1.,(ABS(xlongo)-90.)))
     

      xlongo = xlongo/radian 
      xlato  = llat/radian          !xlongo and xlato in rad.
      dum = -sin(xlongo)*cos(xlato)
      dum = MAX(-1.,MIN(1.,dum))
      xlatn = asin(dum)
      dum = sin(xlato)/cos(xlatn)
      dum = MAX(-1.,MIN(1.,dum))
      xlongn = asin(dum)
      rj = ypivn - offset + xlongn/gridnr  !j - index
      if (icp.GT.1) then
        rj = 180./gridn -(xlongn/gridnr- ypivn + offset)
      endif
                                      !ypivn is the j-index of the
                                      ! true equator   

      dist1 =alog(tan((2.*xlatn+pi)/4.))/gridnr  !Calc. dist.
      ri = xpivn - dist1 - offset     !i - index
                                      !xpivn is the i-index of the 
                                      !model equator.

      return
      end
