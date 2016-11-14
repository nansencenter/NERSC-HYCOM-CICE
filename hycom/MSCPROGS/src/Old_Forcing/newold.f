      SUBROUTINE newold(latn,longn,offi,offj,gridn,ypivo,xpivn,
     +                    ypivn,xold,yold,dxold,dyold)
c
c     Author: R. Bleck and co-workers.
c     Modified by K. Simonsen 2/12 1993.          
c
c --- the purpose of newold is to find the positions of the new (rotated
c --- mercator) grid points in the old (lat/lon) grid. positions are given
c --- in fractional grid units.
c ---
c --- WARNING: The routine is origionaly written for a non rotated grid
c --- and therefor the notation may be a little confusing. However, do
c --- not think about in geogr. pos., only in relativ pos.
c -----------------------------------------------------------------------
c --- input:
c ---   latn,longn: Size of model grid (idm,jdm)
c ---   offi: offset i-index  (h=0,u=.5,v=0) 
c ---   offj: offset j-index  (h=0,u=0,v=.5) 
c ---   gridn:  gridresolution in degrees.
c ---   ypivo:  Model equator meridiane
c ---   xpivn: i-index for the ME
c ---   ypivn: j-index for true equator.
c -----------------------------------------------------------------------
c --- output:
c ---
c ---      xold,yold    --  locations of the new (counterclockwise rotated,
c ---                       mercator) grid points in the old (lat/lon) grid.
c ---      dxold,dyold  --  components of unit vector pointing in the new
c ---                       x direction expressed in terms of the old grid
c ---                       directions.
c -----------------------------------------------------------------------
      REAL offi,offj,gridn,ypivo,xpivn,ypivn
      INTEGER latn,longn 
      REAL xold(latn,longn),yold(latn,longn),
     +         dxold(latn,longn),dyold(latn,longn)
      
      REAL radian,pi,pi2,rjnp,
     +     gridnr,xlatn,xlongn,xlato,xlongo,dist,dum

      INTEGER jnp,i,j
       
      data radian/57.29578/,pi/3.14159265/ 
      pi2 = pi/2.
      gridnr=gridn/radian                       !gridres. in rad       
 
      jnp = INT(ypivn) + int (90./gridn)        !j- index of the Nort Pole.
      rjnp = ypivn + 90./gridn                  !j- index of the Nort Pole.
                                                

      do 2 i=1,latn                             !For all model lat. 
                                                !(long in rot.)
        dist = (-xpivn+float(i)-offi)*gridnr    !Dist. in radians
        xlatn=(2.*atan(exp(dist))-pi/2.)        !Geogr. lat. 
                                                !(=long in rotated)  
        do 2 j=1,longn                          !For all model long
                                                !(lat in rot)
          xlongn=(float(j)-offj-ypivn)*gridnr
          xlato= asin(sin(xlongn)*cos(xlatn))
          dum = sin(xlatn)/cos(xlato)
          dum = MAX(-1., MIN(1., dum))
          xlongo=asin(dum)
          

          if (FLOAT(j).gt.rjnp) then                    !Calc.  
            xlongo = sign(1.,xlato)*pi - xlongo !|pi/2|<long<pi
          endif

          dyold(i,j)=sin(xlongo)*sin(xlongn)    !Grid res. 
          dxold(i,j)=cos(xlongo)/cos(xlatn)     !Grid res. 

          xold(i,j) = xlato*radian              !Geo. lat.
          yold(i,j) = xlongo*radian             !Geo. long
          yold(i,j)= yold(i,j) + ypivo          !ypivo 'equator' meridian.
   
          if (yold(i,j).le. -180.) then         !Check if long< -pi
            yold(i,j) = yold(i,j) + 360.
          endif
          if (yold(i,j).gt. 180.) then          !Check if long > pi
            yold(i,j) = yold(i,j) - 360
          endif 
 2    continue 
      RETURN
      END
