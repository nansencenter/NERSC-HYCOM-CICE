      FUNCTION spherdist(lon1,lat1,lon2,lat2)
c --- -----------------------------------------
c --- Computes the diatance between geo. pos.
c --- lon1,lat1 and lon2,lat2. 
c --- INPUT is in degrees.
c --- -----------------------------------------

      REAL lon1,lat1,lon2,lat2
     +     ,invradian,rearth
     +     ,rlon1,rlat1,rlon2,rlat2
     +     ,x1,y1,z1,x2,y2,z2,dx,dy,dz,dr


      data invradian/0.017453292/,rearth/6371001.0/


      rlon1=lon1*invradian             !lon1 in rad
      rlat1=(90.-lat1)*invradian       !90-lat1 in rad 

      rlon2=lon2*invradian             !lon2 in rad
      rlat2=(90.-lat2)*invradian       !90-lat2 in rad 

      x1= SIN(rlat1)*COS(rlon1)        !x,y,z of pos 1.
      y1= SIN(rlat1)*SIN(rlon1)
      z1= COS(rlat1) 

      x2= SIN(rlat2)*COS(rlon2)        !x,y,z of pos 2.
      y2= SIN(rlat2)*SIN(rlon2)
      z2= COS(rlat2) 

      dx=x2-x1                         !distances in x, y, z 
      dy=y2-y1
      dz=z2-z1

      dr=SQRT(dx*dx+dy*dy+dz*dz)       !distance
 
      spherdist=dr*rearth

      RETURN
      END
