      SUBROUTINE filldepth(nx,ny,iflg,depth)

c --- --------------------------------------------
c --- Fill up in sea point, if depth indicate land
c --- --------------------------------------------
      INTEGER nx,ny,iflg(nx,ny),im,ip,jm,jp
      REAL depth(nx,ny),w,d1,d2,d3,d4,fl,fld,flf,flq


      DO i=1,nx
       im=MAX(i-1,1)
       ip=MIN(i+1,nx)
       DO j=1,ny
        jm=MAX(j-1,1)
        jp=MIN(j+1,ny)

        fld=.5+SIGN(.5,depth(i,j)-1.)   != 1 if depth>1m ie. ok
                                        != 0 if depth<1m 
        flf=FLOAT(MIN(iflg(i,j),1))     != 1 if sea point
                                        != 0 if landpoint
        fl=(1.-fld)*flf                 !=1 if fillup is needed
  
        IF(fl.GT..5) THEN
        w=0. 

        flq=.5+SIGN(.5,depth(im,j))
        d1=depth(im,j)*FLOAT(MIN(iflg(im,j),1))*flq
        w=w+FLOAT(MIN(iflg(im,j),1))*flq

        flq=.5+SIGN(.5,depth(ip,j))
        d2=depth(ip,j)*FLOAT(MIN(iflg(ip,j),1))*flq
        w=w+FLOAT(MIN(iflg(ip,j),1))*flq

        flq=.5+SIGN(.5,depth(i,jm))
        d3=depth(i,jm)*FLOAT(MIN(iflg(i,jm),1))*flq
        w=w+FLOAT(MIN(iflg(i,jm),1))*flq

        flq=.5+SIGN(.5,depth(i,jp))
        d4=depth(i,jp)*FLOAT(MIN(iflg(i,jp),1))*flq
        w=w+FLOAT(MIN(iflg(i,jp),1))*flq

        
        d = (d1+d2+d3+d4)/MAX(w,1.)
 
        depth(i,j)=(1.-fl)*depth(i,j)+fl*(d1+d2+d3+d4)/w
        ENDIF
       ENDDO
      ENDDO
      RETURN
      END
        

         
        


  
