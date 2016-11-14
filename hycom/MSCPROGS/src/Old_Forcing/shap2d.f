      SUBROUTINE shap2d(yd,yf,nx,ny,iflg)
c --- -----------------------------

c --- ---------------------------------------------
c --- INPUT: 
c ---     nx,ny: dim. of data array
c ---     yd: unfiltered data
c ---     yf: dummy array
c --- OUTPUT:
c ---     yd: filtered data    
c --- ---------------------------------------------
      PARAMETER(n=8)    !order of the Shapiro filter
c --- ----------------------------------------------
       REAL n2f,af,jfac,fact(2*n+1)
     &     ,yd(nx,ny),yf(nx,ny)
     &     ,yff,fact1

       INTEGER  i,j,m,nx,ny,jn
     &  ,n2,imin,imax,im,imm,idm,imn,immn,idmm,immmn
     &  ,iflg(nx,ny)

       DO i=1,nx
        DO j=1,ny
         yd(i,j)=MAX(yd(i,j),20.)
        ENDDO
       ENDDO
     

 
       n2=2*n                !Points involved in the filter
       a=MOD(n-1,2)
       fact1=(1.-2.*a)/2.**n2

       n2f=1.                 !Compute 2n!
       DO m=1,n2
        n2f=n2f*FLOAT(m)
       ENDDO
       fact(1)=1.
       m= 0

       jn= 1
       jfac=1.                !Compute j!
       af = 1.                ! 2n!/(2n-m)!
       DO m=1,n2
        jfac=jfac*FLOAT(m)    !Compute j!
        jn=-1*jn              !Compute(-1)^j
        af=af*FLOAT(n2-m+1)
        fact(m+1)=FLOAT(jn)* af/jfac
       ENDDO

       imin=2
       imax=nx-1
       jmin=2
       jmax=ny-1
c       GOTO 800
       DO j=jmin,jmax         !Perform filtering in x direction
        DO i=imin,imax
         yff=0.
         DO m=0,n2
          im=i+n-m 
          imm=MIN(imax,MAX(imin,im))
          idm=imm-im
          imn=im+2*idm 
          immn=MIN(imax,MAX(imin,im+2*idm))
          idmm=immn-imn
          immmn= MIN(imax,MAX(imin,imn+2*idmm))
          yff=yff+ fact(m+1)*yd(immmn,j)
         ENDDO   
         yf(i,j)=yd(i,j)+fact1*yff
        ENDDO
       ENDDO
c       GOTO 900
 800   CONTINUE

       DO i=imin,imax         !Perform filtering in x direction
        DO j=jmin,jmax
         yff=0.
         DO m=0,n2
          im=j+n-m 
          imm=MIN(imax,MAX(imin,im))
          idm=imm-im
          imn=im+2*idm 
          immn=MIN(imax,MAX(imin,im+2*idm))
          idmm=immn-imn
          immmn= MIN(imax,MAX(imin,imn+2*idmm))
          yff=yff+ fact(m+1)*yf(i,immmn)
         ENDDO   
         yd(i,j)=yf(i,j)+fact1*yff
        ENDDO
       ENDDO
       DO i=1,nx
        DO j=1,ny
         yd(i,j)=MAX(yd(i,j),20.)*MIN(FLOAT(iflg(i,j)),1.)
        ENDDO
       ENDDO
 900   RETURN
       END 
      
