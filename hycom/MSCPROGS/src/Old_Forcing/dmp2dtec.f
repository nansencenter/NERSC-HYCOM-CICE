      SUBROUTINE dmp2dtec(fld,ifld,dd,nx,ny,ir)
  
      INTEGER nx,ny,ifld(nx,ny),ii,jj
      REAL fld(nx,ny),dd(0:nx+1,0:ny+1)
      CHARACTER*1 ir

      IF(ir.EQ.'r')THEN
      DO i=0,nx+1
       ii=MAX(1,MIN(nx,i))
       DO j=0,ny+1
        jj=MAX(1,MIN(ny,j))
        dd(i,j)=fld(ii,jj)
c        IF(j.NE.jj.OR.i.NE.ii)
c     &    write(*,*)'j,jj',j,jj,dd(i,j),fld(ii,jj)
       ENDDO
      ENDDO
      ELSE IF(ir.EQ.'i') THEN
      DO i=0,nx+1
       ii=MAX(MIN(nx,i),1)
       DO j=0,ny+1
        jj=MAX(MIN(ny,j),1)
        dd(i,j)=FLOAT(ifld(ii,jj))
       ENDDO
      ENDDO
      ELSE
      write(*,*)'Some mismatch in dmt2dtec'
      ENDIF  
      WRITE(1,100)((dd(i,j),i=0,nx+1),j=0,ny+1)
 100  FORMAT(10(1x,e12.6))
      RETURN
      END
