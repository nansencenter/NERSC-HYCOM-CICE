      SUBROUTINE reddat(nx,ny,nz,field,filen)

c      INCLUDE 'dimensions.h'
c      INCLUDE 'gridpos.h'

      CHARACTER*30 filen
      REAL grdn,ypn,xpn,ypo,field(nx,ny,nz),rfact
     .     ,vmax,vmin,dd,dw,offs
      INTEGER nxl,nyl,nzl,nx,ny,nz,fact,iutil(nx,ny,1)
 
      vmax=-9999.
      vmin= 9999.

      write(*,*)'I read: ',filen 
c      OPEN(50,FILE=filen,FORM='UNFORMATTED')
      OPEN(50,FILE=filen,FORM='FORMATTED')
      READ(50,200)nxl,nyl,nzl,fact
      IF(nx.NE.nxl.OR.ny.NE.nyl.OR.nzl.NE.nz)GOTO 100
      rfact=1./FLOAT(fact)
      j=1
      k=1
      READ(50,201)(iutil(i,j,k),i=1,nxl*nyl*nz)
      CLOSE(50)


      DO i=1,nxl*nyl*nz
       field(i,j,k)=FLOAT(iutil(i,j,k))*rfact
       fl=.5+SIGN(.5,field(i,j,k)+998.)        !=0 if no data
       vmin=MIN(vmin,field(i,j,k)*fl+(1.-fl)*vmin)
       vmax=MAX(vmax,field(i,j,k)*fl+(1.-fl)*vmax)
      ENDDO
      dd=vmax-vmin


      offs=0.
      IF(vmin.LT.0.)offs=-vmin
      IF(vmax.LT.1.AND.dd.LT.1.)dw=100.
      IF(vmax.GE.1.AND.dd.LT.100.)dw=1.  
      IF(dd.GT.100.AND.dd.LT.1000.)dw=.1
      IF(dd.GT.1000.AND.dd.LT.10000.)dw=.01
      IF(dd.GT.10000.AND.dd.LT.100000.)dw=.001
      IF(vmin.GT.200.AND.vmax.LT.400) THEN
        offs=-223.15
        dw=1.
      ENDIF

      WRITE(*,202)filen,vmin,vmax,offs,dw
      idm=nxl
      jdm=nyl
      DO k=1,nz
       WRITE(*,203)k,filen
       DO i=idm,1,-1
        WRITE(*,'(67I2)')(INT((offs+field(i,j,k))*dw),j=jdm-3,2,-1)
       ENDDO
      ENDDO
       
       





 200  FORMAT(4I8)
 201  FORMAT(15I8)
 202  FORMAT('File: ',a10,' Min: ',f9.2,' Max: ',f9.2,
     &       ' Offs: ',f9.2,' Fact ',f9.2)
 203  FORMAT('Month: ',i2,'File: ',a10)
      GOTO 101

 100  WRITE(*,*)'Wrong dimension in file ',filen    
       WRITE(*,*)'File dimension: (nx,ny,nz)',nxl,nyl,nzl
       WRITE(*,*)'Model dimension:(nx,ny,nz)',nx,ny,nz
      CLOSE(50)
      STOP 'reddat'
 101  RETURN 
      END

       
      



