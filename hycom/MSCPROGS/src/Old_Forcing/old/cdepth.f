      SUBROUTINE cdepth(nx,ny,depth,iflg,dh,dd,iw)

      INTEGER iw,nx,ny
     &    ,iflg(nx,ny),L,M

      REAL depth(nx,ny),dd(nx,ny)
     &    ,dh(0:nx,0:ny) 

      IF(iw.GT.0)THEN             !Load Bathymetry
       WRITE(*,*)'Load bathymetry from file '
       OPEN(69,FILE='Data/bath.dat',ERR=250,FORM='formatted')
       READ(69,'(2i4)',ERR=250) L,M
       READ(69,'(15e14.7)',ERR=250)((dh(i,j),i=0,L),j=0,M)
       CLOSE(69)


       IF(iw.LT.9) CALL rmask(iflg,nx,ny)     !Load mask
      ENDIF


      IF(iw.EQ.0)THEN
       DO i=1,nx
        DO j=1,ny
         depth(i,j)=0.
         iflg(i,j)=1
        ENDDO
       ENDDO
       WRITE(*,*)'Depth=0  and iflg=1 in the entire domain'


      ELSE IF(iw.EQ.1)THEN
       WRITE(*,*)'Load bathymetry from file fort.63',nx,ny
       DO i=1,nx
        DO j=1,ny
         IF(iflg(i,j).EQ.2)iflg(i,j)=1    !Include flag=2 areas
         depth(i,j)=dh(i,j)*MIN(FLOAT(iflg(i,j)),1.)
        ENDDO
       ENDDO
                                  !Fillup missing values
       CALL filldepth(nx,ny,iflg,depth)

      ELSE IF(iw.EQ.3)THEN
       DO i=1,nx
        DO j=1,ny
         IF(iflg(i,j).EQ.2)iflg(i,j)=0    !Exclude flag=2 areas
         depth(i,j)=dh(i,j)*MIN(FLOAT(iflg(i,j)),1.)
        ENDDO
       ENDDO
                                  !Ensure numerical stability
       CALL gridstab(nx,ny,iflg,depth)
                                  !Fillup missing values
       CALL filldepth(nx,ny,iflg,depth)
                                  !Run Shapirofilter
       CALL shap2d(depth,dd,nx,ny,iflg)

      ELSE IF(iw.EQ.4)THEN
       DO i=1,nx
        DO j=1,ny
         IF(iflg(i,j).EQ.2)iflg(i,j)=1    !Include flag=2 areas
         depth(i,j)=dh(i,j)*MIN(FLOAT(iflg(i,j)),1.)
        ENDDO
       ENDDO
                                  !Ensure numerical stability
       CALL gridstab(nx,ny,iflg,depth)
                                  !Fillup missing values
       CALL filldepth(nx,ny,iflg,depth)
                                  !Run Shapirofilter
       CALL shap2d(depth,dd,nx,ny,iflg)

      ELSE IF(iw.EQ.5)THEN
       DO i=1,nx
        DO j=1,ny
         IF(iflg(i,j).EQ.2)iflg(i,j)=1    !Include flag=2 areas
         depth(i,j)=dh(i,j)*MIN(FLOAT(iflg(i,j)),1.)
        ENDDO
       ENDDO
                                  !Fillup missing values
       CALL filldepth(nx,ny,iflg,depth)
                                  !Ensure numerical stability
       CALL gridstab(nx,ny,iflg,depth)
                                  !Run Shapirofilter
       CALL shap2d(depth,dd,nx,ny,iflg)
       CALL shap2d(depth,dd,nx,ny,iflg)
                                  !Fillup missing values
       CALL filldepth(nx,ny,iflg,depth)

       WRITE(*,*)'-----------------------------------------'
       WRITE(*,*)'Special treatment in the Greenland Strait'
       WRITE(*,*)'and the Faroe Bank Channel. Only valid for'
       WRITE(*,*)'ESOP-1 version'
       WRITE(*,*)'-----------------------------------------'
   
       j=41                      !Greenland Strait
       DO i=66,68 
        depth(i,j)=MAX(depth(i,j),580.)
       ENDDO
       j=42
       i=42
       depth(i,j)=MAX(depth(i,j),300.)

       i=62              !Faroe Bank Channel
       DO j=18,24
        write(*,*)i,j,depth(i,j)
        depth(i,j)=MAX(depth(i,j),800.)
        write(*,*)depth(i,j)
        
       ENDDO
       i=63              
       DO j=17,25
        depth(i,j)=MAX(depth(i,j),600.)
       ENDDO
       i=61              !Wville Thomsen Ridge           
       DO j=17,22
        depth(i,j)=MIN(depth(i,j),550.)
       ENDDO




      ELSE IF(iw.EQ.9)THEN
       WRITE(*,*)'Load bathymetry from file fort.63'
                                  
                                  !Generate mask file from depthfile
       DO i=nx,1,-1
        write(10,'(71i1)')(INT(.5+SIGN(.5,dh(i,j)-10.)),j=ny,1,-1)
       ENDDO
       WRITE(*,*)'Mask is dumped to file fort.10'
      ENDIF
   
      GOTO 900     
 250  WRITE(*,*)'Error when reading bathymetry in cdepth.f'
      STOP'cdepth'
 900  RETURN
      END
