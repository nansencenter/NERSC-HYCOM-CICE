      SUBROUTINE fillall(fd,nx,ny,nt,mean)
c --------------------------------------------------------------------
c --- Fill up in model grid points with no data.
c --- Author Knud Simonsen, NERSC, 2/12 1993.
c --- ----------------------------------------------------------------
c --- The routine is based on the routine FILLUP in the OPYC code.
c --- In case of no data then the value in the same point in the field
c --- (month) just before is inserted as an initial value for the 
c --- interpolation. If this point do not contains data, then the "mean"
c --- value is inserted as initial interpolation value.
c ---
c --- Input:  fd: Datafile
c ---         nxm,nym: Grid size.
c ---         nt: Number of fields (fi 12 months: nt=12) 
c ---         mean: A value which is set into missing points before
c ---               the interpolations starts.
c --- Output: fd: filled datafile
c --------------------------------------------------------------------     
c      PARAMETER(nxl=100,nyl=137)
c      PARAMETER(nxl=121,nyl=71)
      PARAMETER(nxl=291,nyl=127)
c                                       !nyl=jdm in newf4.f
      REAL fd(nx,ny,nt),we,ww,wn,ws,mean,mean1
      INTEGER nx,ny,i,ii,j,mn,ms,mw,me,k,l,f,mo,
     +        ip,im,jp,jm,mb
      INTEGER kflg(nxl,nyl),skflg
      
      IF(nx.NE.nxl.OR.ny.NE.nyl) THEN
        write(*,*)'nx,ny = ',nx,ny  
        write(*,*)'nxl,nyl = ',nxl,nyl
        STOP 'Check FILLALL'
      ENDIF
                             !Check also "write(1,"-'s
  
      f=10                   !Factor, which the points with
                             !data are weighed more than 
                             !cells without data.
      mb=nt                  !Save index for the last 'month'
      DO 10 mo=1,nt          !For all months do
      skflg=0
      DO i=1,nx
       DO j=1,ny
        IF (fd(i,j,mo).LT.-998.)THEN   !Identify no data grid
          kflg(i,j) = 0                !No data exist flag
          IF(fd(i,j,mb).GT.-998.) THEN !If data exsist in the 
            mean1= fd(i,j,mb)          !field just before, then 
          ELSE                         !this value is used, else 
            mean1=mean                 !the mean value is used
          ENDIF 
          fd(i,j,mo) = mean1
        ELSE
          kflg(i,j) = 1                !Data exist flag 
        ENDIF
       ENDDO
      ENDDO

      DO ii=1,50                       !Number of iter.
       DO k=1,-1,-2                    !Def dir. of j-iter.
        ms = (1+k+(1-k)*ny)/2          !j-start (1 or ny)
        mn = ((1+k)*ny+1-k)/2          !j-stop  (ny or 1) 
        DO l=1,-1,-2                   !Def dir of i-iter
         mw = (1+l+(1-l)*nx)/2         !i-start (1 or nx)
         me = ((1+l)*nx+1-l)/2         !i-stop  (nx or 1) 
            
         DO j=ms,mn,k              !For all MD do
          DO i=mw,me,l         
                                       !If no data exist 
           IF(kflg(i,j).eq.0.) THEN 
                                       ! Ensure that data are
                                       ! witin the MD:
            ip=MIN(i+1,nx)             ! =i+1
            im=MAX(i-1,1)              ! =i-1
            jp=MIN(j+1,ny)             ! =j+1
            jm=MAX(j-1,1)              ! =j-1   
                                       ! Define weights for 
                                       ! surrounding points..
                                       ! If land: w=0
                                       ! If data: w=f+1
                                       ! If no data: w=1
            we=FLOAT(kflg(ip,j)*f+1)
            ww=FLOAT(kflg(im,j)*f+1)
            wn=FLOAT(kflg(i,jp)*f+1)
            ws=FLOAT(kflg(i,jm)*f+1)
            fd(i,j,mo)=(we*fd(ip,j,mo) 
     +                + ww*fd(im,j,mo)
     +                + wn*fd(i,jp,mo) 
     +                + ws*fd(i,jm,mo))
     +                 /MAX((we+ww+wn+ws),0.00001)
           ENDIF
          ENDDO         ! i-loop
         ENDDO         ! j-loop
        ENDDO         ! l-loop
       ENDDO         ! k-loop 
      ENDDO
      mb=mo
  10  CONTINUE  !mo-loop
      RETURN
      END  
