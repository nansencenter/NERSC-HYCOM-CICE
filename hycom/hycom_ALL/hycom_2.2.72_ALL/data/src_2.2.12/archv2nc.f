      PROGRAM a2nc
      use globals
      use date_sub
      
      IMPLICIT NONE
      REAL,DIMENSION(4096)   :: PAD
      INTEGER                :: IOS,locate,ix1,ix2,jy1,jy2   
      INTEGER                :: i,j,k,ii,jj,kk,kz,yy2,mm2,dd2
      LOGICAL                :: transport,dim3,mthin
      LOGICAL,parameter      :: LSPVAL=.TRUE.
      REAL, parameter        :: InvG=1.01978e-4,CtoM=0.1014
      REAL                   :: THBASE
      REAL                   :: xStart,xEnd,yStart,yEnd,hmin,hmax
      INTEGER                :: NPAD,FldIndex,Num3DFlds,itype,LTIndex
      INTEGER                :: iidm,jjdm,kstart,kend,ni
      INTEGER                :: uBaroIndex,vBaroIndex
      CHARACTER*240          :: datafile,grid,depthFile,sfile
      CHARACTER*240          :: tmp,tmp1,tmp2
      CHARACTER(LEN=100)     :: concat,outputFormat   
      INTEGER                :: year, month,day
      CHARACTER (LEN=8)      :: TODAYS_DATE 

      CALL DATE_AND_TIME (DATE=TODAYS_DATE)
      layers=(/(i,i=1,20)/)

      

!!!!!!!!!!!! READ INPUT PARAMETERS FROM SCRIPT !!!!!!!!!!!
      READ *,title                    ! READ netcdf title
      READ *,institution              ! READ source institution
      READ *,expt                     ! READ EXPT no   
      READ *,history                  ! creation history
      READ *,domain_name              ! domain_name
      READ *,field_date               ! date of field
      READ *,julian_date_unit         ! Start of the NetCDF time axis
      READ *,z_positive_direction     ! down or up
      READ *,conventions              ! dataset conventions
      READ *,idm                      ! Longitudinal dimension
      READ *,jdm                      ! Latitudinal  dimension
      READ *,kdm                      ! Veritcal dimension
      READ *,domain                   ! subset or whole domain
      READ *,xstart                   ! Starting x coordinate for subsampling
      READ *,xend                     ! Ending x coordinate for subsampling
      READ *,ystart                   ! Starting y coordinate for subsampling
      READ *,yend                     ! Ending y coordinate for subsampling
      READ *,kstart                   ! Starting vertical coordinate
      READ *,kend                     ! Ending vertical coordinate
      READ *,FldIndex                 ! record number of the variable in the archive file 
      READ *,LTIndex                  ! record number of the layer thickness variable in the file
      READ *,uBaroIndex               ! record number of barotropic uvel
      READ *,VBaroIndex               ! record number of barotropic vvel
      READ *,Num3DFlds                ! number of 3D fields
      READ *,year                     ! Year from input file 
      READ *,month                    ! Month from Input file
      READ *,day                      ! day from Input file
      READ *,name                     ! Name of the variable
      READ *,long_name                ! long name of the variable
      READ *,standard_name            ! standard name of the variable      
      READ *,units                    ! units of the variable
      READ *,unit_long                ! long name of units      
      READ *,dim3                     ! is it a 3D variable?
      READ *,fname                    ! output file name
      READ *,grid                     ! regional grid files
      READ *,depthFile                ! regional depth files
      READ *,sfile                    ! sigma file
      READ *,datafile                 ! name of the file to read data from 


!! READ sigma values from a sigma file
      

       open (unit=ni,file=trim(sfile),form='formatted', status='old')

       read(ni,'(f6.3)') thbase
       allocate(sigma(kdm))
       do k=1,kdm
       read(ni,'(f6.3)') sigma(k)
       enddo



!! calculate julian day of the archive file starting from 1900/12/31

      field_julian_date=ndays(1900,12,31,month,day,year)

!! set history
      history = "converted from mean archive file on: "  // TODAYS_DATE 


!! set filler/spval
       filler=2.00**100
       spval=filler
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Calculate the number of elements padded!!!!!!!!!!!!!!!!!!!!!!!!
       NPAD = 4096 - MOD(IDM*JDM,4096)
       IF     (NPAD.EQ.4096) THEN
        NPAD = 0
       ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Allocate and Read the fields from regional grid


 
       ALLOCATE(lon(IDM,JDM))
       ALLOCATE(lat(IDM,JDM))
       
       SELECT CASE (name)

       CASE('ubaro')
       CALL RAW(lon,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,5)
       CALL RAW(lat,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,6)
        
       CASE('umix')
       CALL RAW(lon,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,5)
       CALL RAW(lat,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,6)
       
       CASE('uvel')
       CALL RAW(lon,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,5)
       CALL RAW(lat,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,6)
       
       CASE('vbaro')
       CALL RAW(lon,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,7)
       CALL RAW(lat,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,8)
       
       CASE('vmix')
       CALL RAW(lon,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,7)
       CALL RAW(lat,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,8)
      
       CASE('vvel')
       CALL RAW(lon,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,7)
       CALL RAW(lat,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,8)

       CASE DEFAULT       


       CALL RAW(lon,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,1)
       CALL RAW(lat,IDM,JDM,PAD,NPAD, LSPVAL,SPVAL,grid,2)

       END SELECT


       !lon=360+lon
       
       ix1=locate(lon(:,1),xstart,IDM)
       ix2=locate(lon(:,1),xend,IDM)
       jy1=locate(lat(1,:),ystart,JDM)
       jy2=locate(lat(1,:),yend,JDM)


       ii=ix2-ix1;jj=jy2-jy1

       iidm=ii;jjdm=jj

       !print *,ix1,ix2,jy1,jy2,iidm,jjdm
       !  stop


       IF(domain.eq.'full') then
       ix1=1;jy1=1
       iidm=IDM;jjdm=JDM
       ii=idm-1;jj=jdm-1  
       ENDIF
       
        ALLOCATE(A2D(IDM,JDM))
        ALLOCATE(A(IDM,JDM,KDM))
        ALLOCATE(depth(IDM,JDM))
        ALLOCATE(ubaro(IDM,JDM))
        ALLOCATE(vbaro(IDM,JDM))
        ALLOCATE(p(IDM,JDM,kDM+1))
        ALLOCATE(dp(IDM,JDM,KDM))
        ALLOCATE(AZ(IDM,JDM,Kz))
        ALLOCATE(WK1(IDM,JDM,KDM))
 

      
               
       CALL RAW(depth,IDM,JDM,PAD,NPAD,LSPVAL,SPVAL,depthFile,1)
       CALL RAW(ubaro,IDM,JDM,PAD,NPAD,LSPVAL,SPVAL,datafile,uBaroIndex)
       CALL RAW(vbaro,IDM,JDM,PAD,NPAD,LSPVAL,SPVAL,datafile,vBaroIndex)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
       
        
       IF (dim3 .eq. .true.) then
       print *, "output in layers"       

        DO k=1,KDM

           CALL RAW(A(:,:,k),IDM,JDM,PAD,NPAD,
     &     LSPVAL,SPVAL,datafile,((k-1)*Num3DFlds+FldIndex))

          
           CALL RAW(dp(:,:,k),IDM,JDM,PAD,NPAD,
     &         LSPVAL,SPVAL,datafile,((k-1)*Num3DFlds+LTIndex))


         enddo

        do k=1,kdm
        DO j=1,JDM
        DO i=1,IDM

        if ( depth(i,j) .lt. spval) then
            dp(i,j,k)=dp(i,j,k)/9806.
               p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
         else
                p(i,j,k+1)=spval
        endif
         ENDDO
         ENDDO
         enddo           
       
         
              do j=1,jdm
                DO i=1,idm
                  p(i,j,1)=0.
             enddo
             enddo
!         do k=1,kdm
! 
!        where (p(:,:,k)+2.0 .gt. p(:,:,k+1))
!         A(:,:,k)=spval
!         endwhere
!         enddo


       
       else 
             
        CALL RAW(A2d,IDM,JDM,PAD,NPAD,LSPVAL,SPVAL,datafile,FldIndex)
       
       END IF 
   

          SELECT CASE (name)

!! note in the following INVG is 1./9806.
!! CtoM is to convert ssh to meters 1.0/9.806.



          CASE('ssh')
          WHERE(A2D < SPVAL)  A2D=A2D*CtoM
         

          CASE('mldpth')
          WHERE(A2D < SPVAL)  A2D=A2D*InvG          
 
          CASE('bldpth')
          WHERE(A2D < SPVAL)  A2D=A2D*InvG

          CASE('umix')
          WHERE(A2D < SPVAL)  A2D=A2D+ubaro
          
          CASE('vmix')
          WHERE(A2D < SPVAL)  A2D=A2D+vbaro

          CASE('thmix')
          WHERE(A2D < SPVAL)  A2D=A2D+THBASE

          CASE('dens')
          WHERE(A  < SPVAL)   A=A+ THBASE

          CASE('lthk')
          WHERE(A <SPVAL)    A=A*InvG

!          CASE('uvel')
          
!          DO k=1,kdm

!          WHERE(ubaro < SPVAL)  A(:,:,k)=A(:,:,k)+ubaro(:,:)
!          ENDDO

!          CASE('vvel')
          
!           DO k=1,kdm
            
!          WHERE(vbaro < SPVAL)  A(:,:,k)=A(:,:,k)+vbaro(:,:)

!          ENDDO

          CASE DEFAULT
          END SELECT


        CALL wrtNcdf(iidm,jjdm,kdm,ix1,jy1,SPVAL,dim3)
        END 
      
             
