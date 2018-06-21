program nestbat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program reads the global and local depths and newpos files
! and creates a new local depths file which at all boundaries are
! consistent with the global file.

! The files read are 
!  gdepths???x???.uf
!  gnewpos.uf
!  ldepths???x???.uf
!  lnewpos.uf
!  grid.info from global grid

! The new depths file is saved on 
!  ndepths???x???.uf   
!
! To use: copy the required files to a new empty directory and
! execute nestbat....
! NB: This version is for hycom 2.1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use mod_confmap
   implicit none

   integer nxl,nyl,nxg,nyg
   real, allocatable :: ldepths(:,:),llon(:,:),llat(:,:),tmp(:,:), &
      depths2(:,:), depths(:,:)
   real, allocatable :: gdepths(:,:),glon(:,:),glat(:,:)
   real, allocatable :: testrel(:,:)
   real*8, allocatable :: io1(:,:), io2(:,:)
   character(len=7) tag7
   character(len=80) a80,matfile
   logical ex
   integer j,i,ipiv,jpiv,inest,l,i2,i0,j2
   real lon_n,lat_n,ba1,ba2,ba3,ba4,testrel2
   integer,parameter :: chunk = 120
   character(len=chunk) cchunk
   logical, allocatable :: flag(:,:),flag1(:,:)

   integer :: iter, numpoints, nlandneighbours

   logical :: lperiodic,inirange,samegrid
   integer :: igrace, jgrace, ipib, jpib

   print *,'Type help for help, type enter to continue'
   read(*,'(a)') a80
   if (trim(a80)=='help') then
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'! This program reads the global and local depths and newpos files     !'
print *,'! and creates a new local depths file which at all boundaries are     !'
print *,'! consistent with the global file (smoothed towards boundary). This is!'
print *,'! necessary for nesting...                                            !'
print *,'!                                                                     !'
print *,'! The files read are                                                  !'
print *,'!  gdepths???x???.uf   - global depths file                           !' 
print *,'!  gnewpos.uf          - global newpos file                           !'
print *,'!  ldepths???x???.uf   - local  depths file                           !'
print *,'!  lnewpos.uf          - local  newpos file                           !'
print *,'!  grid.info           - grid info for global grid                    !'
print *,'!                                                                     !'
print *,'! User input when running nestbat:                                    !'
print *,'!  Grid dimensions for global grid                                    !'
print *,'!  Grid dimensions for local  grid                                    !'
print *,'!  Width of boundary zone (transition zone from global to local grid) !'
print *,'!                                                                     !' 
print *,'! The new depths file is saved in                                     !'
print *,'!  ndepths???x???.uf                                                  !'
print *,'!                                                                     !' 
print *,'! Diagnostic files are saved in                                       !'
print *,'!  nestrelmask.asc -- text file for visual check of the nesting zone  !'
print *,'!  tecnestbat.tec  -- tecplot file with old and new local depths      !'
print *,'!                                                                     !' 
print *,'! To use: copy the required files to a new empty directory and        !'
print *,'! execute nestbat....                                                 !'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop '(nestbat)'
   end if

   write(*,'(a)',advance='no')'Chose grid-dimensions for global grid (nx ny): '
   read(*,*)nxg,nyg
   write(*,'(a)',advance='no')'Chose grid-dimensions for local  grid (nx ny): '
   read(*,*)nxl,nyl
   write(*,'(a)',advance='no')'Width of boundary zone: '
   read(*,*)inest

   write(*,*)
   write(*,'(a)')'The following should be set to true for Mercator nesting, false otherwise'
   write(*,'(a)',advance='no')'T if grids are identical, F if not '
   read(*,*)samegrid

   if (samegrid) then
      print *,'Assuming the grids are the same'
   end if

   allocate(ldepths(nxl,nyl),llon(nxl,nyl),llat(nxl,nyl))
   allocate(gdepths(nxg,nyg),glon(nxg,nyg),glat(nxg,nyg))
   allocate(depths(nxl,nyl),tmp(nxl,nyl),depths2(nxl,nyl))
   allocate(testrel(nxl,nyl))
   allocate(flag(nxl,nyl),flag1(nxl,nyl))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reading local files
   allocate(io1(nxl,nyl))
   allocate(io2(nxl,nyl))
   write(tag7,'(i3.3,a,i3.3)')nxl,'x',nyl
   inquire(file='ldepths'//tag7//'.uf',exist=ex)
   if (.not.ex) stop 'ldepths???x???.uf file does not exist'
   open (unit=10,file='ldepths'//tag7//'.uf',status='old',form='unformatted')
   !read(10)ldepths
   read(10) io1
   close(10)
   ldepths=io1
   where (ldepths > 1e25) ldepths=0.

   inquire(file='lnewpos.uf',exist=ex)
   if (.not.ex) stop 'lnewpos.uf file does not exist'
   open(10,file='lnewpos.uf',form='unformatted',status='old')
   !read(10)llat,llon
   read(10)io1,io2
   close(10)
   llat=io1
   llon=io2
   print *,minval(io1), maxval(io1)
   print *,minval(io2), maxval(io2)
   deallocate(io1,io2)

! reading global files
   allocate(io1(nxg,nyg))
   allocate(io2(nxg,nyg))
   write(tag7,'(i3.3,a,i3.3)')nxg,'x',nyg
   inquire(file='gdepths'//tag7//'.uf',exist=ex)
   if (.not.ex) stop 'gdepths???x???.uf file does not exist'
   open (unit=10,file='gdepths'//tag7//'.uf',status='old',form='unformatted')
   !read(10)gdepths
   read(10)io1
   close(10)
   gdepths=io1
   where (gdepths > 1e25) gdepths=0.

   inquire(file='gnewpos.uf',exist=ex)
   if (.not.ex) stop 'gnewpos.uf file does not exist'
   open(10,file='gnewpos.uf',form='unformatted',status='old')
   !read(10)glat,glon
   read(10)io1,io2
   close(10)
   glat=io1
   glon=io2
   !
   print *,minval(io1), maxval(io1)
   deallocate(io1,io2)



!  Check if global grid is periodic - shown by i-boundary values

   lperiodic=.false.
   if (any(gdepths(1,:) > 0.1) .and. any(gdepths(nxg,:) > 0.1) ) then
      print *,'Periodic global grid '
      lperiodic=.true.
   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   tmp=ldepths

   ! Interpolate global grid to local grid
   if (.not. samegrid) then
   call initconfmap(nxg,nyg)
   do j=1,nyl
   do i=1,nxl
      call oldtonew(llat(i,j),llon(i,j),lat_n,lon_n)
      call pivotp(lon_n,lat_n,ipiv,jpiv)

      !print *,'bef:',i,j,ipiv,jpiv

      if (lperiodic) then
         ipib=mod(ipiv,nxg)+1
         inirange=.true.
      else
         ipib=ipiv+1
         inirange=ipiv>=1 .and. ipiv < nxg
      end if
      jpib=jpiv+1

      ! grace for i
      igrace=min(nxg-ipiv,ipiv-1) ! negative when ipiv < 1 or ipiv > nxl

      ! grace for j
      jgrace=min(nyg-jpiv,jpiv-1) ! negative when jpiv < 1 or jpiv > nxl

      !print *,'bef2:',i,j,ipiv,jpiv,ipib,jpib


      if (inirange .and. jpiv >=1 .and. jpiv < nyg ) then

         !call bilincoeff(glon,glat,nxg,nyg,llon(i,j),llat(i,j),ipiv,jpiv,ba1,ba2,ba3,ba4,lperiodic)
         call bilincoeff(glon,glat,nxg,nyg,llon(i,j),llat(i,j),ipiv,jpiv,ba1,ba2,ba3,ba4)

      ! Outside of global model grid -- but inner model has zero depths
      else if (ldepths(i,j)<0.1) then
         ! Simply set to land (weights zero)
         ipiv=max(1,min(ipiv,nxg))
         jpiv=max(1,min(jpiv,nyg))
         ipib=ipiv
         jpib=jpiv
         ba1=0. ; ba2=0.; ba3=0. ; ba4=0.

      else if (igrace>-50 .and. jgrace > -50 ) then
         print *,'Warning: Pivot point outside model domain - saved by grace'
         print *,i,j,ipiv,jpiv

         ! Simply set to land (weights zero)
         ipiv=max(1,min(ipiv,nxg))
         jpiv=max(1,min(jpiv,nyg))
         ipib=ipiv
         jpib=jpiv
         ba1=0. ; ba2=0.; ba3=0. ; ba4=0.
      else
         print *,'Warning: Pivot point outside model domain - not saved by grace'
         print *,i,j,ipiv,jpiv
         stop
      end if
      !print *,'aft:',i,j,ipiv,jpiv,ipib,jpib

      tmp(i,j)=ba1*gdepths(ipiv,jpiv) + ba2*gdepths(ipib,jpiv)&
              +ba3*gdepths(ipib,jpib) + ba4*gdepths(ipiv,jpib)
      ! Grace 
   enddo
   enddo
   else

      if (nxg/=nxl .or. nyg/=nyl) then
         print *,'Samegrid is true but dimensions are different !'
         call exit(1)
      end if
      tmp=gdepths
   end if



! flag- the logical array to check if the nestmodel domain connects directly to the open
! boundary: to make a smooth transition only in those points.

   flag=.false.
   flag1=.false.
 
   testrel=0.
   testrel(1,:)=1.
   testrel(nxl,:)=1.
   testrel(:,1)=1.
   testrel(:,nyl)=1.

   forall (j=2:nyl-1, ldepths(2,j) > 1.0)
    flag(2,j)=.true.
    testrel(2,j)=((1./float(inest-1))**2.-1.)**4.
   end forall
   forall (j=2:nyl-1, ldepths(nxl-1,j) > 1.0)
     flag(nxl-1,j)=.true.
     testrel(nxl-1,j)=((1./float(inest-1))**2.-1.)**4.
   end forall


   do i=3,20                     !N of itteration
    do i2=3,i                    !points checked for the itteration
     do j=3,nyl-2
      if(.not.flag(i2,j).and.ldepths(i2,j) > 1.0) then
   if(flag(i2-1,j-1).or.flag(i2-1,j).or.flag(i2-1,j+1).or. &
   & flag(i2,j-1).or.flag(i2,j+1).or.                      &
   & flag(i2+1,j-1).or.flag(i2,j).or.flag(i2+1,j+1))  then
         flag1(i2,j)=.true.
         testrel(i2,j)=float(i-1)/float(inest-1)
         testrel(i2,j)=(testrel(i2,j)**2.-1.)**4.                !(i=1)=1,(i=20)=0 
   endif
       endif
     enddo
    enddo
     forall (i2=3:i,j=3:nyl-2, flag1(i2,j)) flag(i2,j)=.true.     
   enddo

   do i=nxl-2,nxl-19,-1
    do i2=nxl-2,i,-1
     do j=3,nyl-2
      if(.not.flag(i2,j).and.ldepths(i2,j) > 1.0) then
   if(flag(i2-1,j-1).or.flag(i2-1,j).or.flag(i2-1,j+1).or. &
   & flag(i2,j-1).or.flag(i2,j+1).or.                      &
   & flag(i2+1,j-1).or.flag(i2,j).or.flag(i2+1,j+1))  then
         flag1(i2,j)=.true.
         testrel(i2,j)=float(nxl-i)/float(inest-1)
         testrel(i2,j)=(testrel(i2,j)**2.-1.)**4.           !(i=1)=1,(i=20)=0
   endif
       endif
     enddo
    enddo
     forall (i2=nxl-2:i:-1,j=3:nyl-2, flag1(i2,j)) flag(i2,j)=.true. 
   enddo


     flag=.false.                        !to include the corners in the analysis
     flag1=.false.

   forall (i=2:nxl-1, ldepths(i,2) > 1.0)
    flag(i,2)=.true.
    testrel(i,2)=((1./float(inest-1))**2.-1.)**4.
   end forall
   forall (i=2:nxl-1, ldepths(i,nyl-1) > 1.0)
     flag(i,nyl-1)=.true.
     testrel(i,nyl-1)=((1./float(inest-1))**2.-1.)**4.
   end forall

   do j=3,20
    do i2=3,j
     do i=3,nxl-2
      if(.not.flag(i,i2).and.ldepths(i,i2) > 1.0) then
   if(flag(i-1,i2-1).or.flag(i,i2-1).or.flag(i+1,i2-1).or. &
   & flag(i-1,i2).or.flag(i+1,i2).or.                      &
   & flag(i-1,i2+1).or.flag(i,i2+1).or.flag(i+1,i2+1))  then
         flag1(i,i2)=.true.
         testrel2=float(j-1)/float(inest-1) 
         testrel(i,i2)=max((testrel2**2.-1.)**4.,testrel(i,i2))     !(i=1)=1,(i=20)=0
   endif
       endif
     enddo
    enddo 
     forall (i=3:nxl-2,i2=3:j, flag1(i,i2)) flag(i,i2)=.true.
   enddo

   do j=nyl-2,nyl-19,-1
    do i2=nyl-2,j,-1
     do i=3,nxl-2
      if(.not.flag(i,i2).and.ldepths(i,i2) > 1.0) then
   if(flag(i-1,i2-1).or.flag(i,i2-1).or.flag(i+1,i2-1).or. &
   & flag(i-1,i2).or.flag(i+1,i2).or.                      &
   & flag(i-1,i2+1).or.flag(i,i2+1).or.flag(i+1,i2+1))  then
         flag1(i,i2)=.true.
         testrel2=float(nyl-j)/float(inest-1)
         testrel(i,i2)=max((testrel2**2.-1.)**4.,testrel(i,i2))     !(i=1)=1,(i=20)=0
   endif
       endif
     enddo
    enddo
     forall (i=3:nxl-2,i2=nyl-2:j:-1, flag1(i,i2)) flag(i,i2)=.true.
   enddo
 

   open(10,file='nestbnd.dat')
      do j=nyl,1,-1
         write(10,'(i4,40f5.2)') j,(testrel(i,j), i=101,140)
      enddo
   close(10)


   ! Here the actual smoothing to the uter grid takes place
!       where(ldepths<=4.)ldepths=0.0
!       where(ldepths>0.0.and.ldepths<=30.)ldepths=20.1
   depths=0.
   do j=2,nyl-1
   do i=2,nxl-1
      depths(i,j)=(1.0-testrel(i,j))*ldepths(i,j)+testrel(i,j)*tmp(i,j)
   enddo
   enddo
   where(depths<=4.)depths=0.0
   where(depths>0.0.and.depths<=30.)depths=20.1

! Check the depths matrix for points with +three neighbours
   depths2=depths
   iter=1
   numpoints=1
   do while (numpoints>0)
      numpoints=0
      do j=2,nyl-1
      do i=2,nxl-1
      if (depths(i,j)>0.1) then
         nlandneighbours=0
         if (depths(i+1,j  )<0.1) nlandneighbours = nlandneighbours + 1
         if (depths(i-1,j  )<0.1) nlandneighbours = nlandneighbours + 1
         if (depths(i  ,j+1)<0.1) nlandneighbours = nlandneighbours + 1
         if (depths(i  ,j-1)<0.1) nlandneighbours = nlandneighbours + 1

         if (nlandneighbours>=3) then
            depths2(i,j)=0.
            numpoints=numpoints+1
         end if
      end if
      enddo
      enddo
      depths=depths2
      print '(a,i4,a,i7,a)','Iteration ',iter,' has ', numpoints ,' modified points'
      iter=iter+1
   end do




! dumping new depths file
   allocate(io1(nxl,nyl))
   io1=depths
   write(tag7,'(i3.3,a,i3.3)')nxl,'x',nyl
   open (unit=10,file='ndepths'//tag7//'.uf',status='unknown',form='unformatted')
   !write(10)depths
   write(10)io1
   close(10)
   deallocate(io1)


   call tecdump(nxl,nyl,nxg,nyg,depths,tmp,ldepths,gdepths,glon,glat,llon,llat)

   !print'(i4,140l1)',(j,(flag(i,j), i=1,nxl), j=nyl,50,-1)
   open(10,file='nestrelmask.asc')
   write(10,*) 'Nesting mask for grid '//tag7
   write(10,*) 'Shows how much nesting pulls local bathymetry towards global file'
   write(10,*) 'Legend: *=land , .=unmodified ocean , [0-9] nestrel strength'
   do i=1,nxl,chunk
   cchunk=''
      write(10,*) 'i chunk ',i,' -- > ',min(nxl,i+chunk)
      write(10,*)
      do j=nyl,1,-1
         do i2=i,min(nxl,i+chunk)
            if (depths(i2,j)<5) then
               cchunk(i2-i+1:i2-i+2)= '*'
            elseif (testrel(i2,j)>.001) then
               write(cchunk(i2-i+1:i2-i+2),'(i1.1)')  &
                  floor(9.999*(testrel(i2,j))/maxval(testrel))
                  !floor(9*(1.-testrel(i2,j))/maxval(1.-testrel))
            else
               cchunk(i2-i+1:i2-i+2)='.'
            end if
         end do
         !print *,cchunk
         write(10,'(a)') cchunk
      end do
      write(10,*)
   end do
   close(10)


   print *,'nestrel mask in nestrelmask.asc (visual ascii file)'
   print *,'new bathymetry in ndepths'//tag7//'.uf'
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!KAL - also dumps to netcdf file
subroutine tecdump(nxl,nyl,nxg,nyg,depths,tmp,ldepths,gdepths,glon,glat,llon,llat)
   use netcdf
   implicit none
   integer, intent(in) :: nxl,nyl,nxg,nyg
   real, intent(in) :: depths(nxl,nyl)
   real, intent(in) :: ldepths(nxl,nyl)
   real, intent(in) :: tmp(nxl,nyl)
   real, intent(in) :: llon(nxl,nyl)
   real, intent(in) :: llat(nxl,nyl)
   real, intent(in) :: gdepths(nxg,nyg)
   real, intent(in) :: glon(nxg,nyg)
   real, intent(in) :: glat(nxg,nyg)
   integer i,j,k

   character*2 tag
   integer :: ierr
   integer :: ncid, var_id, nxlid,nylid,nxgid,nygid
   real :: tmpg(nxg,nyg)
   real :: tmpl(nxl,nyl)

   open(10,file='tecnestbat.tec',status='unknown')
      write(10,*)'TITLE = "NESTBAT"'
      write(10,*)'VARIABLES = "i" "j" "lon" "lat" "depths"'

!Global zone
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK T="Global", I=',nxg,', J=',nyg,', K=1'
      write(10,'(30I4)')((i,i=1,nxg),j=1,nyg)
      write(10,'(30I4)')((j,i=1,nxg),j=1,nyg)
      write(10,900)glon(1:nxg,1:nyg)
      write(10,900)glat(1:nxg,1:nyg)
      write(10,900)gdepths(1:nxg,1:nyg)


! Old local zone
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK T="Local", I=',nxl,', J=',nyl,', K=1'
      write(10,'(30I4)')((i,i=1,nxl),j=1,nyl)
      write(10,'(30I4)')((j,i=1,nxl),j=1,nyl)
      write(10,900)llon(1:nxl,1:nyl)
      write(10,900)llat(1:nxl,1:nyl)
      write(10,900)ldepths(1:nxl,1:nyl)

! tmp local zone
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK T="LocalG", I=',nxl,', J=',nyl,', K=1'
      write(10,'(30I4)')((i,i=1,nxl),j=1,nyl)
      write(10,'(30I4)')((j,i=1,nxl),j=1,nyl)
      write(10,900)llon(1:nxl,1:nyl)
      write(10,900)llat(1:nxl,1:nyl)
      write(10,900)tmp(1:nxl,1:nyl)

! New local zone
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK T="New", I=',nxl,', J=',nyl,', K=1'
      write(10,'(30I4)')((i,i=1,nxl),j=1,nyl)
      write(10,'(30I4)')((j,i=1,nxl),j=1,nyl)
      write(10,900)llon(1:nxl,1:nyl)
      write(10,900)llat(1:nxl,1:nyl)
      write(10,900)depths(1:nxl,1:nyl)

   close(10)
 900 format(10(1x,e12.5))

      if (NF90_CREATE('nestbat.nc',NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(obsstats)'
      end if
      ierr=NF90_DEF_DIM(ncid,'nxl',nxl,nxlid)
      ierr=NF90_DEF_DIM(ncid,'nyl',nyl,nylid)
      ierr=NF90_DEF_DIM(ncid,'nxg',nxg,nxgid)
      ierr=NF90_DEF_DIM(ncid,'nyg',nyg,nygid)

      ierr=NF90_DEF_VAR(ncid,'glon',NF90_Float,(/nxgid,nygid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,glon)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'glat',NF90_Float,(/nxgid,nygid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,glat)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'llon',NF90_Float,(/nxlid,nylid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,llon)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'llat',NF90_Float,(/nxlid,nylid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,llat)

! Global depths
      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'gdepths',NF90_Float,(/nxgid,nygid/),var_id)
      ierr=NF90_PUT_ATT(ncid,var_id,'_FillValue',real(-1e14,kind=4))
      ierr=NF90_PUT_ATT(ncid,var_id,'missing_value',real(-1e14,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpg=gdepths; where (tmpg<.1) tmpg=-1e14
      ierr=NF90_PUT_VAR(ncid,var_id,tmpg)

! Original depths
      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'ldepths',NF90_Float,(/nxlid,nylid/),var_id)
      ierr=NF90_PUT_ATT(ncid,var_id,'_FillValue',real(-1e14,kind=4))
      ierr=NF90_PUT_ATT(ncid,var_id,'missing_value',real(-1e14,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpl=ldepths; where (tmpl<.1) tmpl=-1e14
      ierr=NF90_PUT_VAR(ncid,var_id,tmpl)

! NEstbat depths
      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'ndepths',NF90_Float,(/nxlid,nylid/),var_id)
      ierr=NF90_PUT_ATT(ncid,var_id,'_FillValue',real(-1e14,kind=4))
      ierr=NF90_PUT_ATT(ncid,var_id,'missing_value',real(-1e14,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpl=depths; where (tmpl<.1) tmpl=-1e14
      ierr=NF90_PUT_VAR(ncid,var_id,tmpl)

      ierr=NF90_CLOSE(ncid)
end subroutine tecdump

end program
