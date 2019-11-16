program nestbat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program reads the global and local depths and pos files
! and creates a new local depths file which at all boundaries are
! consistent with the global file.
!
! NB: This version uses regional.grid and regional.depth files, 
!     in accordance with new setup in NERSC HYCOM 2.2
!
! The files read are 
!  global.grid.[ab]   (outer model grid)
!  global.depth.[ab]  (outer model bathymetry)
!  regional.grid.[ab]   (local model grid)
!  regional.depth.[ab]  (local model bathymetry)
!  grid.info from global grid
!
! The new depths file is saved on 
!  nesttbat.depth.[ab]
!
! To use: copy the required files to a new empty directory and
! execute nestbat....
! NB: This version is for hycom 2.2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use mod_parameters
   use mod_xc     , only: nxl => idm, nyl => jdm, xcspmd
   use mod_za     , only: zaiost, zaiopf, zaiowr, zaiocl
   use mod_grid   , only: ldepths=>depths, llon=>plon, llat=>plat, get_grid
   use mod_xc_global, only : nxg=>idm,nyg=>jdm, &
                             xcspmd_global=>xcspmd
   use mod_za_global, only : zaiost_global=>zaiost, &
                             zaiopf_global=>zaiopf, &
                             zaiocl_global=>zaiocl, &
                             zaiord_global=>zaiord, &
                             zaiosk_global=>zaiosk
   use mod_confmap
   implicit none

   real, allocatable :: depths2(:,:), depths(:,:), tmp(:,:)
   real, allocatable :: gdepths(:,:),glon(:,:),glat(:,:)
   real, allocatable :: testrel(:,:)
   integer, allocatable :: imask(:,:)
   character(len=7) tag7
   character(len=80) a80,matfile
   logical ex
   integer j,i,ipiv,jpiv,inest,l,i2,i0,j2
   real lon_n,lat_n,ba1,ba2,ba3,ba4,testrel2, xmin, xmax, xmin2, xmax2
   integer,parameter :: chunk = 120
   character(len=chunk) cchunk
   logical, allocatable :: flag(:,:),flag1(:,:)
   integer :: iter, numpoints, nlandneighbours
   logical :: lperiodic,inirange,samegrid
   integer :: igrace, jgrace, ipib, jpib
   integer :: idist, jdist
   real    :: rdist

print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'! This program reads the global and local depths and grid files       !'
print *,'! and creates a new local depths file which at all boundaries are     !'
print *,'! consistent with the global file (smoothed towards boundary). This is!'
print *,'! necessary for nesting...                                            !'
print *,'!                                                                     !'
print *,'! The files read are                                                  !'
print *,'!  regional.depth.[ab] - local  depths file                           !' 
print *,'!  regional.grid.[ab]  - local  grid   file                           !' 
print *,'!  global.depth.[ab]   - global depths file                           !' 
print *,'!  global.grid.[ab]    - global grid   file                           !' 
print *,'!  grid.info           - grid info for global grid (NB!)              !'
print *,'!                                                                     !'
print *,'! User input when running nestbat:                                    !'
print *,'!  None                                                               !'
print *,'!  nesting boundary is hardcoded to 20 grid cells                     !'
print *,'!                                                                     !' 
print *,'! The new depths file is saved in                                     !'
print *,'!  nestbat.depth.[ab]                                                 !'
print *,'!                                                                     !' 
print *,'! Diagnostic files are saved in                                       !'
print *,'!  nestrelmask.asc -- text file for visual check of the nesting zone  !'
print *,'!  nestbat.tec     -- tecplot file with old and new local depths ++   !'
print *,'!  nestbat.nc      -- netcdf  file with old and new local depths ++   !'
print *,'!                                                                     !' 
print *,'! To use: copy the required files to a new empty directory and        !'
print *,'! execute nestbat.... Note that global grid/depth files have to be    !'
print *,'! renamed from regional.* to global.*                                 !'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'


   samegrid=.false.
   if (samegrid) then
      print *,'Assuming the grids are the same'
   end if

   !inest is hardcoded
   inest=20

   ! NB: this will initialize the regional grid
   print '(a)','Init LOCAL grid'
   call xcspmd()
   call zaiost()
   call get_grid()
   allocate(depths(nxl,nyl),tmp(nxl,nyl),depths2(nxl,nyl))
   allocate(testrel(nxl,nyl))
   allocate(flag(nxl,nyl),flag1(nxl,nyl))

   ! NB: this will initialize the global grid - mod_za_global is set up to read
   !     global.grid.[ab] for initialization
   print *
   print '(a)','Init GLOBAL grid'
   call xcspmd_global()
   call zaiost_global()
   allocate(gdepths(nxg,nyg),glon(nxg,nyg),glat(nxg,nyg),imask(nxg,nyg))
   call zaiopf_global('global.grid.a','old',711)
   call zaiord_global(glon,imask,.false.,xmin,xmax,711)
   call zaiord_global(glat,imask,.false.,xmin,xmax,711)
   call zaiocl_global(711)
   call zaiopf_global('global.depth.a','old',711)
   call zaiord_global(gdepths,imask,.false.,xmin,xmax,711)
   call zaiocl_global(711)
   where (gdepths > 0.5* huge) gdepths=0.
   deallocate(imask)
   !print *,minval(gdepths),maxval(gdepths)
   !print *,minval(glon),maxval(glon)
   !print *,minval(glat),maxval(glat)
   print *

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
      enddo
      enddo
   else
      if (nxg/=nxl .or. nyg/=nyl) then
         print *,'Samegrid is true but dimensions are different !'
         call exit(1)
      end if
      tmp=gdepths
   end if


   ! New (and more intuitive approach). 
   do j=2,nyl-1
   do i=2,nxl-1
      idist=min(i,nxl-i)
      jdist=min(j,nyl-j)

      ! Distance to nearest open boundary
      rdist=inest+100
      do j2=max(2,j-inest),min(nyl-1,j+inest)
      do i2=max(2,i-inest),min(nxl-1,i+inest)
         if      ((j2==2 .or. j2==nyl-1) .and. ldepths(i2,j2)>1.) then
            rdist=min(rdist,sqrt(real(i-i2)**2 + real(j-j2)**2))
         else if ((i2==2 .or. i2==nxl-1) .and. ldepths(i2,j2)>1.) then
            rdist=min(rdist,sqrt(real(i-i2)**2 + real(j-j2)**2))
         end if
      end do
      end do

      ! If distance to nearest open boundary is < inest it can be smoothed
      if (rdist<real(inest)+1e-4 .and.ldepths(i,j)>1.) then
         flag(i,j)=.true.
         testrel(i,j)=rdist/float(inest-1)
         testrel(i,j)=(testrel(i,j)**2.-1.)**4.           !(i=1)=1,(i=20)=0
      end if
   end do
   end do

   ! Here the actual smoothing to the outer grid takes place
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
      print '(a,i4,a,i7,a)','Depth check iteration ',iter,' has ', numpoints ,' modified points'
      iter=iter+1
   end do
   print *




   ! dumping new depths file - note that new setup requires "huge" values
   ! for land points
   allocate(imask(nxl,nyl))
   where(depths<.1) depths=huge
   xmax=maxval(depths,mask=depths<0.5*huge)
   xmin=minval(depths,mask=depths<0.5*huge)
   call zaiopf('nestbat.depth.a','replace',711)
   call zaiowr(depths, imask,.false., xmin2,xmax2,711,.false.)
   call zaiocl(711)
   !print *,xmin,xmax
   !print *,xmin2,xmax2,2.0**99
   open (10,file='nestbat.depth.b', status='replace')
   write(10,'(a)') 'Depth grid modified by nestbat '
   write(10,'(a)') ''
   write(10,'(a)') ''
   write(10,'(a)') ''
   write(10,'(a)') ''
   write(10,'(a,2f10.3)') 'min,max depth = ',xmin,xmax
   close(10)

   ! NB - remove huge values again for diagnostics 
   where(depths>0.5*huge) depths=0.
   call tecdump(nxl,nyl,nxg,nyg,depths,tmp,ldepths,gdepths,glon,glat,llon,llat)

   !print'(i4,140l1)',(j,(flag(i,j), i=1,nxl), j=nyl,50,-1)
   open(10,file='nestrelmask.asc')
   write(10,'(a,2i6)') 'Nesting mask for grid of size ', nxl, nyl
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


   print '(a)','nestrel mask in nestrelmask.asc (visual ascii file)'
   print '(a)','New bathymetry dumped to nestbat.depth.[ab]'
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!KAL - now also dumps to netcdf file
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

   open(10,file='nestbat.tec',status='unknown')
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

! NEstbat depths - original depths
      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'difference',NF90_Float,(/nxlid,nylid/),var_id)
      ierr=NF90_PUT_ATT(ncid,var_id,'comment','nestbat depths minus original depths')
      ierr=NF90_PUT_ATT(ncid,var_id,'_FillValue',real(-1e14,kind=4))
      ierr=NF90_PUT_ATT(ncid,var_id,'missing_value',real(-1e14,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpl=depths-ldepths; where (depths<.1 .and. ldepths<.1) tmpl=-1e14
      ierr=NF90_PUT_VAR(ncid,var_id,tmpl)
      ierr=NF90_CLOSE(ncid)

      print '(a)','Diagnostics in nestbat.nc and nestbat.tec'
end subroutine tecdump

end program
