module mod_sections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module contains routines and variables used to calculate section data.
! The main object of the module is to set up variables describing points along a
! section, specified via input files. The points along the sections are relative
! to a input model grid.
!
! Routines:
!  read_section_in   : parses inputfile (infile below) and sets the start and
!                      end points + names of sections. After reading we 
!                      have subsections - each with an associated section name.
!                      At this stage, the names subsections are not unique wrt 
!                      subsections. Routine also allocates necessary module
!                      variables
!
!  section_nodepoints: Calculates grid nodes along sections along with some other
!                      info - also sets up the directional flags used for
!                      transport calculations. Needs read_section_in to be
!                      called first. Sections are defined by "great circles"
!
!  sections_join     : Joins subsections into final sections. ( A section in the
!                      infile may have more than one segment. This step sets up
!                      the "final" section by joining the subsections)
!
!  save_section_nodes: Basically saves what was produced in section_nodepoints
!
!  read_section_nodes: Basically reads data saved by save_section_nodes
!
!  ncwrite_secdata   : Puts data along sections into netcdf files
!
!  get_node_data     : Input is a grid field, output is points along the
!                      specified section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! - Knut Liseter     2003 - First version based on code from Mats Bentsen
! - Knut Liseter Oct.2008 - Complete rewrite. Mainly for the nodepoints and transport 
!                           setup. Uses a much simpler logic for nodepoint extraction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




use mod_xc
implicit none

! Section definition file
character(len=*), parameter  :: infile2='sections.in' ! Name of input file
real,             parameter  :: cversion=1.2          ! Current version of input file
integer,          parameter  :: max_sdm=1000          ! Max points per section
integer,          parameter  :: max_sec=50            ! Max sections
logical                      :: verbose               ! Set this for debugging (sections.in)

! flags if section is on grid - default is in lonlat - geodesic measures
logical, save :: secOnGrid=.false.

! Section name
character(len=20),dimension(:), allocatable :: secname

! Section information
integer :: nsec                                       ! Number of sections
real,    dimension(:,:), allocatable,save :: seclat   ! First/2nd point on section
real,    dimension(:,:), allocatable,save :: seclon   ! First/2nd point on section
integer, dimension(:),   allocatable,save :: sdm      ! Number of nodes in section
integer, dimension(:,:), allocatable,save :: ndeipiv  ! i pivot point p-cell
integer, dimension(:,:), allocatable,save :: ndejpiv  ! j pivot point p-cell
real,    dimension(:,:), allocatable,save :: ndedist  ! eff. sec length in box
real   , dimension(:,:), allocatable,save :: ndelon   ! node longitue
real   , dimension(:,:), allocatable,save :: ndelat   ! node latitude
integer, dimension(:,:), allocatable,save :: ndeflagu ! node flag for u-velocity
integer, dimension(:,:), allocatable,save :: ndeflagv ! node flag for v-velocity


contains



   ! Routine basically reads sections.in and allocates/sets arrays seclon,
   ! seclat and secname. These are the start/end points of sections
   subroutine read_sections_in()
   implicit none
   character(len=20)  :: ch20
   character(len=80)  :: ch80
   character(len=20)  :: tmp_secname
   logical :: ex,first
   integer :: i,j,sec,ios
   real :: rl(6), tmp,fversion

   allocate(seclon (max_sec,2))
   allocate(seclat (max_sec,2))
   allocate(secname(max_sec))

   if (.not.allocated(seclon) .or. &
       .not.allocated(seclat) .or. &
       .not.allocated(secname))    &
   then
      print *,'Some arrays needed  are unallocated'
      print *, '(section_read)'
      call exit(1)
   end if

     
   ! Read the infile
   print *
   print *,'Reading sections'
   inquire(file=infile2,exist=ex)
   if (ex) then

      ! Header
      open(10,file=infile2)
      read (10,*) fversion, verbose

      ! 1.2 is still the official version. 1.3 allows one to specify
      ! if the start/end sections are on the model grid
      if (abs(fversion-1.3)<.001) then
         read (10,*) ch20
         if (trim(adjustl(ch20))=='grid') then
            secOnGrid=.true.
         end if
      end if

      ch20 ='###'
      do while (ch20(1:1)=='#')
         read (10,*) ch20
      end do
      backspace(10)
      
      j=1
      ios=0
      do while (j<=max_sec.and.ios==0)

         ! Read separation markers  (#)
         read (10,*,iostat=ios) ch80
         do while (ch80(1:1)=='#' .and. ios==0)
            read (10,*,iostat=ios) ch80
         end do

         if (ios==0) then 

            ! Read Section Name
            tmp_secname = trim(ch80)

            ! Read points up to next separation marker (#)
            read(10,'(a80)',iostat=ios) ch80
            first=.true.
            do while (ch80(1:1)/='#' .and. ios==0 .and. j<max_sec)
               !print *,ch80
               read(ch80,*,iostat=ios) rl(1),rl(2)

               if (first.and.ios==0) then 
                  first=.false.
                  seclon(j,1) =rl(1)
                  seclat(j,1) =rl(2)
               else if (ios==0) then
                  secname(j)=tmp_secname
                  seclon(j,2) =rl(1)
                  seclat(j,2) =rl(2)
                  secname(j+1)=tmp_secname
                  seclon(j+1,1) =rl(1)
                  seclat(j+1,1) =rl(2)
                  j=j+1
               end if
               read(10,'(a80)',iostat=ios) ch80
            end do

         end if
      end do
      nsec = j-1

      if (nsec > max_sec) then
         print *,'OBS max number of sections exceeded!'
      end if
   else 
      write (*,'(a)') 'No sections infile '//infile2//' exists. I will dump '
      write (*,'(a)') 'a sample infile between ------ markers below, from       '
      write (*,'(a)') 'which you can cut and paste.                             '
      fversion = cversion ; verbose=.false.
      seclon(1,1) = -10. ; seclon(1,2) = -  0. ; seclat(1,1)=60. ; seclat(1,2) = 60.
      secname(1) = '<Name here>'
      nsec=1
   end if


   ! Dump sections infile to screen
   if(.not.ex) &
   write(*, '(a)') '---------------------- START -- -----------------------------------'
   write(* ,'(f3.1,1xl1,15xa)') fversion, verbose,'# File version and verbosity flag'
   write(*,'(a)')  '#####'
   write(* ,'(a)') '##### Sections:  lon,lat pairs for start and end point.'

   write(*,'(a)') secname(1)
   write(* ,'(2f7.2)') seclon(1,1),seclat(1,1)
   write(* ,'(2f7.2)') seclon(1,2),seclat(1,2)
   do j=2,nsec
      
      if (secname(j-1)==secname(j)) then
         write(* ,'(2f7.2)') seclon(j,2),seclat(j,2)
      else 
         write(*,'(a)') '#'
         write(*,'(a)') secname(j)
         write(* ,'(2f7.2)') seclon(j,1),seclat(j,1)
         write(* ,'(2f7.2)') seclon(j,2),seclat(j,2)
      end if
         
   end do
   if(.not.ex) &
   write(*, '(a)') '----------------------- END ---------------------------------------'

   if (.not.ex) then
      print *
      call exit(1)
   end if
   if (fversion /=cversion .and.  abs( fversion-0.1 - cversion ) >1e-6 ) then
      print *,fversion,cversion
      print *,'Fileversion does not match current version -- update !!'
      print *,'Current version:',cversion,'  infile version:',fversion
      print*, '(secinit)'
      call exit(1)
   end if

   print *,'Sections read'
   print *
   end subroutine read_sections_in




   ! This troutine sets up the node points (p-cell) of the grid along the section
   subroutine section_nodepoints(lon,lat,idm,jdm,periodic)
   use mod_sphere_tools
   !use m_sort
   use m_handle_err
   use netcdf
   implicit none
   logical, intent(in) :: periodic
   integer, intent(in) :: idm,jdm
   real, dimension(idm,jdm), intent(in) :: lon,lat

   real   , dimension(idm,jdm) :: rmask
   logical, dimension(idm,jdm) :: mask
   integer, dimension(idm,jdm) :: flagu, flagv
   integer :: isec,i,j,ib,jb, ipnt,n
   integer, allocatable :: indx(:)
   real,dimension(3) :: nvec, rvec, rvec1, rvec2, cp1, cp2

   character(len=3) :: css
   integer :: idmid, jdmid, ncid, varid


   ! Sanity check
   if (.not.(allocated(seclon) .and. allocated(seclat) .and. &
             allocated(secname) )) then
      print *,'Error -- section arrays are unallocated'
      print *, '(transport_init)'
      call exit(1)
   end if

   allocate(ndeipiv (max_sdm,nsec))
   allocate(ndejpiv (max_sdm,nsec))
   allocate(ndelon  (max_sdm,nsec))
   allocate(ndelat  (max_sdm,nsec))
   allocate(ndeflagu(max_sdm,nsec))
   allocate(ndeflagv(max_sdm,nsec))
   allocate(ndedist (max_sdm,nsec))
   allocate(sdm     (max_sec))
   allocate(indx    (max_sdm))

   ! Loop over all sections
   do isec=1,nsec

      do i=1,max_sdm
         indx(i)=i
      end do

      ! Radius vectors for start and end points
      rvec1=geo2cart(seclon(isec,1),seclat(isec,1))
      rvec2=geo2cart(seclon(isec,2),seclat(isec,2))

      ! Normal vector of plane defined by cross product of
      ! 1) Vector from earth center to start of section
      ! 2) Vector from earth center to end   of section
      nvec=cross_product(rvec1,rvec2)

      ! Now go through grid and mark all points on one side of the sphere
      ! (i.e. all points whose (radius vecor X normal vector) is negative)
      do j=1,jdm
      do i=1,idm
         
         ! radius vector : earth center -> point
         rvec=geo2cart(lon(i,j),lat(i,j))

         ! dot product of radius vector with normal vector sets the mask
         mask(i,j)=dot_product(nvec,rvec)<0.

         ! init flags
         flagu(i,j)=0
         flagv(i,j)=0

      end do
      end do

      ! Now calculate the node points along the hemisphere line by 
      ! using a ``telescopic'' sum
      do j=1,jdm
      do i=1,idm

         if (periodic) then
            ib=mod(i,idm)+1
         else
            ib=min(i+1,idm)
         end if
         jb=min(j+1,jdm)

         if (mask(i,j)) then
            flagu(i,j) = flagu(i,j)+1
            flagv(i,j) = flagv(i,j)+1
            flagu(ib,j) = flagu(ib,j)-1
            flagv(i,jb) = flagv(i,jb)-1
         end if
      end do
      end do


      ! Remove points along boundary
      do i=1,idm
         flagu(i,  1)=0
         flagv(i,  1)=0
         flagu(i,jdm)=0
         flagv(i,jdm)=0
      end do

      if (.not.periodic)  then
      do j=1,jdm
         flagu(1,  j)=0
         flagv(1,  j)=0
         flagu(idm,j)=0
         flagv(idm,j)=0
      end do
      end if


      ! Finally reduce the number of points to those that are between
      ! start and end points of section
      do j=1,jdm
      do i=1,idm

         ! radius vector : earth center -> point
         rvec=geo2cart(lon(i,j),lat(i,j))

         ! Cross product start/end and this point
         cp1=cross_product(rvec1,rvec)
         cp2=cross_product(rvec,rvec2)

         ! These must have the same sign for rvec to be between 
         ! rvec1 and rvec2
         if (dot_product(cp1,cp2)<0.) then
            flagu(i,j)=0
            flagv(i,j)=0
         end if

      end do
      end do


      ! Number of node points
      sdm(isec)= count(flagu /=0 .or. flagv /=0 )

      ! security check on length (max_sdm)
      if (sdm(isec)>max_sdm) then
         print *,'Section length  of '//trim(secname(isec))//' exceeds max_sdm ',sdm(isec)
         call exit(1)
      else if (sdm(isec)==0) then
         print *,'Section length  of '//trim(secname(isec))//' is zero (on this grid). Remove section'
         call exit(1)
      end if

      ! Unsorted node points (for now)
      ipnt=1
      do j=1,jdm
      do i=1,idm
         if (flagu(i,j)/=0 .or. flagv(i,j)/=0) then
            ndelon (ipnt,isec)=lon(i,j)
            ndelat (ipnt,isec)=lat(i,j)
            ndeipiv(ipnt,isec)=i
            ndejpiv(ipnt,isec)=j
            ndeflagu(ipnt,isec)=flagu(i,j)
            ndeflagv(ipnt,isec)=flagv(i,j)
            ndedist(ipnt,isec)=spherdist(lon(i,j),lat(i,j), &
               seclon(isec,1),seclat(isec,1))
            ipnt =ipnt+1
         end if
      end do
      end do

      ! Sort by distance ...
      n=sdm(isec)
      call sort(n,ndedist(1:n,isec),indx(1:n))
      ndelon(1:n  ,isec)=ndelon(indx(1:n)  ,isec)
      ndelat(1:n  ,isec)=ndelat(indx(1:n)  ,isec)
      ndeipiv(1:n ,isec)=ndeipiv(indx(1:n) ,isec)
      ndejpiv(1:n ,isec)=ndejpiv(indx(1:n) ,isec)
      ndeflagu(1:n,isec)=ndeflagu(indx(1:n),isec)
      ndeflagv(1:n,isec)=ndeflagv(indx(1:n),isec)


      write(css,'(i3.3)') isec
      call handle_err(NF90_create('tst'//css//'.nc',NF90_CLOBBER,ncid))
      call handle_err(NF90_DEF_DIM(ncid,'idm',idm,idmid))
      call handle_err(NF90_DEF_DIM(ncid,'jdm',jdm,jdmid))
      call handle_err(NF90_DEF_VAR(ncid,'mask',NF90_Float,(/idmid,jdmid/),varid))
      call handle_err(NF90_ENDDEF(ncid))
      rmask=0. ; where(mask) rmask=1.
      call handle_err(NF90_PUT_VAR(ncid,varid,rmask))
      call handle_err(NF90_CLOSE(ncid))



   end do
   end subroutine section_nodepoints


   subroutine sections_join
   use mod_za
   use mod_grid
   use mod_sphere_tools
   implicit none

   character(len=20), dimension(max_sec) :: jsecname
   integer, dimension(max_sec) :: jsecindex,jsdm
   integer, dimension(max_sdm) :: j_ipiv, j_jpiv, j_flagu, j_flagv
   real   , dimension(max_sdm) :: jdist, jlon, jlat

   integer :: njsec, isec, counter,lw,up, ijsec, n


   ! Get Sections with common name (jsecindex)
   njsec=1
   jsecindex(1) = njsec
   do isec=2,nsec
      if (secname(isec)==secname(isec-1)) then
         jsecindex(isec) = njsec
      else
         njsec=njsec+1
         jsecindex(isec) = njsec
      end if
   end do

   ! Get size of  joint section
   jsdm=0
   do isec=1,nsec
      !if (verbose)print *,jsecindex(isec)
      jsdm(jsecindex(isec)) = jsdm(jsecindex(isec)) + sdm(isec)
      jsecname(jsecindex(isec))= secname(isec)
   end do


   ! Go through each joint section and extend it by its piecewise elements
   do ijsec=1,njsec

      ! security check on length (max_sdm)
      if (jsdm(ijsec)>max_sdm) then
         print *,'Joined section length exceeds max_sdm ',jsdm(ijsec)
         call exit(1)
      end if

      counter=1
      do isec=1,nsec
      if (jsecindex(isec)==ijsec) then

         lw=counter             ! lower index of section into joint section
         up=counter+sdm(isec)-1 ! upper index of section into joint section

         if (verbose) write(*,  &
         '("ijsec=",i3,"  lw=",i4,"  up=",i4,"  sdm(sec)=",i4," jsdm(ijsec)=",i4)') &
         ijsec,lw,up,sdm(isec),jsdm(ijsec)
         
         ! NB! distance is now cummulative ...
         jdist(lw:up) = ndedist(1:sdm(isec),isec)
         if (counter/=1)  jdist(lw:up) = jdist(lw:up) + jdist(lw-1)

         jlon(lw:up) = ndelon(1:sdm(isec),isec)
         jlat(lw:up) = ndelat(1:sdm(isec),isec)
         j_ipiv(lw:up) = ndeipiv(1:sdm(isec),isec)
         j_jpiv(lw:up) = ndejpiv(1:sdm(isec),isec)
         j_flagu(lw:up) = ndeflagu(1:sdm(isec),isec)
         j_flagv(lw:up) = ndeflagv(1:sdm(isec),isec)
         secname(ijsec)= secname(isec)

         counter=up+1

      end if
      end do

      ! Joint section is finished - put it back into section arrays at location 
      ! ijsec. This is safe since ijsec < jsecindex when jsecindex>ijsec
      sdm(ijsec) = counter-1
      n=sdm(ijsec)
      ndelon  (1:n,ijsec) = jlon(1:n)
      ndelat  (1:n,ijsec) = jlat(1:n)
      ndedist (1:n,ijsec) = jdist(1:n)
      ndeipiv (1:n,ijsec) = j_ipiv(1:n)
      ndejpiv (1:n,ijsec) = j_jpiv(1:n)
      ndeflagu(1:n,ijsec) = j_flagu(1:n)
      ndeflagv(1:n,ijsec) = j_flagv(1:n)
   end do

   ! Finally update number of sections
   nsec=njsec

   end subroutine sections_join


   subroutine read_section_nodes()
   implicit none

   character(len=80) :: fil_trans,fil_sec
   character(len= 3) :: css
   character(len= 1) :: c1
   integer :: ios,npt,maxnpt,nfiles,isec, ipnt,i,j
   logical :: lok,ex,ex_sec,ex_trans
   real :: rtmp


   ! Check for existence of files -- count them
   maxnpt=0
   do isec=1,max_sec
      write(css,'(i3.3)') isec
      fil_trans = "transport"//css//".dat"
      fil_sec   = "section"//css//".dat"
      inquire(exist=ex_trans,file=trim(fil_trans))
      inquire(exist=ex_sec  ,file=trim(fil_sec  ))
      if (ex_trans .and. ex_sec) then
         !print *,sec
         open(10,file=trim(fil_sec),action='read',status='old')
         read(10,'("%",i5.5)') npt
         close(10)
         maxnpt=max(npt,maxnpt) ! not really used
         nfiles=isec
      end if
   end do
   nsec=nfiles

   ! Allocate arrays -- sections
   ! TODO : check size fits
   allocate(ndeipiv (max_sdm,nsec))
   allocate(ndejpiv (max_sdm,nsec))
   allocate(ndeflagu(max_sdm,nsec))
   allocate(ndeflagv(max_sdm,nsec))
   allocate(ndelat  (max_sdm,nsec))
   allocate(ndelon  (max_sdm,nsec))
   allocate(ndedist (max_sdm,nsec))
   allocate(secname(nsec))
   allocate(sdm(nsec))

   ! Read sections
   do isec=1,nsec
      write(css,'(i3.3)') isec
      fil_sec   = "section"//css//".dat"
      open(10,file=trim(fil_sec),status='old',action='read')
      read(10,'("%",i5.5)') sdm(isec)

      do ipnt=1,sdm(isec)
         read(10,'(2i5,3e16.8,2x,2e16.8,x,a20)',iostat=ios) &
            ndeipiv(ipnt,isec),ndejpiv(ipnt,isec),          &
            ndelon(ipnt,isec),ndelat(ipnt,isec),ndedist(ipnt,isec), &
            rtmp,rtmp, secname(isec)
         if (ios/=0) then
            print *,'Error reading section ',isec
            stop '(mod_sections:read_section_nodes)'
         end if
      end do
      close(10)
   end do


   ! Read transports
   do isec=1,nsec
      write(css,'(i3.3)') isec
      fil_trans = "transport"//css//".dat"
      open(10,file=trim(fil_trans),status='old',action='read')
      do ipnt=1,sdm(isec)
         ! i,j are pivot points already read
         read(10,'(2i6,2i5)',iostat=ios) i, j, ndeflagu(ipnt,isec), ndeflagv(ipnt,isec)
         !write(6,'(2i6,2i5)',iostat=ios) i, j, ndeflagu(ipnt,isec), ndeflagv(ipnt,isec)
         if (ios/=0) then
            print *,'Error reading section ',isec
            stop '(mod_sections:read_section_nodes)'
         end if
      end do
      close(10)
   end do

   end subroutine read_section_nodes






   subroutine save_section_nodes()
   use mod_grid
   implicit none
   character(len=80) :: fil_int
   character(len= 3) :: css
   integer :: ipnt,isec,i,j


! ---------------------------------------------------------------------------------------
! The rest is for plotting points of sections to a tecplot file



   do isec=1,nsec

      ! For each section, dump the intersection points -- USED LATER
      write(css,"(i3.3)")isec
      fil_int = 'section'//css//".dat"
      open(10,file=trim(fil_int),status='replace')
      write(10,'(a1,i5.5,a)') '%',sdm(isec) , secname(isec)
      do ipnt=1,sdm(isec)
         write(10,'(2i5,3e16.8,2x,2e16.8,x,a20)') ndeipiv(ipnt,isec),ndejpiv(ipnt,isec), &
                                  ndelon(ipnt,isec),ndelat(ipnt,isec),ndedist(ipnt,isec), &
                                  0.,0., secname(isec)
                                  
      end do
      close(10)
      write(6,'(a)',advance='no' ) trim(fil_int)//'  '

      fil_int = 'transport'//css//".dat"
      open(10,file=trim(fil_int),status='replace')
      do ipnt=1,sdm(isec)
         write(10,'(2i6,2i5)') ndeipiv(ipnt,isec),ndejpiv(ipnt,isec),ndeflagu(ipnt,isec),ndeflagv(ipnt,isec)
      end do
      close(10)
      write(6,'(a)' ) trim(fil_int)
   end do


   end subroutine save_section_nodes




   subroutine ncwrite_secdata(cfld,field,k,kdm,var3d,vartime, &
      appendfile,fillvalue,comment)
   use m_handle_err
   use netcdf
   implicit none
   integer, intent(in) :: k,kdm
   character(len=*), intent(in) :: cfld
   real,             intent(in) :: field(idm,jdm)
   logical, intent(in), optional :: var3d
   real   , intent(in), optional :: vartime
   logical, intent(in), optional :: appendfile
   real   , intent(in), optional :: fillvalue
   character(len=*),  intent(in), optional :: comment

   character(len=80) :: ncfil
   character(len= 3) :: css
   logical :: ex,v3d, vtime, v2d, vapp
   integer :: ncid, varid, time_dimension_id, sect_dimension_id, &
      depth_dimension_id, vars_3d(3), vars_2d(2), vars_2d_fixed(1), &
      vars_3d_fixed(2), n, ierr, isec, time_var_id, dist_var_id, &
      kdm_var_id
   real :: tmpsec(max_sdm)
   integer, save :: rdimlen=-1

   v3d=.false.
   if (present(var3d)) v3d=var3d
   v2d=.not.v3d

   vtime=.false.
   if (present(vartime)) vtime=.true.

   vapp=.false.
   if (present(appendfile)) vapp=appendfile

   do isec=1,nsec 
      write(css,'(i3.3)') isec
      ncfil = 'section'//css//".nc"
      inquire(exist=ex,file=trim(ncfil))

      ! Base definition of file if it does not exist, or 
      ! if appendfile is switched off AND this is the first call
      if (.not. ex .or. (.not. vapp .and. rdimlen==-1)) then
         ! Define dimensions
         call handle_err(NF90_create(trim(ncfil),NF90_CLOBBER,ncid))
         call handle_err(NF90_DEF_DIM(ncid,'time'    ,NF90_UNLIMITED,time_dimension_id))
         call handle_err(NF90_DEF_VAR(ncid,'time'    ,NF90_FLOAT,time_dimension_id,time_var_id))
         call handle_err(NF90_DEF_DIM(ncid,'distance',sdm(isec)     ,sect_dimension_id))
         call handle_err(NF90_DEF_VAR(ncid,'distance',NF90_FLOAT,sect_dimension_id,dist_var_id))
         call handle_err(NF90_DEF_DIM(ncid,'depth'   ,kdm           ,depth_dimension_id))
         !call handle_err(NF90_DEF_VAR(ncid,'kdm',NF90_FLOAT,depth_dimension_id,kdm_var_id))
         vars_3d      =(/sect_dimension_id,depth_dimension_id,time_dimension_id/)
         vars_3d_fixed=(/sect_dimension_id,depth_dimension_id/)
         vars_2d      =(/sect_dimension_id,time_dimension_id/)
         vars_2d_fixed=(/sect_dimension_id/)
         if (isec==nsec) rdimlen=1 

         ! Some convenient global attributes
         call handle_err(nf90_put_att(ncid,NF90_GLOBAL, &
            'sectionname',trim(secname(isec))))
         call handle_err(nf90_enddef(ncid))
      else
         ! Get dimensions - we may need them for defining variables
         call handle_err(NF90_open(trim(ncfil),NF90_WRITE,ncid))
         call handle_err(nf90_inq_dimid(ncid, 'time'    , time_dimension_id))
         call handle_err(nf90_inq_varid(ncid, 'time'    , time_var_id))
         call handle_err(nf90_inq_dimid(ncid, 'distance', sect_dimension_id))
         call handle_err(nf90_inq_varid(ncid, 'distance', dist_var_id))
         call handle_err(nf90_inq_dimid(ncid, 'depth'   , depth_dimension_id))
         vars_3d      =(/sect_dimension_id,depth_dimension_id,time_dimension_id/)
         vars_3d_fixed=(/sect_dimension_id,depth_dimension_id/)
         vars_2d      =(/sect_dimension_id,time_dimension_id/)
         vars_2d_fixed=(/sect_dimension_id/)
         if (rdimlen == -1 ) then
            call handle_err(NF90_INQUIRE_DIMENSION(ncid,time_dimension_id,len=rdimlen))
            rdimlen=rdimlen+1
         end if
      end if

      ! Inquire and define (if necessary)  for current variable defined by cfld
      ierr=nf90_inq_varid(ncid,trim(cfld),varid)
      if (ierr/=NF90_NOERR) then

         call handle_err(nf90_redef(ncid))


         ! Define variable 
         if (v3d) then ! always w/time
            call handle_err(nf90_def_var(ncid,trim(cfld),NF90_Float,vars_3d,varid))
         elseif (v2d.and.vtime) then 
            call handle_err(nf90_def_var(ncid,trim(cfld),NF90_Float,vars_2d,varid))
         elseif (v2d) then
            call handle_err(nf90_def_var(ncid,trim(cfld),NF90_Float,vars_2d_fixed,varid))
         else
            print *,'Error - ncwritesectiondata'
            call exit(1)
         end if

         ! Put attributes
         if (present(fillvalue)) then
            call handle_err(nf90_put_att(ncid,varid, &
               '_FillValue',real(fillvalue,4)))
            call handle_err(nf90_put_att(ncid,varid, &
               'missing_value',real(fillvalue,4)))
         end if
         if (present(comment)) then
            call handle_err(nf90_put_att(ncid,varid, &
               'comment',trim(comment)))
         end if



         call handle_err(nf90_enddef(ncid))
      end if


      ! retrieve variables along section nodes
      call get_node_data(field,tmpsec,isec)

      ! Now put variable
      n=sdm(isec)

      if (vtime) call handle_err(nf90_put_var(ncid,time_var_id,vartime,start=(/rdimlen/)))
                 call handle_err(nf90_put_var(ncid,dist_var_id,ndedist(1:sdm(isec),isec),start=(/1,rdimlen/)))

      !print *,n,tmpsec
      if (v3d) then ! always w/time
         call handle_err(nf90_put_var(ncid,varid,tmpsec(1:n),start=(/1,k,rdimlen/)))
      elseif (v2d.and.vtime) then ! always w/time
         call handle_err(nf90_put_var(ncid,varid,tmpsec(1:n),start=(/1,rdimlen/)))
      elseif (v2d) then ! always w/time
         call handle_err(nf90_put_var(ncid,varid,tmpsec(1:n)))
      end if
      call handle_err(nf90_close(ncid))
   end do
   end subroutine
         




   subroutine get_node_data(field,line,isec)
   implicit none
   integer, intent(in) :: isec
   real, intent(in)  :: field(idm,jdm)
   real, intent(out) :: line(max_sdm)
   integer:: j,ipiv,jpiv

   do j=1,sdm(isec)
      ipiv=ndeipiv(j,isec)
      jpiv=ndejpiv(j,isec)
      line(j) = field(ipiv,jpiv)
   end do
   end subroutine


end module mod_sections

