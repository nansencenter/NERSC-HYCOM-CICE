program curviint
   use mod_xc_local , only : lnx=>idm,lny=>jdm
   use mod_xc_global, only : gnx=>idm,gny=>jdm
   use mod_za_local , only : zaiost_local =>zaiost, &
                             zaiopf_local =>zaiopf, &
                             zaiocl_local =>zaiocl, &
                             zaiowr_local =>zaiowr, &
                             zaiord_local =>zaiord
   use mod_za_global, only : zaiost_global=>zaiost, &
                             zaiopf_global=>zaiopf, &
                             zaiocl_global=>zaiocl, &
                             zaiord_global=>zaiord, &
                             zaiosk_global=>zaiosk
   use mod_confmap
   implicit none

   ! Local model - 
   real, dimension(:,:), allocatable ::    &
      llat, llon, ldepths, sumdeep
   integer, dimension(:,:), allocatable ::    pass
   logical, dimension(:,:), allocatable ::    &
    lmsk, donep
   real, dimension(:,:), allocatable ::   &
      lfld,ltemp,sumdepth,lpbot,gpbot,    &
      gpsikk,lpsikk,gthkk,lthkk,loldfld
   integer, dimension(:,:), allocatable :: &
      ipiv,jpiv,ldummy
   integer, dimension(:,:,:), allocatable :: &
      index_filled
   real   , dimension(:,:,:), allocatable :: &
      a
   logical, dimension(:,:), allocatable :: &
      ass

   real*8 , dimension(:,:), allocatable :: &
      gficem,ghicem,ghsnwm,gticem,gtsrfm
   real*8 , dimension(:,:), allocatable :: &
      lficem,lhicem,lhsnwm,lticem,ltsrfm

   ! Global model -
   real, dimension(:,:), allocatable ::   &
      tmp, glat, glon, gdepths
   logical, allocatable, dimension(:,:) :: gmsk
   real,    allocatable, dimension(:,:) :: &
      gfld,gtemp
   integer, dimension(:,:), allocatable :: &
      gdummy

   ! IO arrays
   real*8, dimension(:,:), allocatable :: io1, io2

   character(len=100) grestart
   character(len=80 ) a80,licefname,gicefname
   character(len=11 ) tag7
   character(len=8) :: cfld,oldcfld
   character(len=6) :: cvarin


   integer i,j, i2, j2,lrecl,grecl
   integer :: nirec,nhrec,irec,ios

   real, parameter ::  onem=9806.
   real :: bxmin, bxmax, axmin, axmax
   integer :: itime, lcoord,find,istep,oldlcoord,oldistep
   integer,parameter :: gnop=333, lnop=444, lnop2=555
   real time0
   logical ex
   logical lperiodic
   real :: rvar
   integer :: kdm
   real, dimension(:,:,:,:), allocatable :: dp
   real :: ba1,ba2,ba3,ba4,basum
   integer :: i2pb, j2pb


   character(len=2) tag2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading local depths and newpos file
   write(*,*)' Do you need help at this stage?'
   write(*,*) 'If so, type help here... Otherwise press enter'
   read(*,'(a)') a80
   if (trim(a80)=='help') then
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'! Program for initialization of a new grid using model output from a'
print *,'! another "global" model simulation (not really global but should contain'
print *,'! the model domain of the local model grid you interpolate to). '
print *,'! Note the following:'
print *,'!  0. Always do this in a separate directory.   Copy all files needed to'
print *,'!     this directory to awoid any unwanted effects.'
print *,'!  1. The global model grid is assumed to be orthogonal curvilinear based'
print *,'!     on the CONFGRID program.'
print *,'!  2. The global model restart files are used (ice is optional).'
print *,'!  3. The grid.info file used to generate the global grid must be available'
print *,'!     in the run directory.   '
print *,'!  4. The global depths file must be available. It is named gdepths???x???.uf'
print *,'!  5. The global newpos file must be available. It is named gnewpos.uf'
print *,'!  6. The local newpos.uf must be available.  It is named lnewpos.uf.'
print *,'!  7. The local depths.uf must be available. It is named ldepths???x???.uf'
print *,'!  8. global and local grid dimensions are specified by you in the program.'
print *,'!  9. The input global restart files to be used are specified by you'
print *,'!     in the program. The program uses the new HYCOM restart files (.a .b'
print *,'!     files), specify either the .a - file or the .b - file when the program'
print *,'!     asks for a global restart file.'
print *,'! 10. The output is written to local.(name of global restart files).'
print *,'!     Tecplot files are also produced for diagnostics '
print *,'!'
print *,'! NB: The curviint program assumes the number of levels and the densities are'
print *,'! the same in the local and the global restart files.'
print *,'! GOOD LUCK.'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
stop '(curviint)'
   end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Need to parse local blkdat.input to get kdm, sigma flag, kapflg, thbase
   open (unit=10,file='blkdat.input',status='old',form='formatted')
   cvarin=''
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   do while (cvarin/='kdm   ')
      read(10,*) rvar,cvarin
   end do
   if (cvarin/='kdm   ') then
      write(6,*) 'ERROR: Could not get idm'
      stop '(tracerinsert:blkdat.input parse)'
   else
      kdm=nint(rvar)
   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get grid size dimensions 
   write(*,*)'Input local grid dimensions (nx, ny) :'
   read(*,*) lnx, lny
   allocate(ldepths(lnx,lny))
   allocate(llon   (lnx,lny))
   allocate(llat   (lnx,lny))
   allocate(ldummy (lnx,lny))
   allocate(lmsk   (lnx,lny))
   allocate(donep  (lnx,lny))
   allocate(io1    (lnx,lny))
   allocate(io2    (lnx,lny))
   allocate(lfld   (lnx,lny))
   allocate(lpbot   (lnx,lny))
   allocate(lthkk   (lnx,lny))
   allocate(lpsikk  (lnx,lny))
   allocate(ltemp  (lnx,lny))
   allocate(pass   (lnx,lny))
   allocate(sumdeep(lnx,lny))
   allocate(loldfld(lnx,lny))
   allocate(index_filled(lnx,lny,2))
   lpbot=0.
   lthkk=0.
   lpsikk=0.
   

   if (lnx>999 .or. lny > 999 ) then
      write(tag7,'(i5.5,a,i5.5)')lnx,'x',lny
   else if (lnx>0 .and. lny>0) then
      write(tag7,'(i3.3,a,i3.3)')lnx,'x',lny
   else
      print *,'Negative grid size, I quit..'
      stop '(curvvint)'
   end if

   inquire(file='ldepths'//trim(tag7)//'.uf',exist=ex)
   if (.not.ex) then
      print *,'local depths file does not exist ldepths'//trim(tag7)//'.uf'
      stop '(curviint)'
   endif
   open (unit=10,file='ldepths'//trim(tag7)//'.uf',status='old',form='unformatted')
   read(10)io1
   close(10)
   ldepths=io1

   inquire(file='lnewpos.uf',exist=ex)
   if (.not.ex) stop 'local newpos.uf file does not exist'
   open (unit=10,file='lnewpos.uf',status='old',form='unformatted')
   !read(10)llat,llon
   read(10)io1,io2
   close(10)
   llat=io1
   llon=io2

   lmsk=.false.
   where (ldepths > 0.0) lmsk=.true.

   ! None of the grid points in the local grid have been updated yet
   donep=.false.

   deallocate(io1)
   deallocate(io2)

   call zaiost_local()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading global depths and newpos file

   ! Get grid size dimensions 
   write(*,*)'Input global grid dimensions (nx, ny) :'
   read(*,*) gnx, gny
   allocate(gdepths(gnx,gny))
   allocate(gmsk   (gnx,gny))
   allocate(gdummy (gnx,gny))
   allocate(gtemp  (gnx,gny))
   allocate(glon   (gnx,gny))
   allocate(glat   (gnx,gny))
   allocate(io1    (gnx,gny))
   allocate(io2    (gnx,gny))
   allocate(gfld   (gnx,gny))
   allocate(gpbot  (gnx,gny))
   allocate(gthkk  (gnx,gny))
   allocate(gpsikk (gnx,gny))
   gpbot=0.
   gthkk=0.
   gpsikk=0.
   

   if (gnx>999 .or. gny > 999 ) then
      write(tag7,'(i5.5,a,i5.5)')gnx,'x',gny
   else if (gnx>0 .and. gny>0) then
      write(tag7,'(i3.3,a,i3.3)')gnx,'x',gny
   else
      print *,'Negative grid size, I quit..'
      stop '(curvvint)'
   end if

   write(tag7,'(i3.3,a,i3.3)')gnx,'x',gny
   inquire(file='gdepths'//trim(tag7)//'.uf',exist=ex)
   if (.not.ex) then
      print *,'global depths file does not exist  gdepths',tag7,'.uf'
      stop '(curviint)'
   endif
   open (unit=10,file='gdepths'//trim(tag7)//'.uf',status='old',form='unformatted')
   read(10)io1
   close(10)
   gdepths=io1

   inquire(file='gnewpos.uf',exist=ex)
   if (.not.ex) stop 'global gnewpos.uf file does not exist'
   open (unit=10,file='gnewpos.uf',status='old',form='unformatted')
   !read(10)glat,glon
   read(10)io1,io2
   close(10)
   glat=io1
   glon=io2

   gmsk=.false.
   where (gdepths > 1e28) gdepths=0.
   where (gdepths > 0.0) gmsk=.true.

   deallocate(io1)
   deallocate(io2)

   call zaiost_global()

   call initconfmap(gnx,gny)



   ! Se if global grid id periodic
   lperiodic=.false.
   if (any(gdepths(1,:) > 0.1) .and. any(gdepths(gnx,:) > 0.1) ) then
      print *,'Periodic global grid '
      lperiodic=.true.
   end if





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set up interpolation arrays
   allocate(ipiv    (lnx,lny))
   allocate(jpiv    (lnx,lny))
   allocate(ass     (lnx,lny))
   allocate(a       (lnx,lny,4))

   ! Interpolate details (ipiv, jpiv, a, ass
   donep=.false.
   call bilin_confgrd(llon,llat,lmsk,donep,     & 
                      glon,glat,gmsk,gdepths,   & 
                      ipiv,jpiv,a,ass,lperiodic)

   do j=1,lny
   do i=1,lnx
      if (donep(i,j)) then
         pass(i,j)=1
      else
         pass(i,j)=0
      end if
   end do
   end do
   print *,'bilin_confgrd done'
   !print *,maxval(ipiv)
   !print *,maxval(jpiv)
   !print *,minval(ipiv)
   !print *,minval(jpiv)
   !print *,11,96,lmsk(11,96),llon(11,96),llat(11,96),ipiv(11,96),jpiv(11,96),pass(11,96)


   ! Now remaining wet points except boundary points are filled in by extrapolation
   call extrapolate2(lmsk,donep,ldepths,ipiv,jpiv,a,lperiodic)
   do j=1,lny
   do i=1,lnx
      if (donep(i,j).and.pass(i,j)==0) then
         pass(i,j)=2
      end if
   end do
   end do
   write(*,*)'extrapolate done'
   !print *,11,96,lmsk(11,96),llon(11,96),llat(11,96),ipiv(11,96),jpiv(11,96),pass(11,96)
   !stop '(test)'
   !print *,maxval(ipiv)
   !print *,maxval(jpiv)
   !print *,minval(ipiv)
   !print *,minval(jpiv)


   ! Check if any points are not properly set up
   if (any((.not.donep).and.lmsk)) then
      print '(a,i6,a)','failed to set up ',count((.not.donep).and.lmsk), &
      ' points .. Check for isolated points on the bathymetry grid '
      print '(a)','See error messages above for candidate points'
      stop '(curviint)'
   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read restart file from old grid

   write(*,*)'Input global restart file :'
   read(*,*) grestart
   grestart=adjustl(grestart)
   grestart=trim(grestart)

   ! Check for .a - .b ending
   find=max(index(grestart,'.a')-1,index(grestart,'.b')-1)
   !print *,find
   if  (find<1) find = len_trim(grestart)
   grestart=grestart(1:find)


   inquire(file=trim(grestart)//'.a',exist=ex)
   if (.not.ex) then
      print *,'restart file '//trim(grestart)//'.a'//' does not exist'
      stop '(curviint)'
   endif

   inquire(file=trim(grestart)//'.b',exist=ex)
   if (.not.ex) then
      print *,'restart file '//trim(grestart)//'.b'//' does not exist'
      stop '(curviint)'
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Scan through records in global file - retrieve layer thickness 
   print *,'Global grid size ',kdm,gnx,gny
   allocate (dp(gnx,gny,kdm,2))
   call zaiopf_global(trim(grestart)//'.a','old',gnop)
   open(gnop,file=trim(grestart)//'.b',status='old')

   do irec=1,2
      read (gnop,'(a80)') a80
   end do
   ios =0

   do while (ios==0)
      read (gnop,'(a80)', iostat=ios ) a80
      read(a80,'(a8, 22x, i4,i3)') cfld,lcoord,istep
      read(a80(38:80),*) bxmin,bxmax
      if (trim(cfld)=='dp'.and. lcoord <=kdm .and. istep <=2) then
         call zaiord_global(dp(:,:,lcoord,istep),gdummy,.false.,axmin,axmax,gnop)
         !print *,maxval(dp(:,:,lcoord,istep)/onem)
      else
         call zaiosk_global(gnop)
      end if
   end do

   call zaiocl_global(gnop)
   close(gnop)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Scan through records in input file - First pass fills water column where
!local model is shallower than the global model

   ! Open global restart file
   call zaiopf_global(trim(grestart)//'.a','old',gnop)
   open(gnop,file=trim(grestart)//'.b',status='old')

   ! Open local file
   call zaiopf_local('pass1.a','replace',lnop)
   open(lnop,file='pass1.b',status='replace')



   ! Copy old restart header stuff to new restart
   do irec=1,2
      read (gnop,'(a80)') a80
      write(lnop,'(a80)') a80
   end do


   irec=1
   ios=0
   index_filled=0
   cfld=''
   oldcfld=''
   lcoord=0
   oldlcoord=0
   oldistep=0
   istep=0
   print *,'First Pass:'
   do while (ios==0) 

    
      ! Read record header (.b file)
      read (gnop,'(a80)', iostat=ios ) a80
      if (ios==0) then
         oldcfld=cfld
         oldlcoord=lcoord
         oldistep=istep
         read(a80,'(a8, 22x, i4,i3)') cfld,lcoord,istep
         read(a80(38:80),*) bxmin,bxmax
         !print *,cfld,lcoord,istep


         if (cfld==oldcfld.and.oldistep==istep) then
            write(*,'(a1)',advance='no') '.'
         else
            if (len_trim(oldcfld)/=0) write(*,'(i2.2)') oldlcoord
            print *,cfld,lcoord,istep
            write(*,'(a8,"istep=",i1," level=",i1)',advance='no') cfld,istep,lcoord
         end if
         !print *,'No velrot, No pres corr, var:',cfld, lcoord
         !print *
         !print *,cfld,lcoord,istep

         ! Read record data (.a file)
         call zaiord_global(gfld,gdummy,.false.,axmin,axmax,gnop)
         !print *,'hei'
         !print *,11,96,lmsk(11,96),llon(11,96),llat(11,96),ipiv(11,96),jpiv(11,96)
         !print *,lperiodic
         !stop '(test)'

         ! Interpolate and  extrapolate
         do j=1,lny
         do i=1,lnx
         if (lmsk(i,j)) then
            i2=ipiv(i,j)
            j2=jpiv(i,j)
            if (lperiodic) then
               i2pb=mod(i2,gnx)+1
            else
               i2pb=i2+1
            end if
            j2pb=j2+1

            if (i2pb/=i2+1) then
               print *,'NB:periodic grid set up but needs testing'

               stop
            end if


            !print *,i,j,lnx,lny,i2,j2,cfld,lmsk(i,j)
            
            if (trim(cfld)/='temp' .and. trim(cfld)/='saln' .and.  trim(cfld)/='th3d') then
            !if (.true.) then

               lfld(i,j)=                     &
                  gfld(i2  ,j2  )*a(i,j,1) +  &
                  gfld(i2pb,j2  )*a(i,j,2) +  &
                  gfld(i2pb,j2pb)*a(i,j,3) +  &
                  gfld(i2  ,j2pb)*a(i,j,4) 
               !print *,i2,j2,a(i,j,:)

               if (trim(cfld)=='dp' .and. lfld(i,j)<0.) then
                  print *,'***Got negative dp ... ',i,j,lcoord
                  print *,lfld(i,j)
                  print *,a(i,j,:)
                  print *,gfld(i2,j2)
                  print *,gfld(i2pb,j2)
                  print *,gfld(i2pb,j2pb)
                  print *,gfld(i2,j2pb)
                  !stop '(curviint)'
                  lfld(i,j)=0.
               end if


            else ! temp, salt and density needs to be layer weighted as well

               if (lcoord>kdm .or. istep>2) then
                  print *
                  print *,'global coord exceed those in blkdat!'
                  print *,'lcoord and kdm   ',lcoord,kdm
                  print *,'istep  and nstep ',istep,2
                  stop
               end if

               ba1=a(i,j,1)*dp(i2  ,j2  ,lcoord,istep)
               ba2=a(i,j,2)*dp(i2pb,j2  ,lcoord,istep)
               ba3=a(i,j,3)*dp(i2pb,j2pb,lcoord,istep)
               ba4=a(i,j,4)*dp(i2  ,j2pb,lcoord,istep)
               basum=ba1+ba2+ba3+ba4
               if (basum>=onem) then ! Layer weighted value
                  lfld(i,j)=                     &
                     gfld(i2  ,j2  )*ba1/basum +  &
                     gfld(i2pb,j2  )*ba2/basum +  &
                     gfld(i2pb,j2pb)*ba3/basum +  &
                     gfld(i2  ,j2pb)*ba4/basum 

               else if (basum<onem.and.lcoord>1) then ! Use layer value above
                  lfld(i,j)=loldfld(i,j)
               else if (basum<onem) then ! abort
                  print *,'basum error !'
                  print *,i,j,lcoord,istep
                  print *,i2,i2pb,j2,j2pb,gdepths(i2,j2)
                  print *,'pass ',pass(i,j)
                  print *,'****'
                  print *,a(i,j,1),dp(i2  ,j2  ,lcoord,istep)
                  print *,a(i,j,2),dp(i2pb,j2  ,lcoord,istep)
                  print *,a(i,j,3),dp(i2pb,j2pb,lcoord,istep)
                  print *,a(i,j,4),dp(i2  ,j2pb,lcoord,istep)
                  !
                  print *,a(i,j,1),dp(i2  ,j2  ,lcoord,mod(istep+1,2))
                  print *,a(i,j,2),dp(i2pb,j2  ,lcoord,mod(istep+1,2))
                  print *,a(i,j,3),dp(i2pb,j2pb,lcoord,mod(istep+1,2))
                  print *,a(i,j,4),dp(i2  ,j2pb,lcoord,mod(istep+1,2))
                  print *,maxval(dp(:,:,lcoord,istep))
                  stop '(curviint)'
               end if

               

            end if
         end if
         end do
         end do


         ! Correct pressure fields - prevents layers from going below sea floor
         if (trim(cfld)=='dp' .and. lcoord==1) then
            where (lmsk) sumdeep=max(0.,min(lfld,ldepths*onem))
            where (lmsk) lfld=sumdeep
         else if (trim(cfld)=='dp') then
            !where (lmsk) lfld=min(lfld,max(0.,ldepths*onem-sumdeep))
            where (lmsk) lfld=max(0.,min(lfld,ldepths*onem-sumdeep))
            where (lmsk) sumdeep=sumdeep+lfld
!if (istep==1) then
!print *,'Max diff depths and layer thickness sum:',maxval(ldepths*onem-sumdeep,mask=lmsk)/onem
!print *,'Min diff depths and layer thickness sum:',minval(ldepths*onem-sumdeep,mask=lmsk)/onem
!endif
         end if

         ! Tag final massfilled layer in model
         if (trim(cfld)=='dp') then
            where (lfld>onem .and. lmsk) index_filled(:,:,istep)=lcoord
         end if
         where(.not. lmsk) lfld=0.


         ! Set velocities to zero
         if (trim(cfld)=='u'     .or. &
             trim(cfld)=='v'     .or. &
             trim(cfld)=='ubavg' .or. &
             trim(cfld)=='vbavg'  ) then
             lfld=0.
         end if


         ! Arghh.. Pbot corrections here...
         if (trim(cfld)=='pbot') then
            gpbot=gfld
         !elseif (trim(cfld)=='psikk') then
         !   gpsikk=gfld
         !elseif (trim(cfld)=='thkk') then
         !   gthkk=gfld
         end if

         call zaiowr_local(lfld,ldummy,.false.,axmin,axmax,lnop,.true.)
         write(lnop,4100) cfld,lcoord,istep,axmin,axmax

         irec=irec+1

         if (trim(cfld)=='temp'.and.lcoord==1.and.istep==1) then
            gtemp=gfld
            ltemp=lfld
         end if

         ! Keep old gfld
         loldfld=lfld
      end if
   end do
   write(*,*)
   write(*,*)
   call zaiocl_global(gnop)
   close(gnop)
   call zaiocl_local (lnop)
   close(lnop)


   print *,'Max diff depths and gpbot sum:',maxval(gdepths*onem-gpbot,mask=gmsk)/onem
   print *,'Min diff depths and gpbot sum:',minval(gdepths*onem-gpbot,mask=gmsk)/onem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Second pass fills water column where global model is DEEPER than local model
! (extends deepest layer)

   ! Open local file
   call zaiopf_local('pass1.a','old',lnop)
   open(lnop,file='pass1.b',status='old')
   call zaiopf_local('local.'//trim(grestart)//'.a','replace',lnop2)
   open(lnop2,file='local.'//trim(grestart)//'.b',status='replace')

   open(789,file='curvisec.tec')
   WRITE(789,*)'TITLE="curvitest"'
   WRITE(789,*)'VARIABLES=lon,lat,i,j,depth'

   ! Copy old restart header stuff to new restart
   do irec=1,2
      read (lnop,'(a80)') a80
      write(lnop2,'(a80)') a80
   end do


   irec=1
   ios=0
   cfld=''
   oldcfld=''
   lcoord=0
   oldlcoord=0
   oldistep=0
   istep=0
   print *,'Second pass'
   do while (ios==0) 

    
      ! Read record header (.b file)
      read (lnop,'(a80)', iostat=ios ) a80
      if (ios==0) then
         oldcfld=cfld
         oldlcoord=lcoord
         oldistep=istep
         read(a80,'(a8, 22x, i4,i3)') cfld,lcoord,istep
         read(a80(38:80),*) bxmin,bxmax


         if (cfld==oldcfld.and.oldistep==istep) then
            write(*,'(a1)',advance='no') '.'
         else
            if (len_trim(oldcfld)/=0) write(*,'(i2.2)') oldlcoord
            write(*,'(a8,"istep=",i1," level=",i1)',advance='no') cfld,istep,lcoord
         end if
         !print *,'No velrot, No pres corr, var:',cfld, lcoord

         ! Read record data (.a file)
         call zaiord_local(lfld,gdummy,.false.,axmin,axmax,lnop)


         ! Fill last massfilled layer to ocean floor
         if (trim(cfld)=='dp') then

            ! Accumulated depth
            if (lcoord==1) then
               sumdeep=0.
            end if

            do j=1,lny
            do i=1,lnx
            if (lmsk(i,j)) then

               !print *,i,j,index_filled(i,j,istep)

               ! Fill last mass-filled layer to ocean floor
               if (lcoord==index_filled(i,j,istep)) then
                  lfld(i,j) = max(0.,ldepths(i,j)*onem-sumdeep(i,j))
               else if (lcoord>index_filled(i,j,istep)) then
                  lfld(i,j) = 0.
               end if
            end if
            end do
            end do

            ! Accumulated depth
            if (lcoord==1) then
               sumdeep=lfld
            else
               sumdeep=sumdeep+lfld
            end if
!if (istep==1) then
!print *,'Max diff depths and layer thickness sum:',maxval(ldepths*onem-sumdeep,mask=lmsk)/onem
!print *,'Min diff depths and layer thickness sum:',minval(ldepths*onem-sumdeep,mask=lmsk)/onem
!endif

            ! Dump to section
            if (istep==1) then
               WRITE(789,*)'ZONE I=',lnx,',J=',1,',F=POINT'
               if (lcoord==1) then 
                  do i=1,lnx
                     WRITE(789,*) i,j,llon(i,lny/2),llat(i,lny/2),sumdeep(i,lny/2)/onem
                  end do
               else
                  WRITE(789,*)'D=(1,2,3,4)'
                  do i=1,lnx
                     WRITE(789,*) sumdeep(i,lny/2)/onem
                  end do
               end if
            end if
         end if

         ! Arghh.. Pbot corrections here...
         if (trim(cfld)=='pbot') then
            ! pbot MUST be equal to depths !!
            lfld=ldepths*onem 
            lpbot=lfld
         elseif (trim(cfld)=='psikk') then
            lpsikk=lfld
         elseif (trim(cfld)=='thkk') then
            lthkk=lfld
         end if

         ! constrain mixed layer 
         if (trim(cfld)=='dpmixl') then
           lfld=min(ldepths*onem,lfld)
         end if



         call zaiowr_local(lfld,ldummy,.false.,axmin,axmax,lnop2,.true.)
         write(lnop2,4100) cfld,lcoord,istep,axmin,axmax
      end if
   end do
   call zaiocl_local (lnop)
   close(lnop)
   call zaiocl_local (lnop2)
   close(lnop2)
   close(789)

   ! Sumdeep should now add up to depths...
   print *
   print *
   print *,'Max diff depths and layer thickness sum:',maxval(ldepths*onem-sumdeep,mask=lmsk)/onem
   print *,'Maxloc diff depths and layer thickness sum:',maxloc(ldepths*onem-sumdeep,mask=lmsk)
   !print *,index_filled(203,385,2)
   print *,'Min diff depths and layer thickness sum:',minval(ldepths*onem-sumdeep,mask=lmsk)/onem
   print *,'Max diff depths and pbot sum:',maxval(ldepths*onem-lpbot,mask=lmsk)/onem
   print *,'Min diff depths and pbot sum:',minval(ldepths*onem-lpbot,mask=lmsk)/onem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Third pass (NEW) - replace second time index with values from first.
! Interpolation procedure is nonlinear and strong gradientsmay appear
! because of this
! (extends deepest layer)

   print *
   print *,'NB: Third/fourth/fourth pass not implemented yet'
   print *,'third pass need to initialize thkk,psikk from climatology,'
   print *,'now they are interpolated from global model'
   print *,'fourth pass -- time level 2 = timelevel 1...'



   ! The ice stuff follows
   allocate(gficem(gnx,gny))
   allocate(ghicem(gnx,gny))
   allocate(ghsnwm(gnx,gny))
   allocate(gticem(gnx,gny))
   allocate(gtsrfm(gnx,gny))
   !
   allocate(lficem(lnx,lny))
   allocate(lhicem(lnx,lny))
   allocate(lhsnwm(lnx,lny))
   allocate(lticem(lnx,lny))
   allocate(ltsrfm(lnx,lny))

   print *
   print *
   write(*,*)'Input global restart file for ice :' ; read(*,*) gicefname
   inquire(exist=ex, file=trim(gicefname))
   if (ex) then 
      licefname='local.'//trim(gicefname)
      inquire(iolength=grecl) gficem,ghicem,ghsnwm,gticem,gtsrfm
      inquire(iolength=lrecl) lficem,lhicem,lhsnwm,lticem,ltsrfm
      open(89,file=trim(gicefname),access='direct',status='old',recl=grecl)
      open(90,file=trim(licefname),access='direct',status='replace',recl=lrecl)



      read(89,rec=1) gficem, ghicem, ghsnwm, gticem, gtsrfm

      ! Interpolate and  extrapolate
      do j=1,lny
      do i=1,lnx
      if (lmsk(i,j)) then
         i2=ipiv(i,j)
         j2=jpiv(i,j)
         !print *,i2,j2
         lficem(i,j)=gficem(i2  ,j2  )*a(i,j,1) + gficem(i2+1,j2  )*a(i,j,2) +  &
                     gficem(i2+1,j2+1)*a(i,j,3) + gficem(i2  ,j2+1)*a(i,j,4) 
         lhicem(i,j)=ghicem(i2  ,j2  )*a(i,j,1) + ghicem(i2+1,j2  )*a(i,j,2) +  &
                     ghicem(i2+1,j2+1)*a(i,j,3) + ghicem(i2  ,j2+1)*a(i,j,4) 
         lhsnwm(i,j)=ghsnwm(i2  ,j2  )*a(i,j,1) + ghsnwm(i2+1,j2  )*a(i,j,2) +  &
                     ghsnwm(i2+1,j2+1)*a(i,j,3) + ghsnwm(i2  ,j2+1)*a(i,j,4) 
         ltsrfm(i,j)=gtsrfm(i2  ,j2  )*a(i,j,1) + gtsrfm(i2+1,j2  )*a(i,j,2) +  &
                     gtsrfm(i2+1,j2+1)*a(i,j,3) + gtsrfm(i2  ,j2+1)*a(i,j,4) 
         lticem(i,j)=gticem(i2  ,j2  )*a(i,j,1) + gticem(i2+1,j2  )*a(i,j,2) +  &
                     gticem(i2+1,j2+1)*a(i,j,3) + gticem(i2  ,j2+1)*a(i,j,4) 
      else
         lficem(i,j)=0.
         lhicem(i,j)=0.
         lhsnwm(i,j)=0.
         ltsrfm(i,j)=0.
         lticem(i,j)=0.
      end if
      end do
      end do

      write(90,rec=1) lficem, lhicem, lhsnwm, lticem, ltsrfm

      write(*,*)'Local ice restart file created :',trim(licefname)

   else
      lficem=0.
      lhicem=0.
      lhsnwm=0.
      lticem=0.
      ltsrfm=0.
   end if
   print *
   print *


   open(789,file='curvitest.tec')
   WRITE(789,*)'TITLE="curvitest"'
   WRITE(789,*)'VARIABLES=i,j,lon,lat,depth,pbot,psikk,thkk,ipiv,jpiv,temp,hicem,pass'
   WRITE(789,*)'ZONE I=',lnx,',J=',lny,',F=BLOCK'
   WRITE(789,101)((i,              i=1,lnx),j=1,lny)
   WRITE(789,101)((j,              i=1,lnx),j=1,lny)
   WRITE(789,100)((llon(i,j),      i=1,lnx),j=1,lny)
   WRITE(789,100)((llat(i,j),      i=1,lnx),j=1,lny)
   WRITE(789,100)((ldepths(i,j),   i=1,lnx),j=1,lny)
   WRITE(789,100)((lpbot  (i,j)/onem,   i=1,lnx),j=1,lny)
   WRITE(789,100)((lpsikk (i,j),   i=1,lnx),j=1,lny)
   WRITE(789,100)((lthkk  (i,j),   i=1,lnx),j=1,lny)
   WRITE(789,101)((ipiv   (i,j),   i=1,lnx),j=1,lny)
   WRITE(789,101)((jpiv   (i,j),   i=1,lnx),j=1,lny)
   WRITE(789,100)((ltemp  (i,j),   i=1,lnx),j=1,lny)
   WRITE(789,100)((lhicem (i,j),   i=1,lnx),j=1,lny)
   !WRITE(789,100)((ltemp  (i,j),   i=1,lnx),j=1,lny)
   WRITE(789,101)((pass(i,j)      ,i=1,lnx),j=1,lny)
   WRITE(789,*)'ZONE I=',gnx,',J=',gny,',F=BLOCK'
   WRITE(789,101)((i,              i=1,gnx),j=1,gny)
   WRITE(789,101)((j,              i=1,gnx),j=1,gny)
   WRITE(789,100)((glon(i,j),      i=1,gnx),j=1,gny)
   WRITE(789,100)((glat(i,j),      i=1,gnx),j=1,gny)
   WRITE(789,100)((gdepths(i,j),   i=1,gnx),j=1,gny)
   WRITE(789,100)((gpbot  (i,j)/onem,   i=1,gnx),j=1,gny)
   WRITE(789,100)((gpsikk  (i,j),   i=1,gnx),j=1,gny)
   WRITE(789,100)((gthkk  (i,j),   i=1,gnx),j=1,gny)
   WRITE(789,101)((i,              i=1,gnx),j=1,gny)
   WRITE(789,101)((j,              i=1,gnx),j=1,gny)
   WRITE(789,100)((gtemp  (i,j),   i=1,gnx),j=1,gny)
   WRITE(789,100)((ghicem (i,j),   i=1,gnx),j=1,gny)
   WRITE(789,101)((1              ,i=1,gnx),j=1,gny)
   close(789)
100   FORMAT(10(1x,e12.6)) 
101   FORMAT(30i4)
4100  format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)

  print *,'New restart files created: local.'//trim(grestart)//'.[ab]'
  if (ex) write(*,*)'Local ice restart file created :',trim(licefname)
  print *,'Diagnostic files created (tecplot): curvitest.tec curvisec.tec'

  stop '(normal)'

contains



subroutine bilin_confgrd(llon,llat,lmsk,donep,    &
                         glon,glat,gmsk,gdepths,  &
                         ipiv, jpiv,a,ass,lperiodic)
   use mod_xc_local , only : lnx=>idm,lny=>jdm
   use mod_xc_global, only : gnx=>idm,gny=>jdm
   use mod_confmap
   implicit none
   real,     intent(in)    :: llon(lnx,lny)
   real,     intent(in)    :: llat(lnx,lny)
   logical,  intent(in)    :: lmsk(lnx,lny)
   logical,  intent(out)   :: donep(lnx,lny)
   real,     intent(in)    :: glon(gnx,gny)
   real,     intent(in)    :: glat(gnx,gny)
   real,     intent(in)    :: gdepths(gnx,gny)
   logical,  intent(in)    :: gmsk(gnx,gny)
   logical,  intent(out)   :: ass(lnx,lny)
   real   ,  intent(out)   :: a(lnx,lny,4)
   integer,  intent(out)   :: ipiv(lnx,lny)
   integer,  intent(out)   :: jpiv(lnx,lny)
   logical,  intent( in)   :: lperiodic

   integer l,i,j,ia,ib,k
   integer ipib,jpib
   real aa,bb,atmp(4)
   logical masktmp(4)
   integer :: ip2,jp2
   real lat_n,lon_n
   real dpsum

   real, parameter:: epsil=1.0e-11
   real t1,t2,s1,s2,r1,r2
   real, parameter ::  onem=9806.
   real, parameter ::  onemm=98.06
   integer iloc(2), igrace, jgrace
   logical inirange



   ! For each gridpoint in local grid if wet
   do j=1,lny
   do i=1,lnx
      if (lmsk(i,j)) then 
         t1=llon(i,j)
         t2=llat(i,j)
         !print *,i,j,t1,t2
         call oldtonew(t2,t1,lat_n,lon_n)
         !call newtoold(lat_n,lon_n,s2,s1)
         call pivotp(lon_n,lat_n,ipiv(i,j),jpiv(i,j))

         if (lperiodic) then
            ipib=mod(ipiv(i,j),gnx)+1
            inirange=.true.
         else
            ipib=ipiv(i,j)+1
            inirange=ipiv(i,j)>=1 .and. ipiv(i,j) < gnx
         end if
         jpib=jpiv(i,j)+1

         ! grace for i
         igrace=min(gnx-ipiv(i,j),ipiv(i,j)-1) ! negative when ipiv < 1 or ipiv > nxl
         ! grace for j
         jgrace=min(gny-jpiv(i,j),jpiv(i,j)-1) ! negative when jpiv < 1 or jpiv > nxl

         !! Brutally stop if pivot points are outside of grid
         !if (ipiv(i,j)>gnx-1 .or. ipiv(i,j)<1 .or. &
         !    jpiv(i,j)>gny-1 .or. jpiv(i,j)<1 ) then
         !   print *,'pivot point is outside global grid!'
         !   print *,'pivot point  (i,j) :',ipiv(i,j),jpiv(i,j)
         !   print *,'global grid dim    :',gnx,gny
         !   stop '(bilin_confgrid)'
         !end if




         if (inirange .and. jpiv(i,j)>=1 .and. jpiv(i,j)<gny) then

            if (gmsk(ipiv(i,j),jpiv(i,j)) .and. gmsk(ipib,jpiv(i,j)) .and.&
                gmsk(ipiv(i,j),jpib     ) .and. gmsk(ipib,jpib    )) then
               call bilincoeff(glon,glat,gnx,gny,llon(i,j),llat(i,j), &
                               ipiv(i,j),jpiv(i,j),atmp(1),atmp(2),atmp(3), &
                               atmp(4))
               a(i,j,:)=atmp
               ass(i,j)=.true.
            else
               ass(i,j)=.false. 
               a(i,j,:)=0
               ipiv(i,j)=0; jpiv(i,j)=0
            endif

         else if ( igrace> -3 .and. jgrace>-3) then
            print '(a,2i5,a,2i5,a)','Point ',i,j,'with pivot ',ipiv(i,j),jpiv(i,j),' saved by grace '
            ass(i,j)=.false. 
            a(i,j,:)=0
            ipiv(i,j)=0; jpiv(i,j)=0
         else
            print '(a,2i5,a,2i5,a)','Point ',i,j,'with pivot ',ipiv(i,j),jpiv(i,j),' not saved by grace '
            stop '(bilin_confgrd)'
         endif
      end if
      if (ass(i,j)) then 
         donep(i,j)=.true.
      endif
   enddo
   enddo
end subroutine bilin_confgrd


subroutine extrapolate2(lmsk,donep,ldepths,ipiv,jpiv,a,lperiodic)
use mod_xc_local , only : lnx=>idm,lny=>jdm
implicit none
logical, intent(in)    :: lmsk   (lnx,lny)
logical, intent(inout) :: donep  (lnx,lny)
real   , intent(in)    :: ldepths(lnx,lny)
real   , intent(inout) :: a      (lnx,lny,4)
integer, intent(inout) :: ipiv   (lnx,lny)
integer, intent(inout) :: jpiv   (lnx,lny)
logical, intent(in)    :: lperiodic
integer l,i,j,ia,ib,isign,ifalse,jja,iia,index,k,ja
logical donetmp(lnx,lny)
! Now remaining wet points except bounday points are filled in by extrapolation

   do l=1,200
      if (mod(l,2) == 0) isign=1
      if (mod(l,2) == 1) isign=-1

      ifalse=0
      donetmp=.false.
      do j=1,lny
      do i=1,lnx
         if (lmsk(i,j).and.(.not.donep(i,j))) then   
            ifalse=ifalse+1
            do jja=-1,1
               ja=j+jja*isign

               !if (ja < 1) ja=lny
               !if (ja > lny) ja=1
               if (ja < 1) ja=1
               if (ja > lny) ja=lny

               do iia=-1,1
                  ia=i+iia*isign

                  if (lperiodic) then
                     ia=mod(lnx+ia-1,lnx)+1
                  else
                     if (ia > lnx) ia=lnx
                     if (ia < 1) ia=1
                  end if

                  !if (i==11.and.j==96) print *,ia,ja

                  if (lmsk(ia,ja).and.donep(ia,ja)) then
                     ipiv(i,j)  =ipiv(ia,ja)
                     jpiv(i,j)  =jpiv(ia,ja)
                     a   (i,j,:)= a (ia,ja,:)
                     donetmp(i,j)=.true.
                     ifalse=ifalse-1
                     exit
                  endif
               enddo
               if (donetmp(i,j)) exit
            enddo
         endif
      enddo
      enddo
      print *,'ifalse after iteration',l,ifalse
      where (donetmp) donep=.true.
      if (ifalse == 0) exit
   enddo
   
   if (ifalse /= 0) then
      do j=1,lny
      do i=1,lnx
         if (lmsk(i,j).and.(.not.donep(i,j))) then
            print '(a,2I6)',' Ifalse problem in grip point (i,j): ',i,j
         endif
      enddo
      enddo
   endif
end subroutine extrapolate2



!KALsubroutine chk_pressure(dp,lmsk,donep,depths)
!KAL   use mod_xc , only : lnx,lny,nz
!KAL   implicit none
!KAL   real,    intent(inout) :: dp(lnx,lny,nz)
!KAL   logical, intent(in)    :: lmsk(lnx,lny)
!KAL   real,    intent(in)    :: depths(lnx,lny)
!KAL   logical, intent(in)    :: donep(lnx,lny)
!KAL   real pres(lnx,lny,nz+1)
!KAL   real, parameter ::  onem=9806.
!KAL   integer i,j,k,index
!KAL
!KAL   pres(:,:,1)=0.0
!KAL
!KAL   do j=1,lny
!KAL   do i=1,lnx
!KAL      if (lmsk(i,j).and.donep(i,j)) then   ! wet point
!KAL
!KAL         pres(i,j,1)=0.0
!KAL         do k=1,nz
!KAL            pres(i,j,k+1)=pres(i,j,k)+dp(i,j,k)
!KAL            pres(i,j,k+1)=max(onem,min(depths(i,j)*onem,pres(i,j,k+1)))
!KAL         enddo
!KAL
!KAL!ab
!KAL         do k=1,nz
!KAL            dp(i,j,k)=pres(i,j,k+1)-pres(i,j,k)
!KAL         enddo
!KAL!ab
!KAL
!KAL         if (pres(i,j,nz+1) < depths(i,j)*onem) then
!KAL            write(41,'(a,2I4,60F8.1)')'A',i,j,depths(i,j),pres(i,j,:)/onem
!KAL            index=0
!KAL            do k=nz,1,-1
!KAL               if (pres(i,j,k+1)-pres(i,j,k) > 2.0*onem) then
!KAL                  index=k+1
!KAL                  exit
!KAL               endif
!KAL            enddo
!KAL
!KAL            if (index > 0) then
!KAL               do k=index,nz+1
!KAL                  pres(i,j,k)=depths(i,j)*onem
!KAL               enddo
!KAL
!KAL               do k=1,nz
!KAL                  dp(i,j,k)=pres(i,j,k+1)-pres(i,j,k)
!KAL               enddo
!KAL
!KAL!              write(41,'(a,2I8)')'chk---vertical press processed : ',i,j
!KAL            else
!KAL!               pres(i,j,nz+1)= depths(i,j)*onem
!KAL               write(41,'(a,2I8)')'chk---vertical press problem in: ',i,j
!KAL               write(41,'(a,2I4,60F8.1)')'B',i,j,depths(i,j),pres(i,j,:)/onem
!KAL            endif
!KAL            write(41,'(a,2I4,60F8.1)')'C',i,j,depths(i,j),pres(i,j,:)/onem
!KAL            write(41,*)
!KAL         endif
!KAL      endif
!KAL   enddo
!KAL   enddo
!KAL
!KALend subroutine chk_pressure
    end program curviint
