program curviint
   use mod_xc_local , only : lnx=>idm,lny=>jdm, &
                             xcspmd_local => xcspmd
   use mod_xc_global, only : gnx=>idm,gny=>jdm, &
                             xcspmd_global => xcspmd
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
   use mod_curviint_interp
   use mod_confmap
   use m_parse_blkdat
   implicit none

   real, parameter :: huge=2.0**99

   ! Local model - 
   real, dimension(:,:), allocatable ::    &
      llat, llon, ldepths, sumdeep
   logical, dimension(:,:), allocatable :: &
      lmsk, gmsk
   integer, dimension(:,:), allocatable :: &
      ldummy, gdummy
   integer, dimension(:,:,:), allocatable :: &
      index_filled
   real*8 , dimension(:,:), allocatable :: &
      gficem,ghicem,ghsnwm,gticem,gtsrfm,  &
      lficem,lhicem,lhsnwm,lticem,ltsrfm
   real, dimension(:,:), allocatable ::   &
      tmp, glat, glon, gdepths,           &
      gfld,gtemp,                         &
      lfld,ltemp,sumdepth,lpbot,gpbot,    &
      gpsikk,lpsikk,gthkk,lthkk,loldfld
   real, dimension(:,:,:,:), allocatable :: dp

   character(len=100) grestart, lrestart
   character(len=80 ) a80,licefname,gicefname
   character(len=11 ) tag7
   character(len=8) :: cfld,oldcfld
   character(len=6) :: cvarin

   integer i,j, lrecl,grecl
   integer :: nirec,nhrec,irec,ios, indx
   real, parameter ::  onem=9806.
   real :: bxmin, bxmax, axmin, axmax
   integer :: itime, lcoord,find,istep,oldlcoord,oldistep
   integer,parameter :: gnop=333, lnop=444, lnop2=555
   real time0
   logical ex, exb
   logical lperiodic
   integer :: kdm
#if defined(IARGC)
   integer*4, external :: iargc
#endif

   if (iargc()==1 .or. iargc()==2) then
      call getarg(1,grestart)
      if (iargc()==2) then
         call getarg(2,gicefname)
      end if
   else
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
      print *,'!  4. The global depth file must be available. It is named global.depths.[ab]'
      print *,'!  5. The global grid  file must be available. It is named global.grid.[ab]'
      print *,'!  6. The local  grid must be available.  It is named regional.grid.[ab].'
      print *,'!  7. The local  depth must be available. It is named regional.depth.[ab]'
      print *,'!  8. global and local grid dimensions are extracted from grid files'
      print *,'!  9. The input global restart file to be used is specified as an argument.'
      print *,'!     ICE restart files can also be specified'
      print *,'! 10. The output is written to curviint.(name of global restart files).'
      print *,'!'
      print *,'! NB: The curviint program assumes the number of levels and the densities are'
      print *,'! the same in the local and the global restart files.'
      print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print *,'Usage : '
      print *,'    curviint-2.2 restartfile [ice restart file]'
      stop '(curviint)'
   end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Need to parse local blkdat.input to get kdm, sigma flag, kapflg, thbase
   call parse_blkdat('kdm   ',kdm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Init grid size of local grid
   call xcspmd_local()
   allocate(ldepths(lnx,lny))
   allocate(llon   (lnx,lny))
   allocate(llat   (lnx,lny))
   allocate(ldummy (lnx,lny))
   allocate(lmsk   (lnx,lny))
   allocate(lfld   (lnx,lny))
   allocate(lpbot   (lnx,lny))
   allocate(lthkk   (lnx,lny))
   allocate(lpsikk  (lnx,lny))
   allocate(ltemp  (lnx,lny))
   allocate(sumdeep(lnx,lny))
   allocate(loldfld(lnx,lny))
   allocate(index_filled(lnx,lny,2))
   lpbot=0.
   lthkk=0.
   lpsikk=0.
   
   ! Read local grid 
   call zaiost_local()
   call zaiopf_local('regional.depth.a','old',lnop)
   call zaiord_local(ldepths,ldummy,.false.,axmin,axmax,lnop)
   call zaiocl_local(lnop)
   !
   call zaiopf_local('regional.grid.a','old',lnop)
   call zaiord_local(llon,ldummy,.false.,axmin,axmax,lnop)
   call zaiord_local(llat,ldummy,.false.,axmin,axmax,lnop)
   call zaiocl_local(lnop)

   lmsk=.false.
   where (ldepths > 0.5*huge) ldepths=0.
   where (ldepths > 0.0) lmsk=.true.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading global depths and grid file
   call xcspmd_global()
   allocate(gdepths(gnx,gny))
   allocate(gmsk   (gnx,gny))
   allocate(gdummy (gnx,gny))
   allocate(gtemp  (gnx,gny))
   allocate(glon   (gnx,gny))
   allocate(glat   (gnx,gny))
   allocate(gfld   (gnx,gny))
   allocate(gpbot  (gnx,gny))
   allocate(gthkk  (gnx,gny))
   allocate(gpsikk (gnx,gny))
   gpbot=0.
   gthkk=0.
   gpsikk=0.
   
   ! Read global grid 
   call zaiost_global()
   call zaiopf_global('global.depth.a','old',gnop)
   call zaiord_global(gdepths,gdummy,.false.,axmin,axmax,gnop)
   call zaiocl_global(gnop)
   !
   call zaiopf_global('global.grid.a','old',gnop)
   call zaiord_global(glon,gdummy,.false.,axmin,axmax,gnop)
   call zaiord_global(glat,gdummy,.false.,axmin,axmax,gnop)
   call zaiocl_global(gnop)

   gmsk=.false.
   where (gdepths > 0.5*huge ) gdepths=0.
   where (gdepths > 0.0      ) gmsk=.true.


   call initconfmap(gnx,gny)

   ! Se if global grid is periodic
   lperiodic=.false.
   if (any(gdepths(1,:) > 0.1) .and. any(gdepths(gnx,:) > 0.1) ) then
      print *,'Periodic global grid '
      lperiodic=.true.
   end if


!  Set up interpolation module after grid size and masks are established
   call interp_setup(lmsk,llon,llat,gmsk,glon,glat,lperiodic)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read restart file from old grid Check for .a - .b ending
   find=max(index(grestart,'.a')-1,index(grestart,'.b')-1)
   if  (find<1) find = len_trim(grestart)
   grestart=grestart(1:find)

   inquire(file=trim(grestart)//'.a',exist=ex)
   inquire(file=trim(grestart)//'.a',exist=exb)
   if (.not.ex .or. .not. exb) then
      print *,'restart file '//trim(grestart)//'.[ab]'//' does not exist'
      stop '(curviint)'
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Scan through records in global file - retrieve layer thickness 
   print *,'Global grid size ',kdm,gnx,gny
   allocate (dp(gnx,gny,kdm,2))
   call zaiopf_global(trim(grestart)//'.a','old',gnop)
   open(gnop,file=trim(grestart)//'.b',status='old')

   ! Skip header
   do irec=1,2
      read (gnop,'(a80)') a80
   end do

   ! Search for dp
   ios =0 ; ex=.false.
   do while (ios==0)
      read (gnop,'(a80)', iostat=ios ) a80
      indx=index(a80,'=')
      cfld=a80(1:8)
      read(a80(indx+1:len_trim(a80)),*) lcoord, istep, bxmin, bxmax
      if (trim(cfld)=='dp'.and. lcoord <=kdm .and. istep <=2) then
         ex=.true.
         call zaiord_global(dp(:,:,lcoord,istep),gdummy,.false.,axmin,axmax,gnop)
      else
         call zaiosk_global(gnop)
      end if
   end do
   call zaiocl_global(gnop)
   close(gnop)

   if (.not. ex) then 
      write(6,'(a)') 'Could not find variable dp '
      stop '(curviint)'
   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Scan through records in input file - First pass fills water column where
!local model is shallower than the global model

   ! Open global restart file
   call zaiopf_global(trim(grestart)//'.a','old',gnop)
   open(gnop,file=trim(grestart)//'.b',status='old')

   ! Open temporary restart file
   call zaiopf_local('pass1.a','replace',lnop)
   open(lnop,file='pass1.b',status='replace')

   ! Copy old restart header to new restart
   do irec=1,2
      read (gnop,'(a80)') a80
      write(lnop,'(a80)') a80
   end do

   ios=0
   index_filled=0
   cfld=''
   oldcfld=''
   lcoord=0
   oldlcoord=0
   oldistep=0
   istep=0
   print *
   print *,'First Pass(Horizontal interpolation):'
   do while (ios==0) 

    
      ! Read record header (.b file)
      read (gnop,'(a80)', iostat=ios ) a80
      if (ios==0) then
         oldcfld=cfld
         oldlcoord=lcoord
         oldistep=istep
         indx=index(a80,'=')
         cfld=a80(1:8)
         read(a80(indx+1:len_trim(a80)),*) lcoord, istep, bxmin, bxmax
         if (lcoord>kdm .or. istep>3) then
            print *
            print *,a80
            print *,'global coord exceed those in blkdat!'
            print *,'lcoord and kdm   ',lcoord,kdm
            print *,'istep  and nstep ',istep,2
            stop
         end if

         if (cfld==oldcfld.and.oldistep==istep) then
            write(6,'(a1)',advance='no') '.'
            call flush(6)
         else
            if (len_trim(oldcfld)/=0) write(*,'(i2.2)') oldlcoord
            write(*,'(a8,"istep=",i1," level=",i1)',advance='no') cfld,istep,lcoord
            call flush(6)
         end if

         ! Read record data (.a file)
         call zaiord_global(gfld,gdummy,.false.,axmin,axmax,gnop)

         ! Interpolate/extrapolate
         call bilinear_calc(cfld,lfld,loldfld,lmsk,gfld, &
                             dp(:,:,lcoord,istep),lcoord)

         ! Correct pressure fields - prevents layers from going below sea floor
         if (trim(cfld)=='dp' .and. lcoord==1) then
            where (lmsk) sumdeep=max(0.,min(lfld,ldepths*onem))
            where (lmsk) lfld=sumdeep
         else if (trim(cfld)=='dp') then
            where (lmsk) lfld=max(0.,min(lfld,ldepths*onem-sumdeep))
            where (lmsk) sumdeep=sumdeep+lfld
         end if

         ! Tag final massfilled layer in model
         if (trim(cfld)=='dp') then
            where (lfld>onem .and. lmsk) index_filled(:,:,istep)=lcoord
         end if
         where(.not. lmsk) lfld=0.

         ! Set velocities to zero - model needs spinup at this point anyway
         if (trim(cfld)=='u'     .or. &
             trim(cfld)=='v'     .or. &
             trim(cfld)=='ubavg' .or. &
             trim(cfld)=='vbavg'  ) then
             lfld=0.
         end if

         ! Arghh.. Pbot corrections here...
         if (trim(cfld)=='pbot') then
            gpbot=gfld
         end if

         call zaiowr_local(lfld,ldummy,.false.,axmin,axmax,lnop,.true.)
         write(lnop,4100) cfld,lcoord,istep,axmin,axmax

         if (trim(cfld)=='temp'.and.lcoord==1.and.istep==1) then
            gtemp=gfld
            ltemp=lfld
         end if

         ! Keep old gfld
         loldfld=lfld
      end if
   end do
   write(*,*)
   call zaiocl_global(gnop)
   close(gnop)
   call zaiocl_local (lnop)
   close(lnop)
   print *,'Max diff depths and gpbot sum:',maxval(gdepths*onem-gpbot,mask=gmsk)/onem
   print *,'Min diff depths and gpbot sum:',minval(gdepths*onem-gpbot,mask=gmsk)/onem
   print *


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Second pass fills water column where global model is DEEPER than local model
! (extends deepest layer)

   ! Open local file
   call zaiopf_local('pass1.a','old',lnop)
   open(lnop,file='pass1.b',status='old')
   lrestart='curviint.'//trim(grestart)
   call zaiopf_local(trim(lrestart)//'.a','replace',lnop2)
   open(lnop2,file=trim(lrestart)//'.b',status='replace')

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
   print *,'Second pass(Extend deepest mass-filled layer'
   do while (ios==0) 

    
      ! Read record header (.b file)
      read (lnop,'(a80)', iostat=ios ) a80
      if (ios==0) then
         oldcfld=cfld
         oldlcoord=lcoord
         oldistep=istep
         read(a80,'(a8, 22x, i4,i3)') cfld,lcoord,istep
         read(a80(38:80),*) bxmin,bxmax
         call zaiord_local(lfld,gdummy,.false.,axmin,axmax,lnop)


         if (cfld==oldcfld.and.oldistep==istep) then
            write(*,'(a1)',advance='no') '.'
            call flush(6)
         else
            if (len_trim(oldcfld)/=0) write(*,'(i2.2)') oldlcoord
            write(*,'(a8,"istep=",i1," level=",i1)',advance='no') cfld,istep,lcoord
            call flush(6)
         end if

         ! Fill last massfilled layer to ocean floor
         if (trim(cfld)=='dp') then

            ! Accumulated depth
            if (lcoord==1) then
               sumdeep=0.
            end if

            ! Fill last mass-filled layer to ocean floor
            do j=1,lny
            do i=1,lnx
            if (lmsk(i,j)) then
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

   ! Sumdeep should now add up to depths...
   print *
   print *,'Max diff depths and layer thickness sum:',maxval(ldepths*onem-sumdeep,mask=lmsk)/onem
   print *,'Maxloc diff depths and layer thickness sum:',maxloc(ldepths*onem-sumdeep,mask=lmsk)
   print *,'Min diff depths and layer thickness sum:',minval(ldepths*onem-sumdeep,mask=lmsk)/onem
   print *,'Max diff depths and pbot sum:',maxval(ldepths*onem-lpbot,mask=lmsk)/onem
   print *,'Min diff depths and pbot sum:',minval(ldepths*onem-lpbot,mask=lmsk)/onem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   print *
   print *,'Third pass (initialize thkk,psikk from climatology - NOT IMPLEMENTED)'

!   ! Open local template file
!   call zaiopf_local(trim(template)//'.a','old',lnop2)
!   open(lnop2,file=trim(template)//'.b',status='old')
!   do irec=1,2
!      read (lnop,'(a80)') a80
!   end do
!   ios=0
!   do while (ios==0) 
!      ! Read record header (.b file)
!      read (lnop,'(a80)', iostat=ios ) a80
!      if (ios==0) then
!         read(a80,'(a8, 22x, i4,i3)') cfld,lcoord,istep
!         read(a80(38:80),*) bxmin,bxmax
!         call zaiord_local(lfld,gdummy,.false.,axmin,axmax,lnop)
!         if (trim(cfld)=='thkk') then
!            lthkk=lfld
!         else if (trim(cfld)=='psikk') then
!            lpsikk=lfld
!         end if
!      end if
!   end do
!   call zaiocl_local (lnop2)
!   close(lnop2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   print *,'Fourth pass -- time level 2 = timelevel 1 - NOT IMPLEMENTED'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
   if (trim(gicefname)/='') then
      inquire(exist=ex, file=trim(gicefname))
      if (ex) then 
         licefname='curviint.'//trim(gicefname)
         inquire(iolength=grecl) gficem,ghicem,ghsnwm,gticem,gtsrfm
         inquire(iolength=lrecl) lficem,lhicem,lhsnwm,lticem,ltsrfm
         open(89,file=trim(gicefname),access='direct',status='old',recl=grecl)
         open(90,file=trim(licefname),access='direct',status='replace',recl=lrecl)
         read(89,rec=1) gficem, ghicem, ghsnwm, gticem, gtsrfm

         ! Interpolate/extrapolate
         call bilinear_calc('hicem',lhicem,loldfld,lmsk,ghicem, &
                             dp(:,:,1,1),lcoord)
         call bilinear_calc('ficem',lficem,loldfld,lmsk,gficem, &
                             dp(:,:,1,1),lcoord)
         call bilinear_calc('hsnwm',lhsnwm,loldfld,lmsk,ghsnwm, &
                             dp(:,:,1,1),lcoord)
         call bilinear_calc('ticem',lticem,loldfld,lmsk,gticem, &
                             dp(:,:,1,1),lcoord)
         call bilinear_calc('tsrfm',ltsrfm,loldfld,lmsk,gtsrfm, &
                             dp(:,:,1,1),lcoord)

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
   else
      ex=.false.
   end if


4100  format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)

  print *,'New restart files created: '//trim(lrestart)//'.[ab]'
  if (ex) write(*,*)'Local ice restart file created :',trim(licefname)
  print *,'Restart file can be diagnosed via standard tools'
  stop '(normal)'


end program curviint
