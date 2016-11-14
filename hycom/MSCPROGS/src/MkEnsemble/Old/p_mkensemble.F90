      program mkensemble
! Reads initial model (mem) field from
!   inimem.uf (record 1)
! and initial ice field from
!   iniicemem.uf
! and stores another ensdim-1 fields in records 2,ensdim.
! 1.  make sure inimem.uf and iniicemem.uf exists
! 2.  make sure the directory Intsec exists
! 3.  you should set the parameters rv, rh, dx to suit your application
! 
! LB. inimem.uf is read/written in real4 (hycom_1.2)
!     but not iniicemem.uf (does it make sense?)
      use mod_xc
      use mod_za
      use mod_pseudo
      implicit none 

      integer*4, external :: iargc

      character(len=80) infile,outfile,filebase,cline,rstfile,icefile, bathyfile
      character(len=6) cvarin
      character(len=7) tag7
      character(len=8) var
      character(len=3) cmem
      logical periodic, ex, lok

      real rv,rh,alp,bet, sd_d
      integer iens,i,j,ensdim, n1,n2, layer, level
      integer :: kdm
      integer :: ios,fab

      real :: hmina,hmaxa,hminb,hmaxb,rvar

      character(len=3) rungen
      character(len=50) :: tmparg
      integer k,l,iostat
      real, allocatable, dimension(:,:) :: fld2d,ranfld,accranfld, depths, wrongtopo
      real*8, allocatable, dimension(:,:) :: iofld
      real, allocatable, dimension(:,:,:) :: dptmp
      integer, allocatable :: dummy(:,:)
      real*8, dimension(:,:), allocatable ::  ficem,hicem,hsnwm,ticem,tsrfm

      real, parameter :: onem=9806.




      ! Command line argument - filename 
      if (iargc()==1) then
         call getarg(1,rstfile)
      else
         print *,'usage:'
         print *,'mkensemble restartfile'
         stop '(mkensemble)'
      end if

      fab=max(index(rstfile,'.a'),index(rstfile,'.b'))
      filebase=rstfile(1:fab-1)


      inquire(file='blkdat.input',exist=ex)
      if (.not.ex) then
         write(lp,*) 'ERROR: blkdat.input does not exist'
         stop '(forfun_nersc:blkdat.input parse)'
      endif
      open(10,file='blkdat.input')
      cvarin=''
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*)
      do while (cvarin/='idm   ')
         read(10,*) rvar,cvarin
      end do
      if (cvarin/='idm   ') then
         write(lp,*) 'ERROR: Could not get idm'
         stop '(tracerinsert:blkdat.input parse)'
      else 
         idm=nint(rvar)
      end if

      do while (cvarin/='jdm   ')
         read(10,*) rvar,cvarin
      end do
      if (cvarin/='jdm   ') then
         write(lp,*) 'ERROR: Could not get jdm'
         stop '(tracerinsert:blkdat.input parse)'
      else 
         jdm=nint(rvar)
      end if

      do while (cvarin/='kdm   ')
         read(10,*) rvar,cvarin
      end do
      if (cvarin/='kdm   ') then
         write(lp,*) 'ERROR: Could not get kdm'
         stop '(tracerinsert:blkdat.input parse)'
      else 
         kdm=nint(rvar)
      end if

      call xcspmd()
      call zaiost()
      
      allocate(fld2d    (idm,jdm))
      allocate(dummy    (idm,jdm))
      allocate(ranfld   (idm,jdm))
      allocate(accranfld(idm,jdm))
      allocate(depths   (idm,jdm))
      allocate(wrongtopo(idm,jdm))
      allocate(iofld    (idm,jdm))
      allocate(dptmp    (idm,jdm,kdm))


      ! Get depths
      if (idm>999 .or. jdm > 999) then
         write(tag7,'(i5.5,a,i5.5)')idm,'x',jdm
      else
         write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
      end if
      bathyfile='depths'//trim(tag7)//'.uf'

      ! Open file
      !print *,mnproc,' -- Read depths ..'
      open(unit=10,file=trim(bathyfile),status='old',form='unformatted')
      read(10)iofld
      close(10)
      depths=iofld







      lp=6

      print *,'test'
      print *,filebase
      print *,'test2'
      print *,rstfile
      print *,idm,jdm,kdm





      ! Generation of the initial ensemble
      write(*,'(a)',advance='no')'Choose ensemble size (sugg. 100): '
      read(*,*)ensdim
      if (ensdim <= 1) then
         print *,'ensdim must be larger than 1'
         stop
      endif

      sd_d = 0.1  ! Standard deviation of log(d) no unit
                  ! remember enstmp std Gaussian : Prob(|enstmp|>2)= 5% 
      rv=2.0e00   ! sort of vertical correlation range (in nb of layers)
      rh=25.00   ! horizontal decorrelation scale (number of grid cells)

      !  Generates the vertical correlation of the ensemble.
      alp=exp(-1.0/rv)
      print *,'Alp=',alp
      bet=sqrt(1.0-alp**2) ! keeps var(enstmp)=var(ensmem)
      print *,'Bet=',bet

      !  Prints a section showing the intial layer interfaces
      !   call intface(1,mem,onem) 
      !   call avedp(enserr,mem,onem)

      ! original Hycom state kept as first member in record 1 - 
      ! This ensures all members have the same "ident" in the restart file


      ! FFT dimensions
      n1=2**(ceiling(log(float(idm))/log(2.)))
      n2=2**(ceiling(log(float(jdm))/log(2.)))
      !call pseudo2D(ensmem,idm,jdm,nz,rh,n1,n2)

      print *
      print *
      print *,filebase
      print *,'Creating ensemble members - HYCOM:'
      do iens=2,ensdim  

         ! set file names
         write(cmem,'(i3.3)') iens
         print *,filebase
         infile=trim(filebase)
         outfile=trim(filebase)//'_mem'//cmem

         print *,'Member ',iens, 'New file : ',trim(outfile)
         print *,infile
         print *,outfile

         wrongtopo=0.

         ! Open infile
         call zaiopf(trim(infile)//'.a', 'old', 701)
         open (unit=701,file=trim(infile)//'.b',status='old')

         ! Open outfile
         call zaiopf(trim(outfile)//'.a', 'replace', 801)
         open (unit=801,file=trim(outfile)//'.b',status='replace')


         !print *
         !if (.not.read_restart_mem(rungen,rt,mem,1,1,0.)) then
         !   print *,'Read_restart failed'
         !   stop '(mkensemble)'
         !end if

         ! Start reading infile - 2 first lines are header
         read(701,'(a80)') cline
         write(801,'(a80)') cline

         read(701,'(a80)') cline
         write(801,'(a80)') cline

         ! Start reading - until we reach "dp" fields, this is 
         ! just copying fields from input to output
         do while (.not.(trim(var)=='v'.and.level==2 .and. layer==kdm))
            read(701,'(a80)') cline
            var=cline(1:8)
            i = index(cline,'=')
            read (cline(i+1:),*) layer,level,hminb,hmaxb
            !print *,var,layer,level,hminb,hmaxb
            write(801,4100) var,layer,level,hminb,hmaxb
            call zaiord(fld2d, dummy,.false., hmina,hmaxa,  701)
            call zaiowr(fld2d, dummy,.false., hmina,hmaxa,  801, .false.)
         end do

         ! Process dp layers -- dp must be the same for both time steps,
         ! otherwise temporal gradients will explode.
         do k=1,kdm 

            read(701,'(a80)') cline
            var=cline(1:8)
            i = index(cline,'=')
            read (cline(i+1:),*) layer,level,hminb,hmaxb
            !print *,var,layer,level,hminb,hmaxb
            call zaiord(fld2d, dummy,.false., hmina,hmaxa,  701)

            if (trim(var)=='dp' .and. level == 1 .and. layer==k) then


               ! 3D vertically correlated
               call pseudo2D(   ranfld,idm,jdm,1,rh,n1,n2)
               if (k==1) then
                  call pseudo2D(accranfld,idm,jdm,1,rh,n1,n2)
               else
                  accranfld=alp*accranfld+bet*ranfld
               end if


               !!ensmem=enserr*enstmp
               !!ensmem=0.10*mem%d*enstmp
               !!Lognormal law with expectation = 1
               !ensmem=mem%d*exp(enstmp*sd_d-0.5*sd_d**2)
               !enstmp=mem%d
               print *,'dp test i',k,minval(fld2d)/onem,maxval(fld2d)/onem
               fld2d=fld2d*exp(accranfld*sd_d-0.5*sd_d**2)
               print *,'dp test a',k,minval(fld2d)/onem,maxval(fld2d)/onem
               fld2d=fld2d*exp(accranfld*sd_d-0.5*sd_d**2)



               ! Accumulated depth -- needed later
               wrongtopo=wrongtopo+fld2d

               ! We can afford one 3D variable (about 60 MB for a 800x880x22 grid)
               dptmp(:,:,k)=fld2d

            end if
         end do

         print *,'topo test dpsum  i',minval(wrongtopo)/onem,maxval(wrongtopo)/onem
         print *,'topo test depths i',minval(depths),maxval(depths)

         ! layer sum may be wrong, this corrects it
         do k=1,kdm 
         do j=1,jdm
         do i=1,idm
         if (depths(i,j)>.1) then
            dptmp(i,j,k)=dptmp(i,j,k)*depths(i,j)*onem/wrongtopo(i,j)
         else
            dptmp(i,j,k)=0.
         end if
         end do
         end do
         end do

         print *,'topo test dpsum  a',minval(sum(dptmp,3))/onem,maxval(sum(dptmp,3))/onem
         print *,'topo test depths a',minval(depths),maxval(depths)

         ! Spit new layer thickness into outfile
         do k=1,kdm 
            call zaiowr(dptmp(:,:,k), dummy,.false., hmina,hmaxa,  801, .false.)
            write(801,4100) var,k,level,hmina,hmaxa
         end do

         ! Spit same layer thickness into outfile's second time level
         do k=1,kdm 
            ! Need to skip forward in infile as well
            read(701,'(a80)') cline
            var=cline(1:8)
            i = index(cline,'=')
            read (cline(i+1:),*) layer,level,hminb,hmaxb
            call zaiord(fld2d, dummy,.false., hmina,hmaxa,  701)

            call zaiowr(dptmp(:,:,k), dummy,.false., hmina,hmaxa,  801, .false.)
            write(801,4100) var,layer,level,hmina,hmaxa
         end do




         ! The rest of the outfile is the same as the infile
         ios = 0
         do while (ios==0)

            read(701,'(a80)',iostat=ios) cline

            if (ios==0) then
               var=cline(1:8)
               i = index(cline,'=')
               read (cline(i+1:),*) layer,level,hminb,hmaxb
               !print *,var,layer,level,hminb,hmaxb
               write(801,4100) var,layer,level,hminb,hmaxb
               call zaiord(fld2d, dummy,.false., hmina,hmaxa,  701)
               call zaiowr(fld2d, dummy,.false., hmina,hmaxa,  801, .false.)
            end if
         end do

         ! Close up -- go home 
         call zaiocl(701)
         close(701)
         call zaiocl(801)
         close(801)

         
      enddo


! ICE ensemble
!      print *
!      print *
!      icefile=trim(filebase)//'ICE.uf'
!      inquire(exist=ex, file=trim(icefile))
!      if (.not. ex) then 
!         print *,'No ice file present'
!         stop
!      end if
!      allocate(ficem(idm,jdm))
!      allocate(hicem(idm,jdm))
!      allocate(hsnwm(idm,jdm))
!      allocate(ticem(idm,jdm))
!      allocate(tsrfm(idm,jdm))
!      inquire(iolength=j) ficem,hicem,hsnwm,ticem,tsrfm  
!      open(10,file=trim(icefile),status='old',form='unformatted', &
!              access='direct',recl=j)
!
!
!      print *,'Creating ensemble members - ICE:'
!      do iens=2,ensdim
!         read(10 ,rec=1   ,iostat=ios) ficem,hicem,hsnwm,ticem,tsrfm
!         write(10,rec=iens,iostat=ios) ficem,hicem,hsnwm,ticem,tsrfm
!      enddo

      close(10)

 4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)


      end program mkensemble


