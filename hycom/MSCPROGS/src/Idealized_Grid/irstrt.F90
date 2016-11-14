      program irstrt
      use mod_parameters
      use mod_xc 
      use mod_za
      use mod_sigma
      use m_parse_blkdat
      ! Use this program to create  a restart file corresponding
      ! to a user-specified vertical setup. Layer densities are
      ! read from blkdat.input.

      ! There are some examples of how to set up temp/saln and th3d below.
      ! No input arguments to this routine, if you use it, you should
      ! set up the vertical profile yourself. This program is only meant as a 
      ! starting point.
      implicit none


      character(len=80) gridid, cline, flnmdep, flnmrso, icename
      character(len=11) tag7

      real, allocatable, dimension(:,:,:) ::  saln, th3d, dp, temp,p
      real, allocatable, dimension(:,:)   ::  depths, psikk, thkk, tmp, pbot
      real, allocatable, dimension(:)     ::  theta

      integer, allocatable, dimension(:,:) ::  ip

      real :: centerlon, centerlat, dx, dy
      real :: rad, deg
      real :: xmin, xmax, rdummy
      integer :: idummy
      integer :: i,j,k,l,ios

      integer :: kapflg, thflag
      real    :: thbase, dtime
      integer :: kdm
      integer :: yrflag, iexpt, iversn, nstep=0
      integer :: idmtst, jdmtst

      real*8, allocatable, dimension(:,:)  ::  & ! Hard-coded real*8 for file IO
        io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm
        logical :: case_1, case_2


      flnmdep='regional.depth'
      flnmrso='restart.out'
      dtime=0.0
      iexpt=1
      iversn=21
      nstep=0




      ! First tests -- hardcoded grid size etc
      call parse_blkdat('idm   ','integer',rdummy,idmtst)
      call parse_blkdat('jdm   ','integer',rdummy,jdmtst)
      call xcspmd()
      call zaiost()

      if (idm/=idmtst .or. jdm/=jdmtst) then
         print *,'Grid size inconsistency - i: ',idm,idmtst
         print *,'Grid size inconsistency - j: ',jdm,jdmtst
         print *,'check that blkdat.input and regional.grid.b match'
         stop '()'
      end if


      call parse_blkdat('kapflg','integer',rdummy,kapflg)
      call parse_blkdat('thflag','integer',rdummy,thflag)
      call parse_blkdat('kdm   ','integer',rdummy,kdm   )
      call parse_blkdat('yrflag','integer',rdummy,yrflag)
      call parse_blkdat('thbase','real',thbase,idummy)

      ! First tests - very simple stuff
      allocate(temp  (idm,jdm,kdm))
      allocate(saln  (idm,jdm,kdm))
      allocate(th3d  (idm,jdm,kdm))
      allocate(theta (kdm))
      allocate(tmp   (idm,jdm))
      allocate(dp    (idm,jdm,kdm))
      allocate(p     (idm,jdm,kdm+1))
      allocate(psikk (idm,jdm))
      allocate(thkk  (idm,jdm))
      allocate(pbot  (idm,jdm))
      allocate(depths(idm,jdm))
      allocate(ip    (idm,jdm))
      allocate(io_ficem(idm,jdm))
      allocate(io_hicem(idm,jdm))
      allocate(io_hsnwm(idm,jdm))
      allocate(io_ticem(idm,jdm))
      allocate(io_tsrfm(idm,jdm))
      ip=0

      print *,'Reading densities'
      do k=1,kdm
         call parse_blkdat('sigma ','real',theta(k),idummy,imatch=k)
         !print *,theta(k)
      end do
         

      ! Get model bathymetry
      call zaiopf(flnmdep(1:len_trim(flnmdep))//'.a','old', 12)
      call zaiord(depths, ip,.false., xmin,xmax,  12)
      call zaiocl(12)
      where(depths>0.5*huge) depths=0.


   case_1=.false.
   case_2=.true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  --Case 1 : Wedge of low-salinity water at freezing point !!!!!!
!!!             From surface to  50 in "north"                !!!!!!
!!!             The remainder is 0 degree water               !!!!!!
!!!             Required - two layers                         !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (case_1) then

         if (kdm/=2) then
            print *,'Requires 2 layers'
            stop
         end if

         ! Constant salinity
         saln(:,  1:100,1)=35.
         saln(:,101:200,1)=34.0
         saln(:,   :   ,2)=35.

         ! Top layer freezing point in "north"
         do j=1,jdm
         do i=1,idm
            temp(i,j,1)=tofsig0(theta(1),saln(i,j,1))
         end do
            !print *,theta(1),saln(50,j,1),temp(50,j,1)
         end do
         !print *,minval(temp(:,:,1))
         !print *,maxval(temp(:,:,1))

         do j=1,jdm
         do i=1,idm
            temp(i,j,2)=tofsig0(theta(2),saln(i,j,2))
         end do
         end do

         ! Layer thicknesses
         dp(:,  1:95,2)=(depths(:,1:95)-10.)*onem
         dp(:,105:200,2)=(depths(:,105:200)-100.)*onem
         do j=96,104
         do i=1,idm
            dp(i,j,2)=dp(i,70 ,2) * real(105-  j)/10. + &
                      dp(i,130,2) * real(j  - 95)/10.
         end do
         end do
         dp(:,   :   ,1)=max(0.,depths*onem-dp(:,:,2));

         ! ice thickness and concentration
         io_hicem(:,  1:100)=0. ; io_ficem(:,1  :100)=0.
         io_hicem(:,101:200)=1.5; io_ficem(:,101:200)=1.
         io_hsnwm=0.
         io_tsrfm=273.
         io_ticem=273.

         



         do k=1,kdm
         do j=1,jdm
         do i=1,idm
            if (thflag==0) then
               th3d(i,j,k)=sig0(temp(i,j,k),saln(i,j,k)) - thbase
            elseif (thflag==2) then
               th3d(i,j,k)=sig2(temp(i,j,k),saln(i,j,k)) - thbase
            elseif (thflag==4) then
               th3d(i,j,k)=sig4(temp(i,j,k),saln(i,j,k)) - thbase
            else
               print *,'(unknown thflag)'
               stop '(irstrt)'
            end if
         end do
         end do
         end do

      elseif (case_2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  --Case 2 : Seamount case - kdm layers, salinity = 35 psu !!!!!!
!!!             Densities read from blkdat.input. Temperature !!!!!!
!!!             diagnosed. Layer thicknesses are the same     !!!!!!
!!!             initially, so stratification is determined    !!!!!!
!!!             from blkdat.input!!!!!                        !!!!!!
!!!             No ice in this case                           !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         saln=35.

         ! Constant salinity of 35
         !do k=1,kdm
         !saln(:,:,k)=35.
         !end do
         saln(:,:,:)=35.


         do k=1,kdm
                 do j=1,jdm
                 do i=1,idm
                    temp(i,j,k)=tofsig0(theta(k),saln(i,j,k))
                 end do
                 end do
         end do


         !!! SETTING UP FLAT LAYERS FOR SEAMOUNT TEST!!!
         ! Layer thicknesses in order to make flat layers (Honour to KAL!!!)
         do k=1,kdm
         dp(:,:,k)=(maxval(depths(:,:)*onem)/kdm)
         end do

         ! Setting up pressure surfaces at layer interfaces
          p=0.
          do k=1,kdm
             p(:,:,k+1)=p(:,:,k)+dp(:,:,k)
          end do

         ! Correcting pressure surfaces hindering them to extend deeper than bottom
          do k=1,kdm 
             p(:,:,k+1)=min(depths(:,:)*onem,p(:,:,k+1))
          end do


         ! Calculating new dp's using corrected pressure surfaces
          do k=1,kdm
             dp(:,:,k)=p(:,:,k+1)-p(:,:,k)
          end do


   !KAL      ! Specifying density using function for density as
   !KAL      ! a function of temperature (S=35)
          do k=1,kdm
          do j=1,jdm
          do i=1,idm
          if (thflag==0) then
                th3d(i,j,k)=sig0(temp(i,j,k),saln(i,j,k)) - thbase
             elseif (thflag==2) then
                th3d(i,j,k)=sig2(temp(i,j,k),saln(i,j,k)) - thbase
             elseif (thflag==4) then
                th3d(i,j,k)=sig4(temp(i,j,k),saln(i,j,k)) - thbase
             else
                print *,'(unknown thflag)'
                stop '(irstrt)'
             end if
          end do
          end do
          end do
       end if




      ! Create 
      write(lp,'(a)',advance='yes')'Saving : '//trim(flnmrso)//'.[ab]'
      call zaiopf(flnmrso(1:len_trim(flnmrso))//'.a','replace', 12)
      open (unit=12,file=flnmrso(1:len_trim(flnmrso))//'.b', &
            status='replace',action='write',form='formatted')


      ! Header
      write(12,'(a,3i6)') 'RESTART: iexpt,iversn,yrflag = ', &
                                    iexpt,iversn,yrflag      
      write(cline,*)                nstep,dtime
      write(12,'(a,a)')   'RESTART: nstep,dtime = ', &
                                    cline(1:len_trim(cline))


      ! u - set to zero
      do l=1,2
      do k=1,kdm
         ! Velocities
         call random_number(tmp)
         tmp=(tmp-0.5)*.001
         !tmp=0.
         call zaiowr(tmp, ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'u       ',k,l,xmin,xmax
      end do
      end do

      ! v - set to zero
      do l=1,2
      do k=1,kdm
         ! Velocities
         call random_number(tmp)
         tmp=(tmp-0.5)*.001
         !tmp=0.
         call zaiowr(tmp, ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'v       ',k,l,xmin,xmax
      end do
      end do

      ! dp - set to fixed-levels initially
      do l=1,2
      do k=1,kdm
         call zaiowr(dp(:,:,k) , ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'dp      ',k,l,xmin,xmax
      end do
      end do

      ! temp 
      do l=1,2
      do k=1,kdm
         call zaiowr(temp(:,:,k) , ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'temp    ',k,l,xmin,xmax
      end do
      end do

      ! saln 
      do l=1,2
      do k=1,kdm
         call zaiowr(saln(:,:,k) , ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'saln    ',k,l,xmin,xmax
      end do
      end do

      ! dens
      do l=1,2
      do k=1,kdm
         call zaiowr(th3d(:,:,k) , ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'th3d    ',k,l,xmin,xmax
      end do
      end do

      ! ubavg
      tmp=0.
      k=0
      do l=1,3
         call zaiowr(tmp , ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'ubavg   ',k,l,xmin,xmax
      end do

      ! vbavg
      tmp=0.
      k=0
      do l=1,3
         call zaiowr(tmp , ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'vbavg   ',k,l,xmin,xmax
      end do

      ! pbavg
      tmp=0.
      k=0
      do l=1,3
         call zaiowr(tmp , ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'pbavg   ',k,l,xmin,xmax
      end do

      ! pbot 
      pbot=depths*onem
      k=0
      l=1
      call zaiowr(pbot , ip,.false., xmin,xmax,  12, .true.)
      write(12,4100) 'pbot    ',k,l,xmin,xmax

      ! psikk
      psikk=0.
      k=0
      l=1
      psikk=0.
      call zaiowr(psikk , ip,.false., xmin,xmax,  12, .true.)
      write(12,4100) 'psikk   ',k,l,xmin,xmax

      ! thkk
      k=0
      l=1
      thkk=th3d(:,:,kdm)
      call zaiowr(thkk , ip,.false., xmin,xmax,  12, .true.)
      write(12,4100) 'thkk    ',k,l,xmin,xmax

      ! dpmixl
      k=0
      do l=1,2
         call zaiowr(dp(:,:,1) , ip,.false., xmin,xmax,  12, .true.)
         write(12,4100) 'dpmixl  ',k,l,xmin,xmax
      end do

 4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)



      call ZAIOCL(12)
      close(12)





      ! Create a ice restart file as well
      icename='fileICE.uf'
      write(lp,'(a)',advance='yes')'Saving : '//trim(icename)
      inquire(iolength=j) io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm  
      open(10,file=icename,status='unknown',form='unformatted', &
           access='direct',recl=j)
      write(10,rec=1,err=200,iostat=ios)io_ficem,io_hicem, &
           io_hsnwm,io_ticem,io_tsrfm  
200   continue
      close(10)


      end
