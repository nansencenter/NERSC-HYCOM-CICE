! KAL -- This routine interpolates horizontal mercator fields to a NERSC model
! KAL 
! KAL -- Usable/Inital version: 31.01.2007 - Knut Lisæter



      program vremap_mercator
      use mod_xc
      use mod_za
      use mod_sigma
      use mod_grid, only: hycbathy => depths, hyclon => plon , hyclat => plat, &
                          get_grid
      use mod_hycomfile_io
      use m_parse_blkdat
      use m_layer_remapV2
      use m_layer_mixV1
      use m_pbavg_from_ssh
      use m_rotate2
      implicit none

      real, dimension(:,:), allocatable ::  &
         rtmp,dens, dp, &
         ubavg, vbavg, ssh, psikk, thkk, pbavg, intup

      real, dimension(:,:,:), allocatable ::  &
         oldtemp, oldsaln, newint, newtemp, newsaln, &
         oldutot, oldvtot, newutot, newvtot 

      real, allocatable :: newdens(:), oldint(:), dp0k(:), oldint2(:)

      integer, dimension(:,:), allocatable ::  &
         ip, maxk
      integer :: thflag, kapflg 
      integer :: newkdm,oldkdm, oldkdm2

      character(len=8) :: c8
      character(len=80) :: template
      integer :: k,idummy,ios,kold,i,j,fnd
      integer :: ios2, itime
      real    :: amin, amax,rdummy, bmin, bmax, zlev, thbase
      logical :: lmaxk
#if defined (IARGC)
      integer*4, external :: iargc
#endif
      type(hycomfile) :: hfile

      if (iargc()==1) then
         call getarg(1,template)
         fnd=max(index(template,'.a'),index(template,'.b'))
         if (fnd<=0) then
            print *,'error in restart template (input)'
            call exit(1)
         end if
         template=template(1:fnd-1)
      else
         print *,'usage: vremap_mercator restart_template'
         call exit(1)
      endif


      ! Get HYCOM lon lat and depths
      ! Initialize Arrai IO
      call xcspmd()
      call zaiost()
      call get_grid()

      allocate(ubavg(idm,jdm))
      allocate(vbavg(idm,jdm))
      allocate(ssh(idm,jdm))
      allocate(intup(idm,jdm))
      allocate(psikk(idm,jdm))
      allocate(thkk(idm,jdm))
      allocate(pbavg(idm,jdm))
      allocate(ip    (idm,jdm))
      allocate(rtmp  (idm,jdm))
      print *,'hyc lon:',minval(hyclon),maxval(hyclon)
      print *,'hyc lat:',minval(hyclat),maxval(hyclat)


      ! Get thflag,kapflg,kdm from blkdat file
      call parse_blkdat('kdm   ','integer',rdummy,newkdm)
      call parse_blkdat('thflag','integer',rdummy,thflag)
      call parse_blkdat('kapflg','integer',rdummy,kapflg)
      call parse_blkdat('thbase','real',thbase,idummy)

      ! Get densities from blkdat file
      allocate(newdens(newkdm))
      do k=1,newkdm
         call parse_blkdat('sigma ','real',newdens(k),idummy,imatch=k)
      end do


      ! Go through horizontally interpolated mercator  files -- 
      ! get oldkdm from highest temperature level entry in .b file
      open(10,file='merchint.b')
      oldkdm=1
      ios2=0
      do while (ios2==0)
         read(10,104,iostat=ios2) c8,k,zlev,bmin,bmax
         if (trim(c8)=='temp'.and.ios2==0) oldkdm=max(oldkdm,k)
      end do
      close(10)
      print *,'oldkdm=',oldkdm

      ! Read mercator 3D fields - only temp/saln now
      allocate(oldtemp(idm,jdm,oldkdm))
      allocate(oldsaln(idm,jdm,oldkdm))
      allocate(oldutot(idm,jdm,oldkdm))
      allocate(oldvtot(idm,jdm,oldkdm))
      allocate(maxk   (idm,jdm))
      allocate(oldint (oldkdm))
      allocate(oldint2(oldkdm))
      open(99,file='merchint.b')
      call zaiopf('merchint.a','old',99)
      maxk=0
      ios=0
      lmaxk=.false.
      do while (ios==0)

         read(99,104,iostat=ios) c8,k,zlev,bmin,bmax
         if (ios==0) then
            !print *,'k is ',k
            if(k>0) oldint(k)=zlev
            call zaiord(rtmp,ip,.false.,amin,amax,99) !plon
            if (trim(c8)=='temp') then
               !print *,'temp match',k
               oldtemp(:,:,k)=rtmp
            elseif (trim(c8)=='saln') then
               !print *,'saln match',k
               oldsaln(:,:,k)=rtmp
            elseif (trim(c8)=='utot') then
               !print *,'uvel match',k
               oldutot(:,:,k)=rtmp
            elseif (trim(c8)=='vtot') then
               !print *,'vvel match',k
               oldvtot(:,:,k)=rtmp
            elseif (trim(c8)=='ubavg') then
               !print *,'vbavg match',k
               ubavg=rtmp
            elseif (trim(c8)=='vbavg') then
               !print *,'vbavg match',k
               vbavg=rtmp
            elseif (trim(c8)=='ssh') then
               !print *,'ssh match',k
               ssh=rtmp
            elseif (trim(c8)=='maxk') then
               lmaxk=.true.
               maxk=nint(rtmp)
               !print *,'got maxk'
            end if
         end if
      end do
      !print *,'end'
      close(99)
      call zaiocl(99)
      !stop


      ! 
      if (.not. lmaxk) then
         print *,'Can not get maxk from files'
         print *, '(vremap_mercator)'
         call exit(1)
      end if


      !Set up new layers
      allocate(newint (idm,jdm,newkdm))
      allocate(newtemp(idm,jdm,newkdm))
      allocate(newsaln(idm,jdm,newkdm))
      allocate(newutot(idm,jdm,newkdm))
      allocate(newvtot(idm,jdm,newkdm))
      allocate(dens   (idm,jdm))
      allocate(dp     (idm,jdm))

      ! Cycle through-create new interfaces
      print *,'WARN WARN WARN'
      print *,'WARN WARN WARN'
      print *,'WARN WARN WARN'
      print *,'WARN WARN WARN'
      print *,'WARN WARN WARN'
      print *,'dp0k set specifically in p_vremap_mercator.F90'
      allocate(dp0k(newkdm))
      dp0k=10.
      newint=0.
      do j=1,jdm
      do i=1,idm
         if (maxk(i,j)>0) then

            !print *,i,j 
            oldkdm2=maxk(i,j)
            ! TODO - calc newdp0 !

            call layer_remapV2(oldint(1:oldkdm2),oldtemp(i,j,1:oldkdm2), &
                               oldsaln(i,j,1:oldkdm2),oldkdm2, &
                                newdens,newint(i,j,:),dp0k,newkdm,thflag,.false.)

            where (newint(i,j,:)>hycbathy(i,j)) newint(i,j,:)=hycbathy(i,j)
            where (newint(i,j,:)<hycbathy(i,j) .and.  &
                   abs(newint(i,j,:)-oldint(oldkdm2))<1e-4) newint(i,j,:)=hycbathy(i,j)


            ! "Correct old interfaces" -- needed in procedure below
            
         end if
      end do
      end do


      



          
      ! The vertical interpolation is set up, now mix through layers to get new
      ! temperature and salinities
      newsaln=0.
      newtemp=0.
      do j=1,jdm
      do i=1,idm
         if (maxk(i,j)>0) then
            oldkdm2=maxk(i,j)

            ! Temporary - adheres to hycom bathymetry
            oldint2=oldint
            oldint2=min(oldint2,hycbathy(i,j))
            if (oldkdm2<oldkdm) then
               oldint2(oldkdm2:oldkdm)=hycbathy(i,j)
            else
               oldint2(oldkdm)=hycbathy(i,j)
            end if

            !print *,'hycbathy ,oldkdm: ',hycbathy(i,j),oldkdm,oldkdm2

            call layer_mixV1(oldint2(1:oldkdm2),oldtemp(i,j,1:oldkdm2),oldkdm2, &
                             newint(i,j,:),newtemp(i,j,:),newkdm)

            call layer_mixV1(oldint2(1:oldkdm2),oldsaln(i,j,1:oldkdm2),oldkdm2, &
                             newint(i,j,:),newsaln(i,j,:),newkdm)

            call layer_mixV1(oldint2(1:oldkdm2),oldutot(i,j,1:oldkdm2),oldkdm2, &
                             newint(i,j,:),newutot(i,j,:),newkdm)

            call layer_mixV1(oldint2(1:oldkdm2),oldvtot(i,j,1:oldkdm2),oldkdm2, &
                             newint(i,j,:),newvtot(i,j,:),newkdm)

            !print *,'old int:',oldint(:)
            !print *,'new int:',newint(i,j,:)
            !print *,'old temp:',oldtemp(i,j,:)
            !print *,'new temp:',newtemp(i,j,:)
            !print *


         end if
      end do
      end do

      ! Calculate interface values to dp (needed later
      do k=newkdm,2,-1 ! k=1 is already ok

         newint(:,:,k)=newint(:,:,k)-newint(:,:,k-1)

      end do

      newint=newint*onem !!



      ! Calculate barotropic pressures from ssh - read some values from a rst
      ! file (template)
      call initHF(hfile,trim(template),'restart')
      call HFreadfield(hfile,psikk,idm,jdm,'psikk   ',0,1)
      call HFreadfield(hfile,thkk ,idm,jdm,'thkk    ',0,1)
      call pbavg_from_ssh(newsaln,newtemp,newint,ssh,psikk,thkk,hycbathy,pbavg,idm,jdm,newkdm)

      ! Rotate velocities from e/n to grid as well
      print *,'Rotating velocities'
      do k=1,newkdm
         call rotate2(newutot(:,:,k),newvtot(:,:,k),hyclat,hyclon,idm,jdm,'l2m')
      end do
      call rotate2(ubavg,vbavg,hyclat,hyclon,idm,jdm,'l2m')


      

      call zaiopf('mercvint.a','replace',99)
      open (98,file='mercvint.b',status='replace')


      ! Dump to hycom type restart file
      do itime=1,2
      do k=1,newkdm
         newutot(:,:,k) = newutot(:,:,k)-ubavg ! Baroclinic
         call zaiowr(newutot(:,:,k),ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'u       ',k,itime,amin,amax
      end do
      end do

      do itime=1,2
      do k=1,newkdm
         newvtot(:,:,k) = newvtot(:,:,k)-vbavg ! Baroclinic
         call zaiowr(newvtot(:,:,k),ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'v       ',k,itime,amin,amax
      end do
      end do

      do itime=1,2
      do k=1,newkdm
         !if (k==1) then
         !   dp=newint(:,:,k)
         !else
         !   dp=newint(:,:,k) - newint(:,:,k-1)
         !end if
         call zaiowr(newint(:,:,k) ,ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'dp      ',k,itime,amin,amax
      end do
      end do

      do itime=1,2
      do k=1,newkdm
         call zaiowr(newtemp(:,:,k),ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'temp    ',k,itime,amin,amax
      end do
      end do

      do itime=1,2
      do k=1,newkdm
         call zaiowr(newsaln(:,:,k),ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'saln    ',k,itime,amin,amax
      end do
      end do

      do itime=1,2
      do k=1,newkdm
         do j=1,jdm
         do i=1,idm
            
            if (thflag==2) then
               dens(i,j)=sig2(newtemp(i,j,k),newsaln(i,j,k))
            elseif (thflag==0) then
               dens(i,j)=sig0(newtemp(i,j,k),newsaln(i,j,k))
            else
               print *,'Unknown thflag ',thflag
               call exit(1)
            end if
         end do
         end do
         dens=dens-thbase
         call zaiowr(dens(:,:),ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'th3d    ',k,itime,amin,amax
      end do
      end do

      do itime=1,3
         call zaiowr(ubavg(:,:),ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'ubavg   ',0,itime,amin,amax
      end do

      do itime=1,3
         call zaiowr(vbavg(:,:),ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'vbavg   ',0,itime,amin,amax
      end do

      do itime=1,3
         call zaiowr(pbavg(:,:),ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'pbavg   ',0,itime,amin,amax
      end do


      rtmp=hycbathy*onem
      call zaiowr(rtmp,ip,.false.,amin,amax,99,.false.) !plon
      write(98,4100) 'pbot    ',0,1,amin,amax

      call zaiowr(psikk,ip,.false.,amin,amax,99,.false.) !plon
      write(98,4100) 'psikk   ',0,1,amin,amax

      call zaiowr(thkk,ip,.false.,amin,amax,99,.false.) !plon
      write(98,4100) 'thkk    ',0,1,amin,amax

      ! KAL - final fields in a restart file are dpmixl, we set them to
      ! zero here
      rtmp=0.
      do itime=1,2
         call zaiowr(rtmp,ip,.false.,amin,amax,99,.false.) !plon
         write(98,4100) 'dpmixl  ',0,itime,amin,amax
      end do

      close(98)
      call zaiocl(99)




         

      



   104 format(a8," k=",i3," level=",f10.2,"  "," min/max:",2e14.6)

     ! KAL -- Format used in the final vertically interpolated file (readable by
     ! KAL -- restart routines
4100   format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)

      end program vremap_mercator





         


