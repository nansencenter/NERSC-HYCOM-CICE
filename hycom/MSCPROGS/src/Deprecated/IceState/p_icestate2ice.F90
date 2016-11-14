   program icestate2ice
      use mod_year_info
      use mod_xc
      use mod_za
      implicit none

      ! Program reads icestate files, dumps to ICE format

      character(len=3) :: rungen
      type(year_info)  :: rt
      integer          :: imem  ! ensemble inices
      
      character*80 :: icename,oldicename
      integer      :: j, ios, i,l, irec
      logical      :: ex
      integer :: hl,hk, ifld
      real :: xmin,xmax,xmin2,xmax2

      real*8, dimension(:,:), allocatable :: ficem, hicem, tsrfm,ticem, &
         hsnwm
      real, dimension(:,:), allocatable :: tmpfld
      integer, dimension(:,:), allocatable :: ip
      character*8  :: char8
      integer :: kdm, nthick, nlaymax, iyy,idd, ihh
      integer :: find, nflds 
      integer :: icestate_unit=97
      

      if (iargc()/=1) then
         print *,'Usage: icestate2ice <icestate-file>'
         stop
      else
         call getarg(1,icename)
         find=max(index(icename,'.a'),index(icename,'.b'))
         icename=icename(1:find-1)
         find=index(icename,'ICESTATE')
         oldicename=icename(1:find-1)//'ICE.uf'
      end if

      print *,icename
      print *,oldicename
      imem = 1 ! For now


      inquire(file=trim(icename)//'.b',exist=ex)
      if (.not.ex) then
         if (mnproc==1) write(lp,*) 'File ',trim(icename)//'.[ab] does not exist...'
         stop
      end if

      ! Get header etc etc etc
      open(icestate_unit,file=trim(icename)//'.b',status='old')
      read(icestate_unit,6000) idm
      read(icestate_unit,6000) jdm
      read(icestate_unit,6000) kdm
      read(icestate_unit,6000) nthick
      read(icestate_unit,6000) nlaymax
      read(icestate_unit,6000) iyy
      read(icestate_unit,6000) idd
      read(icestate_unit,6000) ihh
      read (icestate_unit,'(a80)')  
!6000 format(a6,' =',i6)
 6000 format(i6)


      print *,nthick
      print *,nlaymax
      if (nthick>1 .or. nlaymax>1 ) then
         print *,'Routine must be modified to handle nthick>1, nlay>1'
         stop
      end if

      call xcspmd()
      call zaiost()

      ! Start allocating arrays, now that we have dimensions....
      allocate(ficem(idm,jdm))
      allocate(hicem(idm,jdm))
      allocate(hsnwm(idm,jdm))
      allocate(tsrfm(idm,jdm))
      allocate(ticem(idm,jdm))
      allocate(tmpfld(idm,jdm))
      allocate(ip    (idm,jdm))

      call zaiopf(trim(icename)//'.a','old',icestate_unit)

      ! Skip forward  for ensemble runs
      if (imem>1) then 
         do ifld=1,(imem-1)*nflds
            call zaiosk(icestate_unit)
         end do
      end if

      nflds=nlaymax*nthick+8*nthick+3

      ! Read ice variables
      ios=0
      irec=1
      do while (ios==0 .and. irec <= nflds ) 
         read(icestate_unit,117) char8,hk,hl,xmin2,xmax2
         print *,char8,hk,hl

         if (trim(char8)=='nlay') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            !!print *,maxval(tmpfld),minval(tmpfld)
            !icestate%ice(hk)%nlay=1
            !do j=1-margin,jj+margin
            !do l=1,isp(j)
            !do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            !   !print *,i,j,hk,tmpfld(i,j)
            !   icestate(i,j)%ice(hk)%nlay=nint(tmpfld(i,j))
            !end do
            !end do
            !end do
         elseif (trim(char8)=='qstore') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> Not needed for old ice model'
         elseif (trim(char8)=='albs') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> Not needed for old ice model'
         elseif (trim(char8)=='fice') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> ficem needed for old ice model'
            ficem=tmpfld
            print *,maxval(ficem)
         elseif (trim(char8)=='hice') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> hicem needed for old ice model'
            hicem=tmpfld
            print *,maxval(hicem)
         elseif (trim(char8)=='hsnw') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> hsnwm needed for old ice model'
            hsnwm=tmpfld
            print *,maxval(hsnwm)
         elseif (trim(char8)=='rhosnw') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> Not needed for old ice model'
         elseif (trim(char8)=='tsrf') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> tsrfm needed for old ice model (well, not really but)'
            tsrfm=tmpfld
            ticem=tmpfld
            print *,maxval(tsrfm)
         elseif (trim(char8)=='vtp') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> vtp needed for old ice model'
            ticem=tmpfld
            print *,maxval(ticem)
         elseif (trim(char8)=='hml') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> Not needed for old ice model'
         elseif (trim(char8)=='sml') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> Not needed for old ice model'
         elseif (trim(char8)=='tml') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            print *,'   |-----> Not needed for old ice model'
         end if

         if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.  abs(xmax-xmax2)>abs(xmax)*1e-4) then
            if     (mnproc.eq.1) then
               write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
                 'error - .a and .b files not consistent:', &
                 '.a,.b min = ',xmin,xmin2 ,xmin2-xmin , &
                 '.a,.b max = ',xmax,xmax2 ,xmax2-xmax 
               print *,'Variable :',char8,hk,hl
            endif
            call xcstop('(read_ave2)')
            stop '(read_ave2)'
         end if

         ireC=irec+1
      end do

      call zaiocl(icestate_unit)
      if (mnproc==1) close(icestate_unit)


      ! Dump to old ice file
      inquire(iolength=j) ficem,hicem,hsnwm,ticem,tsrfm  !,ticeU,ticeV
      open(87,file=trim(oldicename),status='unknown',form='unformatted', access='direct',recl=j)
      write(87,rec=imem)ficem,hicem, hsnwm,ticem,tsrfm   !,io_iceU,io_iceV
      close(87)




 117  format (a8,' = ',i3,i3,1p2e16.7)

   end program



