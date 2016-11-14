     module mod_icestate_io
     private :: icestate_file
     integer, parameter :: icestate_unit=31 ! Unit identifier
     integer, private :: nflds
     contains





! ===================================================================================
! ============================= save_restart_icestate ===============================
! ===================================================================================
   logical function save_restart_icestate(rungen,rt,imem)
      use mod_icestate_grid
      use mod_year_info
      use mod_icestate
      use mod_za
      implicit none

      character(len=3),intent(in) :: rungen
      type(year_info), intent(in) :: rt
      integer,         intent(in) :: imem  ! ensemble indice

      character*40 :: icename
      character*8  :: char8
      integer      :: j, ios,ifld,hk,hl
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmpfld

      real :: xmin,xmax,xmin2,xmax2


      nflds=nlaymax*nthick+8*nthick+3

      icename=ICESTATE_file(rungen,rt,'restart')
      if (mnproc==1) then
         write(lp,*) 'saving ICESTATE:',icename
         open(icestate_unit,file=trim(icename)//'.b',status='unknown')
         write(icestate_unit,116) itdm,jtdm, kdm,nthick,nlaymax,rt%iyy,rt%idd,rt%ihh
         write (icestate_unit,'(a)') 'char8,hk,hl,xmin,xmax'
      end if

      call zaiopf(trim(icename)//'.a','unknown',icestate_unit)


      ! Skip forward  for ensemble runs
      if (imem>1) then 
         do ifld=1,(imem-1)*nflds
            call zaiosk(icestate_unit)
         end do
      end if

      ! Write ice variables
      do hk=1,nthick
      do ifld = 1,8+nlaymax
         if (    ifld==1) then
            tmpfld=float(icestate%ice(hk)%nlay)
            call zaiowr(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='nlay    ' ; hl=1
         elseif (ifld==2) then
            call zaiowr(icestate%ice(hk)%qstore,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='qstore  ' ; hl=1
         elseif (ifld==3) then
            call zaiowr(icestate%ice(hk)%albs,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='albs    ' ; hl=1
         elseif (ifld==4) then
            call zaiowr(icestate%ice(hk)%fice,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='fice    ' ; hl=1
         elseif (ifld==5) then
            call zaiowr(icestate%ice(hk)%hice,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='hice    ' ; hl=1
         elseif (ifld==6) then
            call zaiowr(icestate%ice(hk)%hsnw,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='hsnw    ' ; hl=1
         elseif (ifld==7) then
            call zaiowr(icestate%ice(hk)%rhosnw,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='rhosnw  ' ; hl=1
         elseif (ifld==8) then
            call zaiowr(icestate%ice(hk)%tsrf,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='tsrf    ' ; hl=1
         elseif (ifld>8) then
            char8='vtp     ' ; hl=ifld-8
            call zaiowr(icestate%ice(hk)%vtp(hl),ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
         end if
         if     (mnproc.eq.1) then
            write (icestate_unit,117) char8,hk,hl,xmin,xmax
            call flush(icestate_unit)
         endif !1st tile
      end do
      end do

      hk=0 ; hl=0
      do ifld = 1,3
         if (    ifld==1) then
            call zaiowr(icestate%hml,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='hml     ' ; hl=1
         elseif (ifld==2) then
            call zaiowr(icestate%sml,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='sml     ' ; hl=1
         elseif (ifld==3) then
            call zaiowr(icestate%tml,ip,.false., &
                        xmin,xmax, icestate_unit, .false.)
            char8='tml     ' ; hl=1
         end if
         if     (mnproc.eq.1) then
            write (icestate_unit,117) char8,hk,hl,xmin,xmax
            call flush(icestate_unit)
         endif !1st tile
      end do

      call zaiocl(icestate_unit)
      if (mnproc==1) close(icestate_unit)

      save_restart_icestate=.true.


 116  format (                                        &
       i5,4x,'''idm   '' = longitudinal array size' / &
       i5,4x,'''jdm   '' = latitudinal  array size' / &
       i5,4x,'''kdm   '' = Vertical     array size' / &
       i5,4x,'''nthick'' = Num ice categories     ' / &
       i5,4x,'''nlay  '' = Max num ice layers     ' / &
       i5,4x,'''dyear '' = Year of integration dump'/ &
       i5,4x,'''dday  '' = Day of integration dump' / &
       i5,4x,'''dhour '' = Hour of integration dump')
 117  format (a8,' = ',i3,i3,1p2e16.7)

   end function save_restart_icestate 



! ===================================================================================
! ============================= read_restart_icestate ===============================
! ===================================================================================
   logical function read_restart_icestate(rungen,rt,imem)
      use mod_year_info
      use mod_icestate
      use mod_za
      implicit none

      character(len=3),intent(in) :: rungen
      type(year_info), intent(in) :: rt
      integer,         intent(in) :: imem  ! ensemble inices
      
      character*40 :: icename
      integer      :: j, ios, i,l, irec
      logical      :: ex
      integer :: hl,hk, ifld
      real :: xmin,xmax,xmin2,xmax2

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmpfld
      character*8  :: char8
      integer :: l_itdm, l_jtdm, &
                 l_kdm, l_nthick, l_nlaymax, l_iyy,l_idd, l_ihh
      
      nflds=nlaymax*nthick+8*nthick+3

      icename=ICESTATE_file(rungen,rt,'restart')
      if (mnproc==1) then
         write(lp,*) 'Reading ICESTATE:',trim(icename)//'.[ab]'
      end if

      inquire(file=trim(icename)//'.b',exist=ex)
      if (.not.ex) then
         if (mnproc==1) write(lp,*) 'File ',trim(icename)//'.[ab] does not exist...'
         read_restart_icestate=.false.
         return
      end if

      call zaiopf(trim(icename)//'.a','old',icestate_unit)
      open(icestate_unit,file=trim(icename)//'.b',status='unknown')
      read(icestate_unit,6000) l_itdm
      read(icestate_unit,6000) l_jtdm
      read(icestate_unit,6000) l_kdm
      read(icestate_unit,6000) l_nthick
      read(icestate_unit,6000) l_nlaymax
      read(icestate_unit,6000) l_iyy
      read(icestate_unit,6000) l_idd
      read(icestate_unit,6000) l_ihh
      read (icestate_unit,'(a80)')  
 6000 format(a6,' =',i6)

      ! Skip forward  for ensemble runs
      if (imem>1) then 
         do ifld=1,(imem-1)*nflds
            call zaiosk(icestate_unit)
         end do
      end if

      ! Write ice variables
      ios = 0
      irec=1
      do while (ios==0 .and. irec <= nflds ) 
         read(icestate_unit,117) char8,hk,hl,xmin2,xmax2
         print *,char8,hk,hl

         if (trim(char8)=='nlay') then
            call zaiord(tmpfld,ip,.false., &
                        xmin,xmax, icestate_unit)
            !print *,maxval(tmpfld),minval(tmpfld)
            icestate%ice(hk)%nlay=1
            do j=1-margin,jj+margin
            do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
               !print *,i,j,hk,tmpfld(i,j)
               icestate(i,j)%ice(hk)%nlay=nint(tmpfld(i,j))
            end do
            end do
            end do
            !icestate%ice(hk)%nlay=nint(tmpfld)
         elseif (trim(char8)=='qstore') then
            call zaiord(icestate(:,:)%ice(hk)%qstore,ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='albs') then
            call zaiord(icestate%ice(hk)%albs,ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='fice') then
            call zaiord(icestate%ice(hk)%fice,ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='hice') then
            call zaiord(icestate%ice(hk)%hice,ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='hsnw') then
            call zaiord(icestate%ice(hk)%hsnw,ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='rhosnw') then
            call zaiord(icestate%ice(hk)%rhosnw,ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='tsrf') then
            call zaiord(icestate%ice(hk)%tsrf,ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='vtp') then
            call zaiord(icestate%ice(hk)%vtp(hl),ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='hml') then
            call zaiord(icestate%hml,ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='sml') then
            call zaiord(icestate%sml,ip,.false., &
                        xmin,xmax, icestate_unit)
         elseif (trim(char8)=='tml') then
            call zaiord(icestate%tml,ip,.false., &
                        xmin,xmax, icestate_unit)
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
      read_restart_icestate=.true.

 117  format (a8,' = ',i3,i3,1p2e16.7)

   end function read_restart_icestate




! ===================================================================================
! ============================= icestate_file =======================================
! ===================================================================================
   function icestate_file(rungen,rt,type)
     use mod_year_info
      implicit none

      character(len=40) icestate_file
      character(len=3) rungen
      character(len=*) type
      type(year_info) rt
      integer :: stlen

      stlen=len(trim(type))
      if (stlen>10) then
         print *,'To long string in icestate_file'
         stop '(icestate_file)'
      end if

      icestate_file=''
      icestate_file(1       :3       ) = rungen(1:3)
      icestate_file(4       :3 +stlen) = type(1:stlen)
      icestate_file(4+ stlen:7 +stlen) = rt%cyy(1:4)
      icestate_file(8+ stlen:8 +stlen) = '_'
      icestate_file(9+ stlen:11+stlen) = rt%cdd(1:3)
      icestate_file(12+stlen:12+stlen) = '_'
      icestate_file(13+stlen:14+stlen) = rt%chh(1:2)
      icestate_file(15+stlen:22+stlen) = 'ICESTATE'
   end function icestate_file

end module mod_icestate_io
