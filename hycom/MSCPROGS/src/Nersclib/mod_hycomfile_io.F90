module mod_hycomfile_io
! Module for reading several types of files produced by hycom,
! including restart, archv, nersc_daily and nersc_weekly. 
!
! Types: 
!   hycomfile - contains info on contents of the hycomfile
!
! Subroutines :
! --subroutine initHF(hycomfile, file, filetype, [template])
!   This normally is one of the first steps for setting up the
!   HYCOM I/O module. It initializes the hycomfile type by parsing
!   the contents of filename, given a hint on what file type this is.
!
! --subroutine HFReadField(hycomfile,field,idm,jdm,cfld,coord,tlevel)
!   Use this routine after calling initHF. It reads a 2D field "field"
!   with horizontal dimensions (idm,jdm), having the name "cfld",
!   at vertical coordinate  "coord" and time level "tlevel"
!
! --subroutine HFReadField3D(df,field,idm,jdm,kdm,cfld,tlevel)
!   Use this routine after calling initHF. It reads a 3D field "field"
!   with 3D  dimensions (idm,jdm,kdm), having the name "cfld",
!   at vertical coordinate  "coord" and time level "tlevel". It
!   uses the routine HFReadField
!
! --subroutine HFReadDPField(df,dp,idm,jdm,coord,tlevel)
!   Use this routine after calling initHF. Since variable names are
!   different across hycom files custom routines are available for some 
!   fields. This routine will read the layer thickness variable 
!   ("dp" in restart files, "thknss" in archive files, etc) of vertical
!   coordinate "coord" at time level "tlevel".
!
! --subroutine HFReadDPField_m(df,dp,idm,jdm,coord,tlevel)
!   Same as HFReadDPField, but always returns layer thickness in meters
!
! --subroutine HFReadDPField_p(df,dp,idm,jdm,coord,tlevel)
!   Same as HFReadDPField, but always returns layer thickness in pressure
!   coordinates
!
! --subroutine HFReaduvtot(df,ut,vt,idm,jdm,vlevel,tlevel)
!   Use this routine after calling initHF. Since variabe names are
!   different across hycom files custom routines are available for some 
!   fields. This routine will read the total velocity of layer 
!   coordinaye "coord" at time level "tlevel".
!   
! --subroutine HFReaduvbaro(df,ub,vb,idm,jdm,tlevel)
!   Use this routine after calling initHF. Since variabe names are
!   different across hycom files custom routines are available for some 
!   fields. This routine will read the barotropic velocity  at
!   time level "tlevel".
!
! --function getfiletype(filename)
!   Function "guesses" file type based on its file name - usually called
!   before initHF to set the filetype.
!
!
! Usage example:
! -------------------------------------------------------------------
!  use mod_hycomfile_io
!  ...
!  type(hycomfile) :: hfile
!  real, dimension(:,:), allocatable :: field
!  real, dimension(:,:,:), allocatable :: field3D
!
!  ! Retrieve grid size first - this uses the other modules in libhycnersc.a
!  call xcspmd()  ! sets idm, jdm
!  call zaiost()
!  call get_grid()
!
!  ! Reading a 2D field
!  fname='TP3restart1990_250_00.a' ! A test file
!  ftype=getfiletype(filename)     ! 'restart' in this case
!  call initHF(hfile, file, filetype)
!  kdm=vDim(hfile)   ! Retrieves vertical dimension in file
!
!  ! Allocate variables to be read
!  allocate(field3D(idm,jdm,kdm))
!  allocate(field(idm,jdm))
!
!  ! Read mixed layer thickness from that file
!  call HFReadField(hfile,field,idm,jdm'dpmixl   ',0,1)
! 
!  ! Read salinity in layer 1 from the same file
!  call HFReadField(hfile,field,idm,jdm'saln     ',1,1)
! 
!  ! Read temperature in all layers from the same file
!  call HFReadField(hfile,field3D,idm,jdm,kdm,'temp     ',1)
! 
! -------------------------------------------------------------------
! TODO: More info
! TODO: Loosened restriction on character length when calling routines (nb check this)
! TODO: ordinal day routine




use mod_parameters
implicit none

type hycomfile
   integer :: iversn   = 0
   integer :: iexpt    = 0
   integer :: nstep    = 0
   integer :: yrflag   = 0
   real*8  :: dtime    = 0
   real*8  :: fyear    = 0
   integer :: count    = 0
   integer :: nrec     = 0
   character(len=80) :: filebase ='', ftype=''
   character(len=8), pointer :: cfld (:)
   integer         , pointer, dimension(:) :: coord, tlevel

   ! Time info
   integer           :: iyear = 0
   integer           :: imonth= 0
   integer           :: iweek = 0
   integer           :: iday  = 0            !day starting at 0 in year
   integer           :: ihour = 0            !hour in day       (0-23)
   integer           :: imin  = 0            !minute in hour    (0-59)
   integer           :: isec  = 0            !second in minute  (0-59)
   real              :: ftime                !floating point year
   character(len=6)  :: ctime  = '000000'    !time (HHMMSS)
   character(len=8)  :: cdate  = '00000000'  !date (YYYYmmSS)

   ! Start info
   integer           :: start_iyear
   integer           :: start_iday
   integer           :: start_ihour = 0            !start hour in day       (0-23)
   integer           :: start_imin  = 0            !start minute in hour    (0-59)
   integer           :: start_isec  = 0            !start second in minute  (0-59)
   character(len=6)  :: start_ctime = '000000'     !start time (HHMMSS)
   character(len=8)  :: start_cdate = '00000000'   !start date (YYYYmmSS)
end type


logical, parameter, private :: silent=.true.

! TODO: Document
! TODO: Make sure DP routine returns in pressure coords always

private :: readHeader, readVarHeaders, readFieldEntry, writeFieldEntry

contains


   subroutine initHF(df,filename,ftype,template)
   implicit none
   type(hycomfile) , intent(inout) :: df
   character(len=*), intent( in) :: filename,ftype
   type(hycomfile) , intent( in),optional :: template
   integer :: ind

   if (present(template)) then
      df=template
      ind=max(index(filename,'.a'),index(filename,'.b'))
      df%filebase=filename(1:ind-1)
      df%ftype=ftype
      return
   end if

   ! Get file base
   ind=max(index(filename,'.a'),index(filename,'.b'))
   df%filebase=filename(1:ind-1)
   df%ftype=ftype

   print *, 'Read Header',df%filebase,df%ftype
   ! Read header
   call ReadHeader(df)
   !print *, 'Read var',df%iyear

   ! Get var info
   call readVarHeaders(df)
   end subroutine



!!!!!!!!!!! Header processing - parameters !!!!!!!!!!!!!



   subroutine readHeader(df)
   use mod_year_info
   implicit none
   type(hycomfile), intent(inout) :: df
   type(year_info)                :: tt1
   character(len=80) :: c80
   character(len=80) :: ctitle(4)
   character(len=5)  :: rforce
   integer :: nop,k,i, ind
   integer :: itime,itime0,year1,mon1,idm1,imon0,idm0
   integer :: iversn, iexpt, yrflag,lidm,ljdm,lkdm, dmonth
   integer :: syear, sday, dyear, dday, itmp, ios(3), indx
   nop=777

   !TODO : check date setup

   df%start_iyear=0
   df%start_iday =0

   open(nop,file=trim(df%filebase)//'.b',status='old')
   if (trim(df%ftype)=='restart') then
      read(nop,'(a)') c80 ; ind=index(c80,'='); 
      read(c80(ind+1:),*) df%iexpt, df%iversn, df%yrflag
      read(nop,'(a)') c80 ; ind=index(c80,'=');
      read(c80(ind+1:),*) df%nstep, df%dtime
      if (df%yrflag==3) then 
         call juliantodate(floor(df%dtime),df%iyear,itmp,itmp,1901,1,1)
         df%iday=floor(df%dtime) - datetojulian(df%iyear,1,1,1901,1,1) ! ordinal day
         df%ihour=nint((df%dtime - floor(df%dtime))*24.)
      else
         ! Due to a timing bug file name is our safest bet is to parse file
         ! names
         read(df%filebase(11:14),'(i4.4)',iostat=ios(1)) df%iyear
         read(df%filebase(16:18),'(i3.3)',iostat=ios(2)) df%iday
         read(df%filebase(20:21),'(i3.3)',iostat=ios(3)) df%ihour
         if (any(ios/=0)) then
            ! Drop it if that failed
            df%iyear=0
            df%iday=0
            df%ihour=0
         end if
      end if
      df%fyear = df%iyear + min((df%iday + df%ihour)/365.,1.)


   else if (trim(df%ftype)=='nersc_daily') then
      read(nop,116) ctitle,df%iversn,df%iexpt,df%yrflag, &
         lidm,ljdm,lkdm,df%start_iyear,df%start_iday ,  &
         df%iyear,df%iday,df%count
         df%ihour=12 ! hardcoded for daily average
         ! Floating point years - only useful as indicator
         df%fyear = df%iyear + min((df%iday + df%ihour)/365.,1.)

   else if (trim(df%ftype)=='nersc_weekly') then
      read(nop,216) ctitle,df%iversn,df%iexpt,df%yrflag, &
         lidm,ljdm,lkdm,dmonth,dyear , df%count
         read(df%filebase(8 :11),'(i4.4)',iostat=ios(1)) df%iyear
         read(df%filebase(13:14),'(i2.2)',iostat=ios(2)) df%imonth
         read(df%filebase(16:16),'(i1.1)',iostat=ios(3)) df%iweek
         df%iday=0
         df%ihour=0
         if (any(ios/=0) .or. df%iyear==9999) then
            ! Drop it if that failed
            df%iyear=0
            df%iday=0
            df%ihour=0
         end if
         !if iweek=9 -> monthly average
         if (df%iweek==9) then
            df%iweek=1
         end if
         !if iweek=9 -> monthly average
         if (df%imonth==99) then
            df%imonth=1
            df%iweek=1
         end if
         ! Floating point years - only useful as indicator
         df%fyear = df%iyear + (df%imonth -1) / 12. + (df%iweek -1. ) /(4.*12.) + 1./96.
         df%iday  = (df%imonth -1) / 12. + (df%iweek -1. ) /(4.*12.)  + 1/96.
         df%iday  = floor(df%iday*365.)
   else if (trim(df%ftype)=='archv'&
       .or.trim(df%ftype)=='archm'&
       .or.trim(df%ftype)=='archs') then
      !!check if old-style archive or new type
      do i=1,10
         read(nop,'(a)') c80
      end do
      close(nop)

      !!reopen b-file
      open(nop,file=trim(df%filebase)//'.b',status='old')
      if (c80(1:5)=='field') then
         !!old-style
         print*,' '
         print*,'*******************************'
         print*,'old style arch[v,m,s] header file'
         print*,'*******************************'
         print*,' '
         read(nop,316) ctitle,df%iversn,df%iexpt,df%yrflag
         !!get dump time from filename
         !!TODO: what if under 1 hour?
         read(df%filebase(7:10),'(i4.4)') df%iyear
         read(df%filebase(12:14),'(i3.3)') df%iday
         read(df%filebase(16:17),'(i2.2)') df%ihour
         print*,'year is ', df%iyear
         print*,'day is ',df%iday
         print*,'hour is ',df%ihour
         df%imin  = 0
         df%isec  = 0
         write(df%ctime,'(i2.2,i2.2,i2.2)') df%ihour,df%imin,df%isec
         !!
         call year_day(real(df%iday),df%iyear,tt1,'ecmwf')
         write(df%cdate,'(i4.4,i2.2,i2.2)')           &
            df%iyear,tt1%imm,1+tt1%idm

         !!set start time to finish time
         df%start_iyear = df%iyear
         df%start_iday  = df%iday
         df%start_ihour = df%ihour
         df%start_imin  = 0
         df%start_isec  = 0
         df%start_cdate = df%cdate
         df%start_ctime = df%ctime

      else
         !!new-style (more information in *.b files, like archv_wav)
         print*,' '
         print*,'********************************************'
         print*,'new style archv header file (like archv_wav)'
         print*,'********************************************'
         print*,' '

         read(nop,416) ctitle,df%iversn,df%iexpt,df%yrflag,lidm,ljdm    &
            ,itime0,df%start_cdate,df%start_ctime                       &
            ,itime,df%cdate,df%ctime

         !!year,month,day
         read(df%cdate,'(i4.4,i2.2,i2.2)') df%iyear,df%imonth,idm1
         df%iday  = datetojulian(df%iyear,df%imonth,idm1,df%iyear,1,1)

         !!time
         read(df%ctime,'(i2.2,i2.2,i2.2)') df%ihour,df%imin,df%isec

         !start date
         read(df%start_cdate,'(i4.4,i2.2,i2.2)') df%start_iyear,imon0,idm0
         df%start_iday  = datetojulian(df%start_iyear,imon0,idm0,df%start_iyear,1,1)

         !!start time
         read(df%start_ctime,'(i2.2,i2.2,i2.2)') df%start_ihour,df%start_imin,df%start_isec
      end if

      ! Floating point years - only useful as indicator
      df%fyear = df%iyear + min((df%iday + df%ihour/24.0)/365.,1.)
   else if (trim(df%ftype)=='archv_wav') then
      read(nop,416) ctitle,df%iversn,df%iexpt,df%yrflag,lidm,ljdm    &
         ,itime0,df%start_cdate,df%start_ctime                       &
         ,itime,df%cdate,df%ctime

      !!year,month,day
      read(df%cdate,'(i4.4,i2.2,i2.2)') df%iyear,df%imonth,idm1
      df%iday  = datetojulian(df%iyear,df%imonth,idm1,df%iyear,1,1)

      !!time
      read(df%ctime,'(i2.2,i2.2,i2.2)') df%ihour,df%imin,df%isec

      !start date
      read(df%start_cdate,'(i4.4,i2.2,i2.2)') df%start_iyear,imon0,idm0
      df%start_iday  = datetojulian(df%start_iyear,imon0,idm0,df%start_iyear,1,1)

      !!start time
      read(df%start_ctime,'(i2.2,i2.2,i2.2)') df%start_ihour,df%start_imin,df%start_isec

      ! Floating point years - only useful as indicator
      df%fyear = df%iyear + min((df%iday + df%ihour/24.0)/365.,1.)

   else
      print *,'readHeader > Unknown file type : '//trim(df%ftype)
      stop
   end if
   close (nop)
! nersc_daily format
116  format (a80/a80/a80/a80/ &
       i5,4x,'''iversn'' = hycom version number x10'/  &
       i5,4x,'''iexpt '' = experiment number x10'/  &
       i5,4x,'''yrflag'' = days in year flag'/  &
       i5,4x,'''idm   '' = longitudinal array size'/  &
       i5,4x,'''jdm   '' = latitudinal  array size'/  &
       i5,4x,'''kdm   '' = Vertical     array size'/  &
       i5,4x,'''syear '' = Year of integration start '/  &
       i5,4x,'''sday  '' = Day of integration start'/  &
       i5,4x,'''dyear '' = Year of this dump      '/  &
       i5,4x,'''dday  '' = Day of this dump     '/  &
       i5,4x,'''count '' = Ensemble counter       ')
! nersc_weekly format
 216  format (a80/a80/a80/a80/ &
     i5,4x,'''iversn'' = hycom version number x10'/ &
     i5,4x,'''iexpt '' = experiment number x10'/ &
     i5,4x,'''yrflag'' = days in year flag'/ &
     i5,4x,'''idm   '' = longitudinal array size'/ &
     i5,4x,'''jdm   '' = latitudinal  array size'/ &
     i5,4x,'''kdm   '' = Vertical     array size'/ &
     i5,4x,'''month '' = Month of this dump     '/ &
     i5,4x,'''year  '' = Year of this dump      '/ &
     i5,4x,'''count '' = Averaging counter      '/ &
     'field       time step  model day', &
     '  k  dens        min              max')

! Archive format
 316  format (a80/a80/a80/a80/                        &
     i5,4x,'''iversn'' = hycom version number x10'/   &
     i5,4x,'''iexpt '' = experiment number x10'/      &
     i5,4x,'''yrflag'' = days in year flag'/          &
     i5,4x,'''idm   '' = longitudinal array size'/    &
     i5,4x,'''jdm   '' = latitudinal  array size'/    &
     'field       time step  model day',              &
     '  k  dens        min              max')

! Archive waves format
 416  format (a80/a80/a80/a80/                        &
     i5,7x,'''iversn'' = hycom version number x10'/   &
     i5,7x,'''iexpt '' = experiment number x10'/      &
     i5,7x,'''yrflag'' = days in year flag'/          &
     i5,7x,'''idm   '' = longitudinal array size'/    &
     i5,7x,'''jdm   '' = latitudinal  array size'/    &
     i7,5x,'''dtime0'' = model start day'/            &
     a,4x, '''date0 '' = start date'/                 &
     a,6x, '''time0 '' = start time (HHMMSS)'/        &
     i7,5x,'''dtime '' = model dump day'/             &
     a,4x, '''date  '' = dump  date'/                 &
     a,6x, '''time  '' = dump  time (HHMMSS)'/        &
     'field       time step  model day',              &
     '  k  dens        min              max')

  end subroutine readHeader


    subroutine skipHeader(ftype,nop)
    implicit none
    character(len=*), intent(in)  :: ftype
    integer         , intent(inout)  :: nop
    character(len=5) :: char5
    integer :: ios
    if (    trim(ftype)=='nersc_weekly' .or. trim(ftype) == 'nersc_daily' &
        .or.trim(ftype)=='archv'.or.trim(ftype)=='archv_wav'&
        .or.trim(ftype)=='archm'.or.trim(ftype)=='archs') then
       ios=0 ; char5=''
       do while (char5/='field' .and. ios==0)
          read(nop,'(a5)',iostat=ios) char5
       end do
    else if (trim(ftype)=='restart' ) then
       read(nop,*) ; read(nop,*) 
    else
       print *,'skipHeader> unknown file type: '//trim(ftype)
       stop
    end if
    end subroutine

!!!!!!!!!!! Header processing - fields !!!!!!!!!!!!!

    subroutine readVarHeaders(df)
    implicit none
    type(hycomfile), intent(inout) :: df
    character(len=5) :: char5
    character(len=8) :: char8
    integer          :: ios, nop, nrec, coord, nstep, indx, irec,tlevel
    real             :: xmin, xmax
    logical          :: ex

    inquire(exist=ex,file=trim(df%filebase)//'.b')
    if (.not. ex) then
       print *,'file does not exist :'
       print *,trim(df%filebase)//'.b'
       stop '()'
    end if
    nop=777
    open(nop,file=trim(df%filebase)//'.b',status='old')

    ! Open input file
    call skipHeader(df%ftype,nop)

    ! Read until we get the index we want
    nrec=0 ; ios=0
    do while(ios==0)
       call readFieldEntry(df%ftype,char8,coord,tlevel,xmin,xmax,nop,ios)
       nrec=nrec+1
       !print *,nrec,char8,coord,tlevel,ios
    end do
    nrec=nrec-1; df%nrec=nrec

    rewind(nop) 
    call skipHeader(df%ftype,nop)
    allocate(df%cfld   (nrec))
    allocate(df%coord  (nrec))
    allocate(df%tlevel (nrec))
    ios=0
    do irec=1,nrec
       call readFieldEntry(df%ftype,char8,coord,tlevel,xmin,xmax,nop,ios)
       !print *,irec,char8,coord,tlevel
       df%cfld (irec)=char8
       df%coord  (irec)=coord
       df%tlevel (irec)=tlevel
    end do
    close(nop)
    end subroutine


!!!!!!!!!!! Header processing - one variable item  !!!!!!!!!!!!

    ! Read one line of variable info and parse it
    subroutine readFieldEntry(ftype,cfld,coord,tlevel,xmin,xmax,nop,ios)
    implicit none
    character(len=*), intent(in) :: ftype
    character(len=8), intent(out) :: cfld
    integer         , intent(out) :: coord,tlevel,ios
    real            , intent(out) :: xmin, xmax
    integer         , intent(in)  :: nop
    integer :: nstep
    real    :: rday,dens

    ! TODO: make sure this works properly - make it more robust
    if (trim(ftype)=='restart') then
       read(nop,4100,iostat=ios) cfld,coord,tlevel,xmin,xmax
    else if (trim(ftype)=="nersc_daily" .or. trim(ftype)=="nersc_weekly") then
       read(nop,117,iostat=ios) cfld,nstep,rday,coord,dens,xmin,xmax
       tlevel=1
    else if (trim(ftype)=="archv".or.trim(ftype)=="archv_wav"&
             .or.trim(ftype)=="archm".or.trim(ftype)=="archs") then
       read(nop,118,iostat=ios) cfld,nstep,rday,coord,dens,xmin,xmax
       tlevel=1
    else
       print *,'readFieldEntry> unknown file type: '//trim(ftype)
       stop
    end if
    4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
    117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
    118  format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)
    end subroutine




    ! Read one line of variable info and parse it
    subroutine writeFieldEntry(ftype,cfld,coord,tlevel,xmin,xmax,nop,ios)
    implicit none
    character(len=*), intent(in)  :: ftype
    character(len=*), intent(in)  :: cfld
    integer         , intent(in)  :: coord,tlevel
    real            , intent(in)  :: xmin, xmax
    integer         , intent(in)  :: nop
    integer         , intent(out) :: ios
    character(len=8) :: cfld2

    cfld2=cfld
    if (trim(ftype)=='restart') then
       write(nop,4100,iostat=ios) cfld2,coord,tlevel,xmin,xmax
    else if (trim(ftype)=="nersc_daily" .or. trim(ftype)=="nersc_weekly") then
       write(nop,117,iostat=ios) cfld2,0,0.,coord,0.,xmin,xmax
    else if (trim(ftype)=="archv".or.trim(ftype)=="archv_wav"&
             .or.trim(ftype)=="archm".or.trim(ftype)=="archs") then
       write(nop,118,iostat=ios) cfld2,0,0.,coord,0.,xmin,xmax
    else
       print *,'writeFieldEntry> unknown file type: '//trim(ftype)
       stop
    end if
    4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
    117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
    118  format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)
    end subroutine



!!!!!!!!!!! Header processing - write !!!!!!!!!!!!!

   subroutine HFwriteHeader(df,kdm)
   use mod_xc, only : idm, jdm
   implicit none
   type(hycomfile), intent(in) :: df
   integer        , intent(in) :: kdm
   character(len=80) :: c80, cline
   character(len=80) :: ctitle(4)
   integer :: nop,k,i, ind
   integer :: iversn, iexpt, yrflag,lidm,ljdm,lkdm, dmonth
   integer :: syear, sday, dyear, dday
   nop=777
   ctitle(1)='Generated by hycave/ensave'
   ctitle(2)='Generated by hycave/ensave'
   ctitle(3)='Generated by hycave/ensave'
   ctitle(4)='Generated by hycave/ensave'

   open(nop,file=trim(df%filebase)//'.b',status='replace')
   if (trim(df%ftype)=='restart') then
      write(nop,'(a,3i6)') 'RESTART: iexpt,iversn,yrflag = ', df%iexpt,df%iversn,df%yrflag
      write(cline,*)                df%nstep,df%dtime
      write(nop,'(a,a)')   'RESTART: nstep,dtime = ',cline(1:len_trim(cline))
      close(nop)
   else if (trim(df%ftype)=='nersc_daily') then
      write(nop,116) ctitle,df%iversn,df%iexpt,df%yrflag, &
         idm,jdm,kdm,df%start_iyear,df%start_iday,  &
         df%iyear,df%iday,df%count
         do k=1,1000
            write(nop,118) k,0
         end do
      write(nop,119)
   else if (trim(df%ftype)=='nersc_weekly') then
      write(nop,216) ctitle,df%iversn,df%iexpt,df%yrflag, &
         idm,jdm,kdm,df%iyear,df%imonth, df%count
   else if (trim(df%ftype)=='archv'&
       .or.trim(df%ftype)=='archm'&
       .or.trim(df%ftype)=='archs') then
      !TODO add archv_wav?
      write(nop,316) ctitle,df%iversn,df%iexpt,df%yrflag, &
         idm,jdm
   else 
      print *,'HFwriteHeader> unknown file type : '//trim(df%ftype)
      stop
   end if
   close (nop)
116  format (a80/a80/a80/a80/ &
       i5,4x,'''iversn'' = hycom version number x10'/  &
       i5,4x,'''iexpt '' = experiment number x10'/  &
       i5,4x,'''yrflag'' = days in year flag'/  &
       i5,4x,'''idm   '' = longitudinal array size'/  &
       i5,4x,'''jdm   '' = latitudinal  array size'/  &
       i5,4x,'''kdm   '' = Vertical     array size'/  &
       i5,4x,'''syear '' = Year of integration start '/  &
       i5,4x,'''sday  '' = Day of integration start'/  &
       i5,4x,'''dyear '' = Year of this dump      '/  &
       i5,4x,'''dday  '' = Day of this dump     '/  &
       i5,4x,'''count '' = Ensemble counter       ')
 216  format (a80/a80/a80/a80/ &
     i5,4x,'''iversn'' = hycom version number x10'/ &
     i5,4x,'''iexpt '' = experiment number x10'/ &
     i5,4x,'''yrflag'' = days in year flag'/ &
     i5,4x,'''idm   '' = longitudinal array size'/ &
     i5,4x,'''jdm   '' = latitudinal  array size'/ &
     i5,4x,'''kdm   '' = Vertical     array size'/ &
     i5,4x,'''month '' = Month of this dump     '/ &
     i5,4x,'''year  '' = Year of this dump      '/ &
     i5,4x,'''count '' = Averaging counter      '/ &
     'field       time step  model day', &
     '  k  dens        min              max')
118  format ('member ',i5.5,' = ',i1,' Ensemble member flag')
119  format('field         time step  model day', &
            '  k  dens        min              max')
! Archive format
 316  format (a80/a80/a80/a80/ &
     i5,4x,'''iversn'' = hycom version number x10'/ &
     i5,4x,'''iexpt '' = experiment number x10'/ &
     i5,4x,'''yrflag'' = days in year flag'/ &
     i5,4x,'''idm   '' = longitudinal array size'/ &
     i5,4x,'''jdm   '' = latitudinal  array size'/ &
     'field       time step  model day', &
     '  k  dens        min              max')
  end subroutine HFwriteHeader


   subroutine HFReadField3D(df,field,idm,jdm,kdm,cfld,tlevel)
   implicit none
   type(hycomfile),  intent(in)  :: df
   integer,          intent(in)  :: idm,jdm,kdm,tlevel
   real,             intent(out) :: field(idm,jdm,kdm)
   character(len=*), intent(in)  :: cfld
   integer :: k
   do k=1,kdm
      call HFReadField(df,field(1,1,k),idm,jdm,cfld,k,tlevel)
   end do
   end subroutine
      

   recursive subroutine HFReadField(df,field,idm,jdm,cfld,coord,tlevel)
   implicit none
   type(hycomfile),  intent(in)  :: df
   integer,          intent(in)  :: idm,jdm,coord,tlevel
   real,             intent(out) :: field(idm,jdm)
   character(len=*), intent(in)  :: cfld

   character(len=8) :: cfld2
   integer, dimension(idm,jdm) :: mask
   real, dimension(idm,jdm) :: dp
   real*4 :: A(idm,jdm), AMN, AMX, spval
   real :: xmin,dens,xmax,dens0,rday0
   integer :: indx,nstep0

   ! Some (unfortunate) corrections for ice variables
   cfld2=cfld
   if     (trim(df%ftype)=="nersc_daily" .and. trim(cfld2)=='hicem') then
      cfld2='hice'
   elseif (trim(df%ftype)=="nersc_daily" .and. trim(cfld2)=='ficem') then
      cfld2='fice'
   end if

   call indexFromHeader(df,cfld2,coord,tlevel,indx)
   if (indx/=-1) then
      spval=undef
      call READRAW(A,AMN,AMX,IDM,JDM,.false.,spval,trim(df%filebase)//'.a',indx)
      field=A
   else
      print '(a,i3.3)', 'Could not get field "'//cfld//'" at coordinate',coord
      field=undef
   end if
   where(field > 0.5*huge) field=0.

   if (trim(df%ftype)=="nersc_daily") then
      field=field/df%count
   else if (trim(df%ftype)=="nersc_weekly") then
      if (is3DVar(df,cfld2,tlevel) .and. .not. isDPVar(df,cfld2)) then
         call HFReadField(df,dp,idm,jdm,'pres    ',coord,tlevel)
         field=field/(max(dp,1e-4)*df%count) ! HFReadField averages, field is not averaged

         ! TODO - blank empty layers (<1e-4) ?
      else
         field=field/df%count
      end if
   end if
   end subroutine HFReadField


   subroutine HFReadDPField(df,dp,idm,jdm,coord,tlevel,units)
   implicit none
   type(hycomfile),  intent(in)  :: df
   integer,          intent(in)  :: idm,jdm,coord,tlevel
   real,             intent(out) :: dp(idm,jdm)
   real, dimension(idm,jdm) :: lw, up
   character(len=*), intent(in), optional :: units
   character(len=10) :: units2

   ! By default we use whatever is in the raw fields
   units2='native'
   if (present(units)) then
      if (trim(units) /= 'meter' .and.  trim(units) /= 'pressure' .and.  &
          trim(units) /= 'native' ) then
          print *,'Invalid unit sent to HFReadDPField'
          stop '(mod_hycomfile_io:HFReadDPField)'
       else
          units2=trim(units)
       end if
    end if


   if (trim(df%ftype)=='restart') then
      call HFReadField(df,dp,idm,jdm,'dp      ',coord,tlevel)
      if (trim(units2) == 'meter') dp=dp/onem
   elseif (trim(df%ftype)=="nersc_daily") then
      if (coord==1) then 
         up=0.
      else
         call HFReadField(df,up,idm,jdm,'pres    ',coord-1,tlevel)
      end if
      call HFReadField(df,lw,idm,jdm,'pres    ',coord,tlevel)
      dp=lw-up
      if (trim(units2) == 'meter') dp=dp/onem
   elseif (trim(df%ftype)=="nersc_weekly") then
      call HFReadField(df,dp,idm,jdm,'pres    ',coord,tlevel)
      if (trim(units2) == 'pressure') dp=dp*onem
   elseif (trim(df%ftype)=="archv".or.trim(df%ftype)=="archv_wav" .or. &
           trim(df%ftype)=="archm".or.trim(df%ftype)=="archs") then
      call HFReadField(df,dp,idm,jdm,'thknss  ',coord,tlevel)
      if (trim(units2) == 'meter') dp=dp/onem
   else
      print *,'HFReadDPField> unknown file type '//trim(df%ftype)
      stop
   end if
   end subroutine HFReadDPField

   subroutine HFReadDPField_p(df,dp,idm,jdm,coord,tlevel)
   implicit none
   type(hycomfile),  intent(in)  :: df
   integer,          intent(in)  :: idm,jdm,coord,tlevel
   real,             intent(out) :: dp(idm,jdm)
   call HFReadDPField(df,dp,idm,jdm,coord,tlevel,units='pressure')
   end subroutine HFReadDPField_p

   subroutine HFReadDPField_m(df,dp,idm,jdm,coord,tlevel)
   implicit none
   type(hycomfile),  intent(in)  :: df
   integer,          intent(in)  :: idm,jdm,coord,tlevel
   real,             intent(out) :: dp(idm,jdm)
   real, dimension(idm,jdm) :: lw, up
   call HFReadDPField(df,dp,idm,jdm,coord,tlevel,units='meter')
   end subroutine HFReadDPField_m



   ! Special routine for reading total velocity
   subroutine HFReaduvtot(df,ut,vt,idm,jdm,vlevel,tlevel)
   implicit none
   type(hycomfile),  intent(in)  :: df
   integer,          intent(in)  :: idm,jdm,tlevel, vlevel
   real,             intent(out) :: ut(idm,jdm)
   real,             intent(out) :: vt(idm,jdm)
   real, dimension(idm,jdm) :: ub, vb
   ub(:,:)=0.
   vb(:,:)=0.

   if (trim(df%ftype)=='restart') then
      call HFReadField(df,ub,idm,jdm,'ubavg   ',0,tlevel)
      call HFReadField(df,vb,idm,jdm,'vbavg   ',0,tlevel)
      call HFReadField(df,ut,idm,jdm,'ut      ',vlevel,tlevel)
      call HFReadField(df,vt,idm,jdm,'vt      ',vlevel,tlevel)
!      ut=ut+ub
!      vt=vt+vb
   elseif (trim(df%ftype)=="nersc_daily") then
      call HFReadField(df,ut,idm,jdm,'utot    ',vlevel,1)
      call HFReadField(df,vt,idm,jdm,'vtot    ',vlevel,1)
   elseif (trim(df%ftype)=="nersc_weekly") then
      call HFReadField(df,ut,idm,jdm,'utot    ',vlevel,1)
      call HFReadField(df,vt,idm,jdm,'vtot    ',vlevel,1)
   elseif(trim(df%ftype)=="archv"&
           .or.trim(df%ftype)=="archs") then
      call HFReadField(df,ut,idm,jdm,'u-vel.  ',vlevel,1)
      call HFReadField(df,vt,idm,jdm,'v-vel.  ',vlevel,1)
      call HFReadField(df,ub,idm,jdm,'u_btrop ',0,1)
      call HFReadField(df,vb,idm,jdm,'v_btrop ',0,1)
   elseif(trim(df%ftype)=="archm") then
      !Alfati. archm already contains total velocity in u-vel.&v-vel.
      call HFReadField(df,ut,idm,jdm,'u-vel.  ',vlevel,1)
      call HFReadField(df,vt,idm,jdm,'v-vel.  ',vlevel,1)
   elseif (trim(df%ftype)=="archv_wav") then
      call HFReadField(df,ut,idm,jdm,'u-vel.  ',vlevel,1)
      call HFReadField(df,vt,idm,jdm,'v-vel.  ',vlevel,1)
      call HFReadField(df,ub,idm,jdm,'u_btrop ',0,1)
      call HFReadField(df,vb,idm,jdm,'v_btrop ',0,1)
   else
      print *,'HFReaduvtot> unknown file type '//trim(df%ftype)
      stop
   end if
   ut=ut+ub
   vt=vt+vb
   end subroutine HFReaduvtot


   ! Special routine for reading barotropic velocity
   subroutine HFReaduvbaro(df,ub,vb,idm,jdm,tlevel)
   implicit none
   type(hycomfile),  intent(in)  :: df
   integer,          intent(in)  :: idm,jdm,tlevel
   real,             intent(out) :: ub(idm,jdm)
   real,             intent(out) :: vb(idm,jdm)
   real, dimension(idm,jdm) :: ut, vt, dp, dpsumu,dpsumv

   integer :: i,j,k, kdm
   real :: dpu, dpv

   kdm=vDim(df)
   if (trim(df%ftype)=='restart') then
      call HFReadField(df,ub,idm,jdm,'ubavg   ',0,1)
      call HFReadField(df,vb,idm,jdm,'vbavg   ',0,1)
   elseif (trim(df%ftype)=="nersc_daily") then

      ub=0.
      vb=0.
      dpsumu=0.
      dpsumv=0.
      do k=1,kdm
         call HFReadField(df,ut,idm,jdm,'utot    ',k,tlevel)
         call HFReadField(df,vt,idm,jdm,'vtot    ',k,tlevel)
         call HFReadField(df,dp,idm,jdm,'pres    ',k,tlevel)

         ! TODO - need periodic fix
         do j=2,jdm
         do i=2,idm
            dpu=0.5*(dp(i,j)+dp(i-1,j))
            dpv=0.5*(dp(i,j)+dp(i,j-1))
            ub(i,j)=ub(i,j)+dpu*ut(i,j)
            vb(i,j)=vb(i,j)+dpv*vt(i,j)
            dpsumu(i,j)= dpsumu(i,j) + dpu
            dpsumv(i,j)= dpsumv(i,j) + dpv
         end do
         end do
      end do

      do j=1,jdm
      do i=1,idm
         ub(i,j)=ub(i,j)/max(0.1,dpsumu(i,j))
         vb(i,j)=vb(i,j)/max(0.1,dpsumv(i,j))
      end do
      end do
      !print *,minval(ub),maxval(ub)
      !print *,minval(vb),maxval(vb)
   elseif (trim(df%ftype)=="nersc_weekly") then
      call HFReadField(df,ub,idm,jdm,'ubavg   ',0,1)
      call HFReadField(df,vb,idm,jdm,'vbavg   ',0,1)
   elseif (trim(df%ftype)=="archv"&
           .or.trim(df%ftype)=="archm"&
           .or.trim(df%ftype)=="archs") then
      call HFReadField(df,ub,idm,jdm,'u_btrop ',0,1)
      call HFReadField(df,vb,idm,jdm,'v_btrop ',0,1)
   elseif (trim(df%ftype)=="archv_wav") then
      call HFReadField(df,ub,idm,jdm,'u_btrop ',0,1)
      call HFReadField(df,vb,idm,jdm,'v_btrop ',0,1)
   else
      print *,'HFReaduvbaro> unknown file type '//trim(df%ftype)
      stop
   end if
   end subroutine HFReaduvbaro




   subroutine HFWriteField(df,field,idm,jdm,cfld,coord,tlevel,indx)
   implicit none
   type(hycomfile) , intent(in) ::df
   integer,          intent(in) :: idm,jdm,coord,indx,tlevel
   real,             intent(in) :: field(idm,jdm)
   character(len=*), intent(in) :: cfld
   integer, dimension(idm,jdm) :: mask
   real*4 :: A(idm,jdm), AMN, AMX,spval
   real :: xmin,xmax
   integer :: nop,ios

   A=field
   spval=undef
   call WRITERAW(A,AMN,AMX,IDM,JDM,.false.,spval,trim(df%filebase)//'.a',indx)
   xmax=AMX ; xmin=AMN
   open(nop,file=trim(df%filebase)//'.b',action='write',form='formatted',status='old',position='append',iostat=ios)
   call writeFieldEntry(df%ftype,cfld,coord,tlevel,xmin,xmax,nop,ios)
   close(nop)
   end subroutine HFWriteField




    subroutine HFHeaderFromIndex(df,indx,ierr,varname,level,timelevel)
    implicit none
    type(hycomfile) , intent(in)  :: df
    integer         , intent(in)  :: indx
    integer         , intent(out) :: ierr
    character(len=*), intent(out), optional :: varname
    integer         , intent(out), optional :: level
    integer         , intent(out), optional :: timelevel
    ierr=0
    if (indx> df%nrec .or. indx <1)  then
       ierr=-1
    else
       if (present(varname)) varname        =df%cfld (indx)
       if (present(level  )) level          =df%coord(indx)
       if (present(timelevel  )) timelevel  =df%tlevel(indx)
    end if
    end subroutine
      

    ! Retrieve a-file index for variable  
    subroutine indexFromHeader(df,cfld,coord,tlevel,indx)
    implicit none
    type(hycomfile) , intent(in)  :: df
    character(len=*), intent(in)  :: cfld
    integer         , intent(in)  :: coord,tlevel
    integer         , intent(out) :: indx
    integer :: irec
    indx=-1
    do irec=1,df%nrec
       if (trim(cfld)==trim(df%cfld(irec)) .and. coord==df%coord(irec) .and. df%tlevel(irec)==tlevel ) then
          indx=irec
       end if
    end do
    end subroutine


   subroutine HFUpdateAverage(dfave,df1)
   implicit none
   type(hycomfile), intent(in)    :: df1
   type(hycomfile), intent(inout) :: dfave
   integer :: k

   ! Consistency check and "addition" of restart average files
   if (df1%iversn/=dfave%iversn) then
      print *, 'iversn differ '
      stop '(df_plus_df)'
   end if

   if (df1%iexpt/=dfave%iexpt) then
      print *, 'iexpt differ '
      stop '(df_plus_df)'
   end if

   if (df1%yrflag/=dfave%yrflag) then
      print *, 'yrflag differ '
      stop '(df_plus_df)'
   end if

   ! Ambigous
   dfave%nstep=9999
   dfave%dtime=9999

   !dfave%count=df1%count+dfave%count
   dfave%count=1
   end subroutine


   logical function is3DVar(df,cfld,timelevel) 
   implicit none
   type(hycomfile) , intent(in) :: df
   character(len=*), intent(in) :: cfld
   integer         , intent(in) :: timelevel
   character(len=8) char8
   char8=adjustl(cfld)
   !print *,count( df%cfld == char8 .and. df%tlevel==timelevel ) 
   if (cfld=='chla_nor') then
     is3DVar=.true.
   else if (cfld=='chla_eco') then
     is3DVar=.true.
   else if(cfld=='attc_nor') then
     is3DVar=.true.
   else if(cfld=='attc_eco') then
     is3DVar=.true.
   else if(cfld=='nit_nor') then
     is3DVar=.true.
   else if(cfld=='nit_eco') then
     is3DVar=.true.
   else if(cfld=='pho_nor') then
     is3DVar=.true.
   else if(cfld=='sil_eco') then
     is3DVar=.true.
   else if(cfld=='pho_eco') then
     is3DVar=.true.
   else if(cfld=='pbio_nor') then
     is3DVar=.true.
   else if(cfld=='pbio_eco') then
     is3DVar=.true.
   else if(cfld=='zbio_nor') then
     is3DVar=.true.
   else if(cfld=='zbio_eco') then
     is3DVar=.true.
   else if(cfld=='oxy_nor') then
     is3DVar=.true.
   else if(cfld=='oxy_eco') then
     is3DVar=.true.
   else if(cfld=='salt1000') then
     is3DVar=.true.
   else if(cfld=='detvflux') then
     is3DVar=.true.
! _FABM__caglar_
   else if(cfld=='chla') then
     is3DVar=.true.
   else if(cfld=='nitrate') then
     is3DVar=.true.
   else if(cfld=='silicate') then
     is3DVar=.true.
   else if(cfld=='phosphat') then
     is3DVar=.true.
   else if(cfld=='pbiomass') then
     is3DVar=.true.
   else if(cfld=='zbiomass') then
     is3DVar=.true.
   else if(cfld=='oxygen') then
     is3DVar=.true.
   else if(cfld=='primprod') then
     is3DVar=.true.
   else if(cfld=='attcoeff') then
     is3DVar=.true.
   else if(cfld=='dic') then
     is3DVar=.true.
   else if(cfld=='ph') then
     is3DVar=.true.
   else if(cfld=='spco2') then
     is3DVar=.true.
! _FABM__caglar_
   else if(cfld=='utotl' .and. trim(df%ftype)=='archm') then
     is3DVar=.true.
   else if(cfld=='utotl' .and. trim(df%ftype)=='archv') then
     is3DVar=.true.
   else if(cfld=='wtotl') then
     is3DVar=.true.
   else        
     is3DVar=count( df%cfld == char8 .and. df%tlevel==timelevel ) > 1
   end if
   end function

   logical function isDPVar(df,cfld) 
   implicit none
   type(hycomfile) , intent(in) :: df
   character(len=*), intent(in) :: cfld
   if (trim(df%ftype)=='restart') then
      isDPVar=trim(cfld)=='dp'
   else if (trim(df%ftype)=='nersc_daily' .or. trim(df%ftype)=='nersc_weekly') then
      isDPVar=trim(cfld)=='pres'
   elseif (trim(df%ftype)=='archv'&
           .or.trim(df%ftype)=='archm'&
           .or.trim(df%ftype)=='archs') then
      isDPVar=trim(cfld)=='tknss'
   elseif (trim(df%ftype)=='archv_wav') then
      isDPVar=trim(cfld)=='tknss'
   else
      print *,'Unknown file type '//trim(df%ftype)
      stop '(mod_hycomfile_io:isDPVar)'
   end if
   end function

   integer function vDim(df) 
   implicit none
   type(hycomfile) , intent(in) :: df
   if (trim(df%ftype)=='restart') then
      vDim=count( df%cfld == 'dp      ' .and. df%tlevel==1 ) 
   else if (trim(df%ftype)=='nersc_daily' .or. trim(df%ftype)=='nersc_weekly') then
      vDim=count( df%cfld == 'pres    ' .and. df%tlevel==1 ) 
   elseif (trim(df%ftype)=='archv'&
           .or.trim(df%ftype)=='archm'&
           .or.trim(df%ftype)=='archs') then
      vDim=count( df%cfld == 'thknss  ' .and. df%tlevel==1 ) 
   elseif (trim(df%ftype)=='archv_wav') then
      vDim=count( df%cfld == 'thknss  ' .and. df%tlevel==1 ) 
   else
      print *,'Unknown file type '//trim(df%ftype)
      stop '(mod_hycomfile_io:isDPVar)'
   end if
   end function


   real function fyear(df) 
   implicit none
   type(hycomfile) , intent(in) :: df
   fyear=df%fyear
   end function


   real function jday1950(df) 
   use mod_year_info
   implicit none
   type(hycomfile) , intent(in) :: df

   ! datetojulian doesnt handle negative jday. Assumes iday is set to 0 on year flipover
   if (df%iyear>=1950) then
      jday1950=datetojulian(df%iyear,1,1,1950,1,1)  + df%iday
   else
      jday1950=datetojulian(1950,1,1,df%iyear,1,1)  + df%iday
   end if
      jday1950=jday1950+df%ihour/24.
   end function



   ! Returns file type based on file name
   function getfiletype(filename)
   implicit none
   character(len=*), intent(in) :: filename
   character(len=20) :: getfiletype
   integer :: findhdr,findab,finddaily,findweek,findrst              &
  &   ,findarchv,findarchv_wav,findarchm,findarchs

   ! Check for type ...
   findhdr  =index(filename,'.hdr')
   findab=max(index(filename,'.a'),index(filename,'.b'))
   if (.not. findab>0 .and. .not. findhdr>0) then 

      print *,'No .ab or .hdr files'
      stop '(mod_hycomfile_io:getfiletype)'

   else

      ! We have a .ab-file. Now figure out what type.
      findrst  =index(filename,'restart')
      finddaily=index(filename,'DAILY')
      findweek =index(filename,'AVE')
      findhdr  =index(filename,'.hdr')
      findarchv_wav=index(filename,'archv_wav')
      findarchv=index(filename,'archv.')
      findarchm=index(filename,'archm.')
      findarchs=index(filename,'archs.')
      if (findrst==4) then
         getfiletype='restart'
      elseif (finddaily==4) then
         getfiletype='nersc_daily'
      elseif (findweek==4) then
         getfiletype='nersc_weekly'
      elseif (findarchv_wav>0) then
         getfiletype='archv_wav'
      elseif (findarchv>0) then
         getfiletype='archv'
      elseif (findarchm>0) then
         getfiletype='archm'
      elseif (findarchs>0) then
         getfiletype='archs'
      elseif (findhdr>0) then
         getfiletype='pak'
         print *,'pak files no longer supported in this version'
         stop '(mod_hycomfile_io:getfiletype)'
      else
         print *,'Can not deduce file type from  file name'
         write(*,*)filename,findarchm
         stop '(mod_hycomfile_io:getfiletype)'
      end if
   end if

   end function getfiletype



      subroutine netcdfInfo(vnamein,gridrotate,stdname,longName,units,vname,cellmethod,limits)
      implicit none
      character(len=*), intent(in)  :: vnamein    ! name as in extract.[daily|weekly]
      logical         , intent(in)  :: gridrotate
      character(len=*) ,intent(out) :: stdName    ! standard_name attrib in netcdf file
      character(len=*) ,intent(out) :: longName   ! long_name attrib in netcdf file
      character(len=*) ,intent(out) :: units      ! Units attribute in netcdf file
      character(len=*) ,intent(out) :: vname      ! Name of variable in netcdf file
      character(len=*) ,intent(out) :: cellmethod ! Name of cell method
      real             ,intent(out) :: limits(2)  ! Lower and upper limits

      ! Standard name and unit lookup table for a given var name
      units=''
      cellmethod='area: mean'
      stdname=''
      longname=''
      limits=(/0,0/)
      select case (trim(vnamein))
      case ('model_depth')
         stdname='sea_floor_depth_below_sea_level' ; units='meter'; vname='model_depth'
         limits=(/0.,10001./)
      case ('saln','salin')
         stdname='sea_water_salinity' ; units='1e-3' ; vname='so'
         longname='Salinity'
         limits=(/0,45/)
      case ('temp') 
         stdname='sea_water_potential_temperature' ; units='degrees_C' ; vname='thetao'
         longname='Sea Temperature'
         limits=(/-3,50/)
      case ('levsaln')
         stdname='sea_water_salinity' ; units='1e-3' ; vname='levitus_salinity'
         limits=(/0,45/)
      case ('levtemp') 
         stdname='sea_water_potential_temperature' ; units='degrees_C' ; vname='levitus_temperature'
         limits=(/-3,50/)
      case ('ssh','srfhgt') 
         stdname='sea_surface_height_above_geoid' ; units='m' ; vname='zos'
         longname='Sea surface height'
         limits=(/-5,5/)
      case ('bsf','strmf') 
         stdname='ocean_barotropic_streamfunction' ; units='m3 s-1' ; vname='stfbaro'
         limits=(/-1e10,1e10/)
      case ('utot','utotl') 
         if (.not.gridrotate) then
            stdname='eastward_sea_water_velocity'
            vname='uo'
         else
            stdname='sea_water_x_velocity' 
            vname='vxo'
         end if
         units='m s-1' 
         limits=(/-3,3/)
      case ('vtot','vtotl') 
         if (.not.gridrotate) then
            stdname='northward_sea_water_velocity' 
            vname='vo'
         else
            stdname='sea_water_y_velocity' 
            vname='vyo'
         end if
         units='m s-1' 
         limits=(/-3,3/)
      case ('wtotl') 
         stdname='upward_sea_water_velocity' 
         units='m day-1' ; vname='wo'
         limits=(/-40,40/)
      case ('u','u-vel.') 
         if (.not.gridrotate) then
            stdname='baroclinic_eastward_sea_water_velocity' 
         else
            stdname='baroclinic_x_sea_water_velocity' 
         end if
         units='m s-1' ; vname='ubaroclin'
         limits=(/-3,3/)
      case ('v','v-vel.') 
         if (.not.gridrotate) then
            stdname='baroclinic_northward_sea_water_velocity' 
         else
            stdname='baroclinic_y_sea_water_velocity' 
         end if
         units='m s-1' ; vname='vbaroclin'
         limits=(/-3,3/)
      case ('hice','hicem','hi','hi_d') 
         stdname='sea_ice_thickness' ; units='m' ; vname='sithick'
         cellmethod='area: mean where sea_ice'
         limits=(/0,20/)
      case ('hsnw','hsnwm','hsnow','hs','hs_d') 
         stdname='surface_snow_thickness' ; units='m' ; vname='sisnthick'
         cellmethod='area: mean where sea_ice'
         limits=(/0,3/)
      case ('fice','ficem','aice','aice_d') 
         stdname='sea_ice_area_fraction' ; units='1' ; vname='siconc'
         limits=(/0,1/)
      case ('fy_age','iage_d')
          stdname='age_of_sea_ice' ; units='day' ; vname='siage'
          limits=(/0,36500/)
       case ('fy_frac','FYarea_d')
          stdname='sea_ice_classification' ; units='1' ; 
          longname = 'sea ice area fraction of first year ice'
          vname='siconc_fy'
          limits=(/0,1/)
      case ('ubavg','u_btrop') 
         if (.not.gridrotate) then
            stdname='barotropic_eastward_sea_water_velocity' 
         else
            stdname='barotropic_x_sea_water_velocity' 
         end if
         units='m s-1' ; vname='ubarotrop'
         limits=(/-3,3/)
      case ('vbavg','v_btrop') 
         if (.not.gridrotate) then
            stdname='barotropic_northward_sea_water_velocity' 
         else
            stdname='barotropic_y_sea_water_velocity' 
         end if
         units='m s-1' ; vname='vbarotrop'
         limits=(/-3,3/)
      case ('uice','uvel','uvel_d') 
         if (.not.gridrotate) then
            stdname='eastward_sea_ice_velocity' 
         else
            stdname='sea_ice_x_velocity' 
         end if
         cellmethod='area: mean where sea_ice'
         units='m s-1' ; vname='vxsi'
         limits=(/-3,3/)
      case ('vice','vvel','vvel_d') 
         if (.not.gridrotate) then
            stdname='northward_sea_ice_velocity' 
         else
            stdname='sea_ice_y_velocity' 
         end if
         cellmethod='area: mean where sea_ice'
         units='m s-1' ; vname='vysi'
         limits=(/-3,3/)
      case ('taux') 
         if (.not.gridrotate) then
            stdname='surface_downward_eastward_stress'
         else
            stdname='surface_downward_x_stress'
         end if
         units='pascal' ; vname='taux'
         limits=(/-3,3/)
      case ('tauy') 
         if (.not.gridrotate) then
            stdname='surface_downward_northward_stress'
         else
            stdname='surface_downward_y_stress'
         end if
         units='pascal' ; vname='tauy'
         limits=(/-3,3/)
     case ('tauxi') 
         if (.not.gridrotate) then
            stdname='downward_eastward_stress_at_sea_ice_base'
         else
            stdname='downward_x_stress_at_sea_ice_base'
         end if
         units='pascal' ; vname='tauxi'   
      case ('tauyi') 
         if (.not.gridrotate) then
            stdname='downward_northward_stress_at_sea_ice_base'
         else
            stdname='downward_y_stress_at_sea_ice_base'
         end if
         units='pascal' ; vname='tauyi'
      case ('qtot','surflx')
         units='W m-2' ; vname='qtot'
         limits=(/-3000.,3000./)
         stdname='surface_downward_heat_flux_in_sea_water'
      case ('swflx')
         units='W m-2' ; vname='swflx'
         stdname='net_downward_shortwave_flux_in_air'
      case ('emnp','fwflux')
         units='kg m-2 s-1' ; vname='fwflux'
         limits=(/-1e-3,1e-3/)
         stdname='water_flux_into_ocean'
      case ('stepmldT')
         units='m' ; vname='stepmldT'
         limits=(/0.,5000./)
         stdname='ocean_mixed_layer_thickness_defined_by_sigma_theta'
      case ('stepmld')
         units='m' ; vname='stepmld'
         limits=(/0.,5000./)
         stdname='ocean_mixed_layer_thickness_defined_by_sigma_theta'
!Alfati. More accurate method for computing MLD using density         
      case ('mld')
         units='m' ; vname='mlotst'
         limits=(/0.,3500./)
         stdname='ocean_mixed_layer_thickness_defined_by_sigma_theta'
      case ('dpmixl','dp_mixl','dpmix') 
         vname='dpmix'
         units='m'
         limits=(/0.,5000./)
         stdname='ocean_mixed_layer_thickness_defined_by_mixing_scheme'
!KAL20151204 - Adding bottom temperature
      case ('btemp') 
         vname    = 'bottomT'
         units    = 'degrees_C'
         limits   = (/-3,50/)
         stdname  = 'sea_water_potential_temperature_at_sea_floor'
         longname = 'Sea floor potential temperature'
!KAL20181210 - Adding bottom current coponents 
      case ('butot') 
         if (.not.gridrotate) then
            stdname='eastward_current_velocity' 
         else
            stdname='sea_current_u_velocity' 
         end if
         vname    = 'butot'
         units    = 'm s-1'
         limits   = (/-3,3/)
         stdname  = 'sea_water_u_component_at_sea_floor'
         longname = 'Sea floor u component'

      case ('bvtot') 
         if (.not.gridrotate) then
            stdname='northward_current_velocity' 
         else
            stdname='sea_current_v_velocity' 
         end if
         vname    = 'bvtot'
         units    = 'm s-1'
         limits   = (/-3,3/)
         stdname  = 'sea_water_v_component_at_sea_floor'
         longname = 'Sea floor v component'

      case ('bkinet') 
         vname    = 'bkinett'
         units    = '(m s-1)^2'
         limits   = (/-3,3/)
         stdname  = 'sea_water_kinetic_energy_at_sea_floor'
         longname = 'Sea floor kinetic energy'

!AS06092011 - adding biological variables for MyOcean
      case ('chla_nor') 
         vname='chla'
         units='kg m-3'
         limits=(/0.,1.e-4/)
         stdname='mass_concentration_of_chlorophyll_a_in_sea_water'
      case ('chla_eco')
         vname='chla'
         units='kg m-3'
         limits=(/0.,1.e-4/)
         stdname='mass_concentration_of_chlorophyll_a_in_sea_water'
      case ('attc_nor')
         vname='attcoef'
         units='m-1'
         limits=(/0.,0.5/)
         stdname='volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water'
      case ('attc_eco')
         vname='attcoef'
         units='m-1'
         limits=(/0.,0.5/)
         stdname='volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water'
      case ('nit_nor')
         vname='nitrat'
         units='mole m-3'
         limits=(/0.,5.e-2/)
      case ('nit_eco')
         vname='nitrat'
         units='mole m-3'
         limits=(/0.,5.e-2/)
         stdname='mole_concentration_of_nitrate_in_sea_water'
      case ('pho_nor')
         vname='phosphat'
         units='mole m-3'
         limits=(/0.,1e-2/)
         stdname='mole_concentration_of_phosphate_in_sea_water'
      case ('pho_eco')
         vname='phosphat'
         units='mole m-3'
         limits=(/0.,1e-2/)
         stdname='mole_concentration_of_phosphate_in_sea_water'
      case ('sil_eco')
         vname='silicate'
         units='mole m-3'
         limits=(/0.,100./)
         stdname='mole_concentration_of_silicate_in_sea_water'
      case ('pbio_nor')
         vname='pbiomass'
         units='mole m-3'
         limits=(/0.,0.01/)
         stdname='mole_concentration_of_phytoplankton_expressed_as_nitrogen_in_sea_water'
      case ('pbio_eco')
         vname='pbiomass'
         units='mole m-3'
         limits=(/0.,0.01/)
         stdname='mole_concentration_of_phytoplankton_expressed_as_nitrogen_in_sea_water'
      case ('zbio_nor')
         vname='zbiomass'
         units='mole m-3'
         limits=(/0.,0.01/)
         stdname='mole_concentration_of_zooplankton_expressed_as_nitrogen_in_sea_water'
      case ('zbio_eco')
         vname='zbiomass'
         units='mole m-3'
         limits=(/0.,0.01/)
         stdname='mole_concentration_of_zooplankton_expressed_as_nitrogen_in_sea_water'
      case ('oxy_nor')
         vname='oxygen'
         units='kg m-3'
         limits=(/0.,0.05/)
         stdname='mass_concentration_of_oxygen_in_sea_water'
      case ('oxy_eco')
         vname='oxygen'
         units='kg m-3'
         limits=(/0.,0.05/)
         stdname='mass_concentration_of_oxygen_in_sea_water'
      case ('salt1000')
         vname='salt1000'
         units='psu / 1000'
         limits=(/0.,0.045/)
         stdname='standard_salinity_divideby_1000'
      case ('detvflux')
         vname='expc'
         units='mol m-2 d-1'
         limits=(/0.0,1500.0/)
         stdname='sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water'
      case ('pp_d_nor')
         vname='pp_depth'
         units='kg m-2 s-1'
         limits=(/0.,1.e-7/)
         stdname='gross_primary_productivity_of_carbon'
      case ('pp_d_eco')
         vname='pp_depth'
         units='kg m-2 s-1'
         limits=(/0.,1.e-7/)
         stdname='gross_primary_productivity_of_carbon'
      case ('npp_euph')
         vname='npp_euph'
         units='g m-2 day-1'
         limits=(/0.,20./)
         stdname='net_primary_productivity_of_carbon_euphotic_depth'
      case ('chl_opti')
         vname='chl_opti'
         units='mg m-2'
         limits=(/0.,200./)
         stdname='depth_integrated_chlorophyll_one_optical_depth'
      case ('chlo_eco')
         vname='chl_opti'
         units='mg m-2'
         limits=(/0.,200./)
         stdname='depth_integrated_chlorophyll_one_optical_depth'
!AS06092011
      case ('albedo','albsni_d')
         vname='sialb'
         units='1'
         limits=(/0.,1./)
         stdname='sea_ice_albedo'

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !TW: waves-in-ice quantities
      case ('dfloe')
         vname='dmax'
         units='m'
         limits=(/0.,300./)
         stdname='sea_ice_max_floe_size'
         cellmethod='area: max where sea_ice'
      case ('swh')
         vname='swh'
         units='m'
         limits=(/0.,20./)
         stdname='combined_wind_swell_significant_wave_height'
         cellmethod='area: mean'
      case ('mwp')
         vname='mwp'
         units='s'
         limits=(/0.,40./)
         stdname='primary_mean_wave_period'
         cellmethod='area: mean'
      case ('mwd')
         vname='mwd'
         units='degrees'
         limits=(/0.,360./)
         stdname='mean_wave_from_direction'
         cellmethod='area: mean'
! _FABM__caglar_
         case ('chla')
         vname='chl'
         units='mg m-3'
         limits=(/0.,100./)
         stdname='mass_concentration_of_chlorophyll_a_in_sea_water'
         case ('nitrate')
         vname='no3'
         units='mmol m-3'
         limits=(/0.,50./)
         stdname='mole_concentration_of_nitrate_in_sea_water'
         case ('silicate')
         vname='si'
         units='mmol m-3'
         limits=(/0.,250./)
         stdname='mole_concentration_of_silicate_in_sea_water'
         case ('phosphat')
         vname='po4'
         units='mmol m-3'
         limits=(/0.,10./)
         stdname='mole_concentration_of_phosphate_in_sea_water'
         case ('pbiomass')
         vname='phyc'
         units='mmol m-3'
         limits=(/0.,500./)
         stdname='mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water'
         case ('zbiomass')
         vname='zooc'
         units='mmol m-3'
         limits=(/0.,500./)
         stdname='mole_concentration_of_zooplankton_expressed_as_carbon_in_sea_water' 
         case ('oxygen')
         vname='o2'
         units='mmol m-3'
         limits=(/0.,1000./)
         stdname='mole_concentration_of_dissolved_molecular_oxygen_in_sea_water'
         case ('pp_depth')
         vname='npp'
         units='mg m-2 d-1'
         limits=(/0.,8460./)
         stdname='net_primary_productivity_of_biomass_expressed_as_carbon'
         case ('primprod')
         vname='nppv'
         units='mg m-3 day-1'
         limits=(/0.,2000./)
         stdname='net_primary_production_of_biomass_expressed_as_carbon_per_unit_volume_in_sea_water'
         case ('attcoeff')
         vname='kd'
         units='m-1'
         limits=(/0.,1./)
         stdname='volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water'
         case ('dic')
         vname='dissic'
         units='mole m-3'
         limits=(/1.,3./)
         stdname='mole_concentration_of_dissolved_inorganic_carbon_in_sea_water'
         case ('ph')
         vname='ph'
         units='1'
         limits=(/7.,10./)
         stdname='sea_water_ph_reported_on_total_scale'
         case ('spco2')
         vname='spco2'
         units='Pa'
         limits=(/0.,100./)
         stdname='surface_partial_pressure_of_carbon_dioxide_in_sea_water'
! _FABM__caglar_
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !TW: surface currents - if don't want to dump full 3d velocities
      case ('usurf') 
         if (.not.gridrotate) then
            stdname='surface_eastward_sea_water_velocity' 
            longname='eastward_sea_water_velocity_at_3m_depth' 
            ! TODO
            ! *with uvsurf_opt=2 in m_archiv_select.F of HYCOM code,
            ! this is the bottom of 1st layer
            ! - read dp00 from blkdat.input to get corresponding depth?
            ! - extract top layer thickness?
            !
            ! *with uvsurf_opt=1 in m_archiv_select.F of HYCOM code,
            ! this is the average over the top 'surf_depth' meters
            ! - how to get these parameters?
         else
            stdname='surface_x_sea_water_velocity' 
            longname='x_sea_water_velocity_at_3m_depth' 
            !TODO see previous TODO
         end if
         units='m s-1'
         vname='usurf'
         limits=(/-3,3/)
      case ('vsurf') 
         if (.not.gridrotate) then
            stdname='surface_northward_sea_water_velocity' 
            longname='northward_sea_water_velocity_at_3m_depth' 
            !TODO see previous TODO
         else
            stdname='surface_y_sea_water_velocity' 
            longname='y_sea_water_velocity_at_3m_depth' 
            !TODO see previous TODO
         end if
         units='m s-1'
         vname='vsurf'
         limits=(/-3,3/)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case default
         vname=trim(vnamein)
      end select
      end subroutine netcdfInfo



      !TODO: Fix for different yrflag
      subroutine forecastDate(hfile,rt)
      use mod_year_info
      implicit none
      type(hycomfile), intent(in)  :: hfile
      type(year_info), intent(out) :: rt
      real rday
      integer diy
      integer  :: ihour,imin,isec

      if (trim(hfile%ftype)=='nersc_daily') then
         call year_day(real(hfile%iday),hfile%iyear,rt,'ecmwf')
      else if (trim(hfile%ftype)=='nersc_weekly') then
         ! Only approxmate, day in year
         diy=(hfile%iweek-1)*7 + (hfile%imonth-1)*30 + floor(4.*hfile%imonth/12.)
         call year_day(real(diy),hfile%iyear,rt,'ecmwf')
      else if (trim(hfile%ftype)=='restart') then
         !TW: restarts can be any time in the day
         rday  = real(hfile%iday) + real(hfile%ihour)/24.0
         call year_day(rday,hfile%iyear,rt,'ecmwf')
         !call year_day(real(hfile%iday),hfile%iyear,rt,'ecmwf')
      else if (trim(hfile%ftype)=='archv'&
          .or.trim(hfile%ftype)=='archm'&
          .or.trim(hfile%ftype)=='archs') then
         !rday  = real(hfile%iday) + real(hfile%ihour)/24.0
         !rday  = hfile%iday+(3600*hfile%ihour+60*hfile%imin+hfile%isec)/(24.*3600.)
!Alfati   substract 1 to fix the date fro TOPAZ5         
         rday  = hfile%iday+(3600*hfile%ihour+60*hfile%imin+hfile%isec)/(24.*3600.)-1
         call year_day(rday,hfile%iyear,rt,'ecmwf')
      else if (trim(hfile%ftype)=='archv_wav') then
         read(hfile%start_ctime,'(i2.2,i2.2,i2.2)') ihour,imin,isec
         rday  = hfile%iday+(3600*hfile%ihour+60*hfile%imin+hfile%isec)/(24.*3600.)
         call year_day(rday,hfile%iyear,rt,'ecmwf')
      else
         write(6,*) 'forecastDate> Unknown file type '//trim(hfile%ftype)
         call exit(1)
      end if
      end subroutine forecastDate
         

      !TODO: Fix for different yrflag
      subroutine startDate(hfile,rt)
      use mod_year_info
      implicit none
      type(hycomfile), intent(in)  :: hfile
      type(year_info), intent(out) :: rt
      integer  :: ihour,imin,isec
      real*8   :: dtime

      if (trim(hfile%ftype)=='nersc_daily') then
         call year_day(real(hfile%start_iday),hfile%start_iyear,rt,'ecmwf')
      else if (trim(hfile%ftype)=='nersc_weekly') then
         call forecastDate(hfile,rt)
      else if (trim(hfile%ftype)=='restart') then
         call forecastDate(hfile,rt)
! CAGLAR - please check this - not sure about including archm here
!          as below there are time conversions
!          - I am very focused on creating outputs that I will leave
!          this without inspection at the moment
      else if (trim(hfile%ftype)=='archv'&
               .or.trim(hfile%ftype)=='archm'&
               .or.trim(hfile%ftype)=='archs') then

         read(hfile%start_ctime,'(i2.2,i2.2,i2.2)') ihour,imin,isec
         !dtime = hfile%start_iday+(3600*ihour+60*imin+isec)/(24.*3600.)
!Alfati   substract 1 to fix the date fro TOPAZ5         
         dtime = hfile%start_iday+(3600*ihour+60*imin+isec)/(24.*3600.)-1
         call year_day(dtime,hfile%start_iyear,rt,'ecmwf')
         !call forecastDate(hfile,rt)

      else if (trim(hfile%ftype)=='archv_wav') then

         read(hfile%start_ctime,'(i2.2,i2.2,i2.2)') ihour,imin,isec
         dtime = hfile%start_iday+(3600*ihour+60*imin+isec)/(24.*3600.)
         call year_day(dtime,hfile%start_iyear,rt,'ecmwf')

      else
         write(6,*) 'startDate> Unknown file type '//trim(hfile%ftype)
         call exit(1)
      end if
      end subroutine startDate
         




! Modified from Alan Wallcraft's RAW routine by Knut Liseter @ NERSC
! So far only the "I" in "IO" is present
      SUBROUTINE READRAW(A,AMN,AMX,IDM,JDM,LSPVAL,SPVAL,CFILE1,K)
      IMPLICIT NONE
!
      REAL*4     SPVALH
      PARAMETER (SPVALH=1e30)
!
      REAL*4,           INTENT(OUT) :: A(IDM,JDM)
      REAL*4,           INTENT(OUT) :: AMN,AMX
      INTEGER,          INTENT(IN)  :: IDM,JDM
      LOGICAL,          INTENT(IN)  :: LSPVAL
      REAL*4,           INTENT(INOUT)  :: SPVAL
      INTEGER,          INTENT(IN)  :: K
      CHARACTER(len=*), INTENT(IN)  :: CFILE1
!
      REAL*4 :: PADA(4096)
!     MOST OF WORK IS DONE HERE.
!

      INTEGER      I,J,IOS,NRECL
      INTEGER NPAD
!
      IF(.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
!
!!! Calculate the number of elements padded!!!!!!!!!!!!!!!!!!!!!!!!
      NPAD=GET_NPAD(IDM,JDM)
      INQUIRE( IOLENGTH=NRECL) A,PADA(1:NPAD)
!     
!     
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='old', &
               ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
!
      READ(11,REC=K,IOSTAT=IOS) A
      close(11)
!
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record ',K, &
                   ' from '//CFILE1(1:LEN_TRIM(CFILE1))
        CALL EXIT(4)
      ENDIF
!
      AMN =  SPVALH
      AMX = -SPVALH
      DO J= 1,JDM
      DO I=1,IDM
         IF     (A(I,J).LE.SPVALH) THEN
            AMN = MIN( AMN, A(I,J) )
            AMX = MAX( AMX, A(I,J) )
         ELSEIF (LSPVAL) THEN
            A(I,J) = SPVAL
         ENDIF
      END DO
      END DO
!                 
      RETURN
      END SUBROUTINE



! Modified from Alan Wallcraft's RAW routine by Knut Liseter @ NERSC
      SUBROUTINE WRITERAW(A,AMN,AMX,IDM,JDM,LSPVAL,SPVAL,CFILE1,K)
      IMPLICIT NONE
      REAL*4     SPVALH
      PARAMETER (SPVALH=1e30)
      REAL*4,        INTENT(INOUT) :: A(IDM,JDM)
      REAL*4,        INTENT(OUT)   :: AMN,AMX
      INTEGER,       INTENT(IN)    :: IDM,JDM
      LOGICAL,       INTENT(IN)    :: LSPVAL
      REAL*4,        INTENT(INOUT) :: SPVAL
      INTEGER,       INTENT(IN)    :: K
      CHARACTER(len=*), INTENT(IN) :: CFILE1
!
      REAL*4 :: PADA(4096)
!
!     MOST OF WORK IS DONE HERE.
!

      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,IOS,NRECL
      INTEGER NPAD
!
      IF(.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
!
!!! Calculate the number of elements padded!!!!!!!!!!!!!!!!!!!!!!!!
      NPAD=GET_NPAD(IDM,JDM)
      PADA=0.
      INQUIRE( IOLENGTH=NRECL) A,PADA(1:NPAD)
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='unknown', &
               ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
!
      WRITE(11,REC=K,IOSTAT=IOS) A,PADA(1:NPAD)
      close(11)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t write record ',K, &
                   ' from '//CFILE1(1:LEN_TRIM(CFILE1))
        CALL EXIT(4)
      ENDIF
!
      AMN =  SPVALH
      AMX = -SPVALH
      DO J= 1,JDM
      DO I=1,IDM
         IF     (A(I,J).LE.SPVALH) THEN
            AMN = MIN( AMN, A(I,J) )
            AMX = MAX( AMX, A(I,J) )
         ELSEIF (LSPVAL) THEN
            A(I,J) = SPVAL
         ENDIF
      END DO
      END DO
!                 
      RETURN
      END SUBROUTINE




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







      INTEGER FUNCTION GET_NPAD(IDM,JDM)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IDM,JDM
         GET_NPAD = 4096 - MOD(IDM*JDM,4096)
         GET_NPAD = mod(GET_NPAD,4096)
      END FUNCTION


end module
