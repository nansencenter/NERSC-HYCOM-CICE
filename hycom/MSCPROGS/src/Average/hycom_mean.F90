program hycom_mean
  !
  ! F90 version of hycom_All/hycom_mean
  !
  use mod_mean  ! HYCOM mean array interface
  use mod_za    ! HYCOM array I/O interface
  use mod_blkin
  implicit none

  ! --- Form the mean (or mean-squared) of a sequence of HYCOM archive files.
  ! --- Layered means weighted by layer thickness.

  character(len=81) :: label
  character(len=18) :: text
  character(len=80) :: flnm
  !logical :: meansq,single,trcout,icegln,hisurf
  integer :: mntype,iweight,i,j,n
  integer :: narchs,iarch
  integer :: iexpt,jexpt,kkin,yrflag,kpalet,mxlflg
  real    :: thbase
  real(kind=8) :: time_min,time_max,time_ave,time(3)

  call xcspmd()  !define idm,jdm
  call zaiost()
  !lp=6

  iexpt  = 0
  ii     = idm
  jj     = jdm
  iorign = 1
  jorign = 1
  
  ! --- 'ntracr' = number of tracers to include in the mean
  ! ---              optional, default is 0
  ! --- 'kk    ' = number of layers involved (1 for surface only means)

  call blkini2(i,j,  'ntracr','kk    ')  !read ntracr or kk as integer
  if (j.eq.1) then !'ntracr'
     ntracr = i
     call blkini(kk,    'kk    ')
  else !kk
     ntracr = 0
     kk     = i
  endif

  hisurf = .true.
  trcout = ntracr.gt.0
  icegln = .true.

  ! --- array allocation and initialiation

  if (ntracr_bgc_3d.gt.0 .or. ntracr_bgc_2d.gt.0) then
     bgcout = .true.
  endif
  
  call mean_alloc

  ! --- land masks.

  call getdepth('regional.depth')
  call bigrid(depths)
  call mean_depths

  ! --- read and sum a sequence of archive files.

  time_min =  huge(time_min)
  time_max = -huge(time_max)
  time_ave =  0.0
  nstep_m = 0

  ! --- 'single' = treat input mean archive as a single sample
  ! ---              optional, default is meansq
  ! --- 'meansq' = form meansq (rather than mean)

  call blkini2(i,j,'single','meansq')  !read single or meansq as integer
  if (j.eq.1) then !'single'
     single = i .ne. 0  !0=F
     call blkinl(meansq,'meansq')
     if (meansq .and. .not.single) then
        write (lp,'(a)') 'error: meansq=T requires single=T'
        call flush(lp)
        stop
     endif !error
  else !'meannsq'
     meansq = i .ne. 0  !0=F
     single = meansq
  endif

  do ! loop until input narchs==0
     ! ---   'narchs' = number of archives to read (==0 to end input)
     call blkini(narchs,'narchs')
     if (narchs.eq.0) then
        exit
     endif
     do iarch= 1,narchs
        read (*,'(a)') flnm
        write (lp,'(2a)') ' input file: ',trim(flnm)
        call getdat(flnm,time,iweight,mntype,iexpt,yrflag,thbase)
        if (iweight.eq.1) then
           call mean_velocity  ! calculate full velocity and KE
        endif
        if (single .and. mntype.eq.1) then
           iweight = 1  !treat mean as a standard archive
        endif
        if (iweight.eq.1 .and. meansq) then
           call mean_addsq(iweight)
        else
           call mean_add(iweight)
        endif
        time_min = min( time_min,  time(1) )
        time_max = max( time_max,  time(2) )
        time_ave =      time_ave + time(3) * iweight
        nstep_m  =      nstep_m  + nstep   * iweight
     enddo
  enddo

  ! --- output the mean or mean square

  call mean_end
  time_ave = time_ave/nmean
  nstep_m = nstep_m/nmean

  read (*,'(a)') flnm
  write (lp,'(2a)') 'output file: ',trim(flnm)
  call flush(lp)
  if (meansq) then
     mntype = 2
  else
     mntype = 1
  endif
  jexpt=iexpt
  call putdat(flnm,time_min,time_max,time_ave,mntype,iexpt,jexpt,yrflag,thbase)

end program hycom_mean
   
