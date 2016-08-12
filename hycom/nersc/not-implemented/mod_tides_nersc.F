      module mod_tides_nersc
      use mod_xc 
      integer, parameter     :: nconst=8
      integer, parameter     :: nastro=8
      logical, save          :: ltides  ! Tidal flag
      logical, save          :: ltideuv ! Include tidal currents in addition to sea level
      character(len=3), save :: ctide    ! Tidal data (FES or CSR)

c --- nconst: Number of tidal constituents
c --- nastro: Number of tidal astronomical arguments,
c ---         Please consult subroutines tide_ast.f
c ---         and tiderSB.f before you even think about
c ---         to change nastro!!

      type tides
         real sbc(2*(itdm+jtdm),2*nconst)
         real sbc_u(2*(itdm+jtdm),4*nconst)
         real freq(nconst)
         real h(2*(itdm+jtdm))
         real u(2*(itdm+jtdm))
         real v(2*(itdm+jtdm))
         real astro_arg(nconst+1)
         integer iastro(nastro+1)
         logical astro_calc
         logical uv
      end type
      type(tides), public :: tide
c
      logical, public :: lslowstart=.true.   ! Only used for tidal spinup
      logical, public :: lastrotide=.true.   ! swithc for synoptic tides
      real, public :: astr_dt=5.0 
      integer, public :: astr_time
      real*8, save :: nersctide_time
c
      contains
c
c --- -------------------------------------------------------
c --- Author: K. Simonsen, June, 1999, neRSC
c --- ------------------------------------------------
c --- This routine reads frequencies, amplitude and phases for
c --- namph tidal constituents from file xxx/AMPPHA.USR
c --- for use in tidal boundary conditions. The input file
c --- is prepared by program CSR2mod.f.
c --- 
c --- -------------------------------------------------------
      subroutine tiderSB(amp_cnv,v_cnv) !,ctide,ltideuv)
      implicit none
      integer, parameter :: apn=2*nconst    !16
      integer, parameter :: apnu=4*nconst   !32
c
      real dec_conv
      real, intent(in) :: amp_cnv
      real, intent(in) :: v_cnv
c
      integer ii,jj
      integer iif,jjf,kkf,ncf,mtot,i,j,k,m,mm
      real ap(apn),ap_u(apnu),freqa,freq_cnv,phase_cnv
      character dummy*80,infname*21,hname*4,astro_name(nastro)*4
      logical astro_loc
      logical ex
      character(len=22) :: fname_elev ! elevation file
      character(len=22) :: fname_curr ! current file
c
      if (mnproc==1) then
      print '(a)','####################################################'
      print '(a)','Prepare for tides "tiderSB"'
      end if
c
c --- ----------------------------- Constituents names. ------
c ---    Please consult tide_ast.f
c ---    before the order below is changed.
      astro_name(1)='Q1  '
      astro_name(2)='O1  '
      astro_name(3)='P1  '
      astro_name(4)='K1  '
      astro_name(5)='N2  '
      astro_name(6)='M2  '
      astro_name(7)='S2  '
      astro_name(8)='K2  '
      if (lastrotide) then
         tide%astro_calc=.true.
      else
         tide%astro_calc=.false.
      endif
      astro_loc=.false.
      freq_cnv = datan(1.D0)*8.                  !2pi   rad
c
c --- ------------------------------ Initialize arrays ---------
      tide%sbc=-999.
      tide%sbc_u=-999.
      tide%h=0.0
      tide%u=0.0
      tide%v=0.0
      tide%freq=0.
      tide%iastro=nconst+1
c
c --- Load Sea Level data from file 
      !fname_elev(1:7)='./Data/'
      !fname_elev(8:10)=ctide(1:3)
      !fname_elev(11:22)='obc_elev.dat'
      fname_elev=ctide(1:3)//'obc_elev.dat'
      inquire(file=trim(fname_elev),exist=ex)
      if (ex) then
         if (mnproc==1)print '(a)','Loading file: '//trim(fname_elev)
      else
         if (mnproc==1)
     &      print '(a)','Tidal file does not exist: '//trim(fname_elev)
         call xcstop('(mod_tides_nersc)')
      endif
      open(10,file=trim(fname_elev),STATUS='OLD')
      read(10,'(a80)')dummy                      !Read text line
c
c --- Check model grid dimensions
      read(10,'(A36,2I4)')dummy(1:36),iif,jjf  
      if ((iif /= itdm).OR.(jjf /= jtdm)) then
         if (mnproc==1) then
         write(lp,'(''Model grid dim.    (itdm,jtdm):'',2I4)')itdm,jtdm
         write(lp,'(''OBC data grid dim. (jj,kk):'',2I4)')jjf,kkf
         end if
         call xcstop('(mod_tides_nersc)')
      endif
      ii=itdm
      jj=jtdm
c
c --- Check no off constituents in file
      read(10,'(A36,I2)')dummy(1:36),ncf 
      if (ncf /= nconst) then
         if (mnproc==1) then
         write(lp,'(a,2i5)')'Wrong number of tidal constituents in file'
     .            ,ncf,nconst
         end if
        call xcstop('(mod_tides_nersc)')
        stop
      endif
      read(10,'(a80)')dummy                      !Read text line
c
c --- Read header info about constituents
      if (mnproc==1) write(lp,'(''tiderSB: No Name Period (hours) '')')
      do i=1,ncf
         read(10,'(A24,A4,f12.9)')dummy(1:24),hname,freqa  !Read freq,name
         if (mnproc==1)write(lp,'(I2,tr2,a3,1x,f12.9)')i,hname,24./freqa
c
c ---    Convert from cyc/day to 2pi /day
         tide%freq(i)=freqa*freq_cnv      
c
c ---    Check if astr. arg. is available
         astro_loc=.false.
         do j=1,nastro                       
            if (astro_name(j) == hname) then
               tide%iastro(i)=j
               astro_loc=.true.
            endif
         enddo
c
         if (.not.astro_loc) then
            if (mnproc==1) 
     &      write(lp,'(a)')'tiderSB:Astro. arg. is not available for '//
     &                    hname
            tide%astro_calc=.false.
         endif
      enddo
c
c --- If not all astro argements are available set all indices to zero argument
      if (.not.tide%astro_calc) then     
         do i=1,nconst
            tide%iastro(i)=nconst+1    
         enddo
         if (mnproc==1) 
     &      write(lp,'(a)')'Astro arguments are put to zero!! '
      endif
      read(10,'(a80)')dummy                      !Read text line
      mtot=2*(ii+jj)
      phase_cnv = datan(1.D0)/45.                !2pi/360
c
      if (ctide(1:3) == 'FES') then
         dec_conv=1.0
      else
         dec_conv=1.0
      endif
c
      do i=1,mtot
         read(10,'(i4,16f8.2)') m,ap
         do j=1,nconst
            tide%sbc(i,2*j-1) = ap(2*j-1)*amp_cnv                !Amplitude
            tide%sbc(i,2*j)   = ap(2*j)  *phase_cnv*dec_conv     !Phase
            !write(*,*)i,j,tide%sbc(i,2*j),tide%sbc(i,2*j-1),ap(2*j)
         enddo
         if(i.ne.m)then
            if (mnproc==1) 
     &         write(lp,'(''tiderSB: Some mismatch in OBC indices'')')
            call xcstop ('(tiderSB)')
            stop '(tiderSB)'
         endif
      enddo
501   close(10)    
      if (mnproc==1)print '(a)',
     &   '####################################################'
c
c
c --- Load current data from file 
      tide%uv=ltideuv
      if (tide%uv) then
         fname_curr(1:7)='./Data/'
         fname_curr(8:10)=ctide(1:3)
         fname_curr(11:22)='obc_curr.dat'
         inquire(file=fname_curr,exist=ex)
         if (ex) then
            if (mnproc==1) print *,'Loading file: ',fname_curr
            OPEN(10,FILE=fname_curr,STATUS='OLD')
            DO i=1,5+ncf
               READ(10,'(a80)')dummy
            ENDDO
            DO i=1,mtot
               READ(10,'(i4,32f8.2)') m,ap_u
               DO j=1,nconst
                  tide%sbc_u(i,4*j-3) = ap_u(4*j-3)*v_cnv        !V_max
                  tide%sbc_u(i,4*j-2) = ap_u(4*j-2)*v_cnv        !V_min
                  tide%sbc_u(i,4*j-1) = ap_u(4*j-1)*phase_cnv    !Net_phase
                  tide%sbc_u(i,4*j)   = ap_u(4*j)  *phase_cnv    !Inclination
               ENDDO
               IF(i.NE.m)THEN
                  if (mnproc==1)
     &            WRITE(*,'(''tiderSB: Some mismatch in OBC indices'')')
                  call xcstop('(tiderSB)')
                  stop '(tiderSB)'
               ENDIF
            ENDDO
            close(10)
         else
            if (mnproc==1)print '(a)',
     &         'tiderSB: Tidal file does not exist: ',fname_curr
         endif
      else
         if (mnproc==1) 
     &      print '(a)','tiderSB: Tidal boundary currents are not used.'
      endif
      end subroutine tiderSB
c
c
c --- -------------------------------------------------------
c --- Author: K. Simonsen, June, 1999, neRSC
c --- -------------------------------------------------------
c --- Calculates tidal heights into array h_tide from:
c ---   freq:    tidal frequncies (rad/day)  (1)
c ---   amphaSB: elvation amp and Greenwich Phase Lags   (1)
c ---   time:    model time (days) 
c ---   astr_dt: Interval (days) between calculation of
c ---            astr. arg. See tideast.f for details.
c --- (1) see subroutine tiderSB for dimension and details
c --- on these data. 
c --- -------------------------------------------------------
      subroutine tide_hgt(time)
      implicit none
      real   , intent(in) :: time
CKAL  integer, intent(in) :: astr_time 
c --- astr_time contained in module, no need to pass as arg
c
      integer mtot,m,n,i
      real hgt,tide_time,floc(16)
      real ugt,vgt
      real, parameter :: pi=3.1415927
      logical, save :: lfirst=.true.
      real, save    :: tfirst
      real, save :: slowstart
      real v1,v2,vcos,vsin
c
      mtot=2*(itdm+jtdm)
c
c --- Tidal time
      if (tide%astro_calc) then
         tide_time = time-astr_time*1.0              !Time from midinterval
         do n=1,nconst
            floc(n)=tide%astro_arg(n)+tide_time*tide%freq(n)
         enddo
      else
         do n=1,nconst                      
           floc(n)=tide%freq(n)*time
         enddo
      endiF
c
c --- Slowstart option
      if (lfirst) then
         if (mnproc==1) print '(a)','tide_hgt: lfirst'
         if (lslowstart) then
            slowstart=0.0
         else
            slowstart=1.0
         endif
         lfirst=.false.
         tfirst=time
      endif
c
      if (lslowstart) then
         if (time-tfirst <= 1.0) then
            slowstart=sin((time-tfirst)*pi/2.0)
         else
            slowstart=1.0
         endif
      endif
c
c --- Add up heights at boundary points from all constituents
      do i=1,mtot
         hgt=0.
         ugt=0.
         vgt=0.
         do n=1,nconst
            hgt=hgt + tide%sbc(i,2*n-1)*COS(floc(n)-tide%sbc(i,2*n))
            if (tide%uv) then
               v1= tide%sbc_u(i,4*n-3)* COS(floc(n)+tide%sbc_u(i,4*n-1))
               v2= tide%sbc_u(i,4*n-2)* SIN(floc(n)+tide%sbc_u(i,4*n-1))
               vcos = COS(tide%sbc_u(i,4*n))
               vsin = SIN(tide%sbc_u(i,4*n))

               ugt=ugt + v1*vcos - v2 *vsin
               vgt=vgt + v1*vsin + v2 *vcos
            endif
         enddo
         tide%h(i)=hgt*slowstart
         if (tide%uv) then
            tide%u(i)=ugt*slowstart
            tide%v(i)=vgt*slowstart
         endif
      enddo
c
cdiag open(10,file='tide_hgtA.dat',status='unknown',position='append')
cdiag write(10,'(100f15.5)')time,tide_time,tide%h(5)/98060.,floc(1:nconst)
cdiag close(10)
cdiag open(10,file='tide_hgtB.dat',status='unknown',position='append')
cdiag write(10,'(100f15.5)')time,tide%sbc(5,1:2*nconst)
cdiag close(10)
cdiag open(10,file='tide_hgtC.dat',status='unknown',position='append')
cdiag write(10,'(100f15.5)')time,tide%astro_arg(1:nconst)
cdiag close(10)
c
      end subroutine tide_hgt


c --- --------------------------------------------------------------
c --- Calculate the astronomical argument at given day and year
c --- partly following 
c --- Manual for tidal heights analysis and predictionæ, by
c --- M. G. G. Foreman, Rep. 77-10, IOS, 1977
c ---
c --- --------------------------------------------------------------
CKAL  subroutine tide_ast(day,year)
CKAL  implicit none
CKAL  integer, intent(in) :: day  ! ordinal day rel. year
CKAL  integer, intent(in) :: year ! actual year
      subroutine tide_ast(dtime,yrflag)
      implicit none
      real*8, intent(in) :: dtime
      integer, intent(in):: yrflag
c
      integer year, day, ihour
      integer, parameter :: apn=16
      integer kd,kd0,k
      real h0,s0,p0,twopi,d1,d2
      real f
      real aarg(nastro)
c     
      integer icc,iyy,i
c
c --- Update astronomical time
      astr_time=int(dtime)
c
c --- Calculate year, ordinal day from dtime
      call forday(dtime, yrflag, year,day,ihour)
      day=day-1
c
      twopi = datan(1.D0)*8.D0
c
c --- ---------------------------- Gregorian Calender -----
c --- Partly adopted from the Foreman package 
c
c --- The time applied is based on the 'universal' Gregorian
c --- calender. A historical curiosum is that the Gregorian
c --- reform of the Julian calender omitted 10 days in 1582
c --- (the day after Oct. 4, 1582 became Oct. 15, 1582)
c --- in order to restore the date of the vernal equinox to
c --- March 21 and revised the leap year rule so that centurial 
c --- years not divisible by 400 were not leap years.
c
c --- Note that day kd=1 corresponds to Jan. 1, year 0000.!!!!!
c --- -----------------------------------------------------

      icc=year/100                       !Centuries 
      iyy=year-icc*100                   !Year in the century

      k=min(iyy,1)*MIN(icc-(icc/4)*4,1)  !Leap year correction
      kd=-k
      kd=kd+icc*36524 + (icc+3)/4        !Last day of last century
      kd=kd+iyy*365+(iyy+3)/4            !Last day in last year

                                         !Input day is day in actual year
      kd=kd+day                          !which corresponds to 
                                         !Gregorian day number kd.

c --- -----------------------------------------------------
c The formulae for the ephermies are  relative to 12.00 UT, 31.12.1899,
c which corresponds to gregorian day: kd0=693961 - 0.5
c --- 
      kd0=693961
      d1=float(kd-kd0)-0.5d0            !Gregorian days since kd0.
      if (mnproc==1)  then
         print '(a,2i5,i10)',
     &      'tide_ast: century, year, day (rel century) ',icc,iyy,kd
         print '(a,f14.2)','tide_ast: day (rel 31.12.1899) ',d1
      end if


c --- ------------------------------- Astronomical argument. -----
c --- Background:
c
c --- The astronomical argument:
c --- V = i tau + j S + k H + l P + m Enp + n Pp + semi 
c --- where (i,j,k,l,m,n) are the Doodson numbers, semi is a   
c --- phase correction , (S, H, P, Enp, Pp) are the ephermies:
c
c --- S: mean longitide of the moon
c --- H: mean longitide of the sun
c --- P: mean longitide of lunar perigee
c --- Enp: negative of the longitude of the mean ascending node
c --- PP: mean longitide of the solar perigee (perihelion)
c
c --- and finaly tau= hour/24 + H - S.
c --- Insertion of tau provides
c
c --- V = i hour/24 + (j-i)S + (k+i) + l P m Enp + n Pp  + semi
c
c --- which is the expression used  by Schwiderski(1980, 1986)
c --- for hour= 0 (midnight) (he did not incl. Enp and Pp),
c --- and also adopted in this rutine. The changes  in the
c --- astronomircal argument are relatively slow  and it is 
c --- common practice to calculate them only once for a given
c --- interval, Ti,  at the time of the midinter val, Tm, and  
c --- and then estimate the argument from:
c
c --- V(t) = V(Tm) + (t-Tm)sigma  for t=Tm-0.5 Ti; Tm+0.5 Ti
c
c --- and the tidal height ht(t)
c
c --- h(t) = SUM_i f_i(Tm) A_i COS{ V(t) - Grwpl_i  + u_i(Tm) }
c
c --- where
c ---   i           is constituent index.
c ---  A_i, Grepl_i are the amplitude and Greenwich Phase Lag 
c ---               for the give location
c ---  f_i,u_i      are nodel modulation and phase correction
c ---               due to influence from the satellite frequnecies
c ---               to time Tm. This is NOT included here.
c ---  V(t)         is the astronomical argument as defined  above              
c ---
c --- In this routine V(Tm) is calculated, while (t-Tm)sigma is 
c --- calculated in routine tide_hgt.f
c ---
c --- The Doodson numbers, phase correction (semi)  and number of
c --- satellite constituents (nj) for the eight main frequencies are:
c ---
c ---         i  j  k  l  m  n  semi  nj
c --- Q1      1 -2  0  1  0  0 -0.25  10
c --- O1      1 -1  0  0  0  0 -0.25   8
c --- P1      1  1 -2  0  0  0 -0.25   6
c --- K1      1  1  0  0  0  0 -0.75  10
c --- N2      2 -1  0  1  0  0  0.0    4
c --- M2      2  0  0  0  0  0  0.0    9
c --- S2      2  2 -2  0  0  0  0.0    3
c --- K2      2  2  0  0  0  0  0.0    5
c
c --- Note: semi is expressed in cycles
c 
c --- --- ----------------------------------------------------------------
c 
c --- The formulae for calculating this ephermies are from the
c --- Foreman package, who adopted this from 
c --- Explanatory Supplement to the Astronomical Ephermies
c --- and the American Ephermis and Nautical Almanac, 1961,
c --- pages 98 and 107. The day d1 is defined above.
c 
c --- Due to uncertainies in numerical model, and the
c --- relative small contribution of constituents not listed
c --- here, only these eight main constituents are included. 
c ---
c --- Since m,n= 0 for these 8 constituents Enp and Pp are
c --- not included in the implementation  of V.
c 
c --- ----------------------------------------------------------------
c
      d2=d1*1.d-4
      h0=279.696678d0+.9856473354d0*d1+.00002267d0*d2*d2

      s0=270.434164d0+13.1763965268d0*d1-.000085d0*d2*d2+
     .   .000000039d0*d2**3

      p0=334.329556d0+.1114040803d0*d1-.0007739d0*d2*d2-
     .   .00000026d0*d2**3
c
c
c     !pp=281.220833d0+.0000470684d0*d1+.0000339d0*d2*d2+.00000007d0*d2**3
c     !np=-259.183275d0+.0529539222d0*d1-.0001557d0*d2*d2-.00000005d0*d2**3
c
      f=360.d0
c
c --- Convert from degrees to cycles
      h0=h0/f
      s0=s0/f
      p0=p0/f
c     !pp=pp/f
c     !np=np/f
c
c --- ------------------------------ Astr. arguments at midneight (cycles)
      aarg(1)= (   h0 - 3.*s0 + p0 -.25)       !Q1
      aarg(2)= (   h0 - 2.*s0      -.25)       !O1
      aarg(3)= (  -h0              -.25)       !P1
      aarg(4)= (   h0              -.75)       !K1
      aarg(5)= (2.*h0 - 3.*s0 + p0     )       !N2
      aarg(6)= (2.*h0 - 2.*s0          )       !M2
      aarg(7)= (0.                     )       !S2
      aarg(8)= (2.*h0                  )       !K2
c
c --- Convert to rad [0;2pi]
      do i=1,nconst
         tide%astro_arg(i)=MOD(aarg(tide%iastro(i))*twopi,twopi)
      enddo
      tide%astro_arg(nconst+1) = 0.
c
      end subroutine tide_ast
            

c --  Diagnostic routine
      subroutine tide_diag(n,dtime)
      use mod_xc
      use mod_forcing_nersc
      use mod_hycom_nersc
      implicit none
      integer, intent(in)             :: n
      real, intent(in) :: dtime

      character(len=15) fname
      character(len= 3) tag3, cdd 
      character(len= 2) chh
      character(len= 4) cyy, css
      integer, save :: irec=0
      integer j,l,i,m, iday, ihour, iyear,iss
      logical, save :: lfirst=.true.
      real, dimension(itdm,jtdm) :: 
     . modlon,modlat
      integer, parameter :: maxpkt=1000
      integer, save :: npkt
      integer, save, dimension(maxpkt) :: ix,jx
      logical ex
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)  :: 
     &  slptile
      logical :: dumpboundry=.true.
      include 'common_blocks.h'
c
c --- Retrieve diag points
      if (lfirst) then
         ix=0; jx=0
         j=1
         inquire(file='tide_diag.in',exist=ex)
         if (ex) then
            open(10,file='tide_diag.in')
            do j=1,maxpkt
               read(10,*,end=100)ix(j),jx(j)
               ix(j)=max(1,min(itdm,ix(j)))
               jx(j)=max(1,min(jtdm,jx(j)))
               if (mnproc==1)print *,'tide:',j
            enddo
100         close(10)
            npkt=j-1
         else
            npkt=-1
         endif
c
c ---    Dump lonlat positions to diag. file
         call xcaget(modlon,plon,0)
         call xcaget(modlat,plat,0)
         if (npkt>0) then
            if (mnproc==1) then
               open(10,file='tidediag.dat')
               do j=1,npkt
                  write(10,'(i3.3,2f13.4)')j,
     &               modlon(ix(j),jx(j)),modlat(ix(j),jx(j))
               enddo
               close(10)
            end if
         end if
      endif ! first
c
      slptile=  slp(:,:,l0)*w0+slp(:,:,l1)*w1
     &         +slp(:,:,l2)*w2+slp(:,:,l3)*w3
c
c --- Dump some variables to diag file -
      if (npkt>0) then
      do l=1,npkt
         if (ix(l)==0) exit
         i=ix(l)
         j=jx(l)
c ---    Each MPI task is response for its local points
         if ( i>i0 .and. i<=i0+ii .and. j>j0 .and. j <=j0+jj) then
            write(tag3,'(i3.3)')l
            open(10,file=rungen//'tide_diag3_'//tag3//'.dat',
     .              status='unknown',position='append')
            write(10,'(9f15.5)') time,plon(i-i0,j-j0),plat(i-i0,j-j0),
     .         srfhgt(i-i0,j-j0)/onem, montg1(i-i0,j-j0)/(thref*onem), 
     .         pbavg(i-i0,j-j0,n)/onem,(slp0-slptile(i-i0,j-j0))/100., 
     .         ubavg(i-i0,j-j0,n),vbavg(i-i0,j-j0,n)
            close(10)
         end if
      enddo
      endif ! npkt > 0
c
c --- generate tecplot scatterplots of tidal height along boundary. Switched
c --- on by presence of tide_diag.in
      if (dumpboundry.and.npkt>0) then
      call forday(dtime, yrflag, iyear,iday,ihour)
      iday=iday-1
      iss=(dtime-floor(dtime))*86400.d0 -ihour*3600. 
      write(cyy,'(i4.4)') iyear
      write(cdd,'(i3.3)') iday
      write(chh,'(i2.2)') ihour
      write(css,'(i4.4)') iss
      call xcaget(modlon,plon,0)
      call xcaget(modlat,plat,0)
      if (mnproc==1.and.dumpboundry.and.npkt>0) then
         open(10,file=rungen//'tide_boundry_'//cyy//'_'//cdd//
     &      '_'//chh//'_'//css//'.tec')
         write(10,'(a)') 'TITLE="Tide boundary data '
     &      //cyy//cdd//chh//'"'
         write(10,'(a)') 'VARIABLES="i" "j" "lon" "lat" "tide_hgt"'
         write(10,'(a)') 'ZONE F=POINT'
         do i=1,itdm
            j=2
            m=i
            write(10,'(2i5,3e15.4)') i,j,modlon(i,j),modlat(i,j),
     &         tide%h(m)
         end do
         do i=1,itdm
            j=jtdm-1
            m=2*itdm+jtdm+1-i
            write(10,'(2i5,3e15.4)') i,j,modlon(i,j),modlat(i,j),
     &         tide%h(m)
         end do
         do j=1,jtdm
            i=itdm-1
            m=itdm+j
            write(10,'(2i5,3e15.4)') i,j,modlon(i,j),modlat(i,j),
     &         tide%h(m)
         end do
         do j=1,jtdm
            i=2
            m=2*(itdm+jtdm)+1-j
            write(10,'(2i5,3e15.4)') i,j,modlon(i,j),modlat(i,j),
     &         tide%h(m)
         end do
         close(10)
      end if ! mnproc==1
      end if ! dumpboundr & npkt>1
c
      lfirst=.false.
      end subroutine tide_diag
      end module mod_tides_nersc
