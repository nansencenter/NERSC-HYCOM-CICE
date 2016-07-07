      subroutine getdat(flnm,time,artype,initl,icegln,
     &                  iexpt,yrflag,kkin)
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      character        flnm*(*)
      double precision time(3)
      logical          initl,icegln
      integer          artype,iexpt,yrflag,kkin
c
c --- read model fields and extract portion of global fields.
c
      integer l
c
      character, allocatable :: util(:)*2
      real,      allocatable :: work(:,:)
c
      allocate( util(idm*jdm+14) )
      allocate( work(idm,jdm)    )
c
      l = len_trim(flnm)
c
      if     (flnm(l-1:l).eq.'.a' .or. flnm(l-1:l).eq.'.b') then
c ---   HYCOM 2.0 array I/O archive file.
        call getdata(flnm,time,artype,.true.,
     &               initl,icegln,iexpt,yrflag,kkin,      work)
      else
c ---   HYCOM 1.0 pakked archive file.
        call getdatp(flnm,time(1),
     &               initl,icegln,iexpt,yrflag,kkin, util,work)
        time(2) = time(1)
        time(3) = time(1)
        artype  = 1
      endif
      deallocate( util, work )
      return
      end

      subroutine getdatb(flnm,time,artype,initl,icegln,
     &                   iexpt,yrflag,kkin)
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      character        flnm*(*)
      double precision time(3)
      logical          initl,icegln
      integer          artype,iexpt,yrflag,kkin
c
c --- read model fields and extract portion of global fields.
c --- ignore ".a" vs ".b" range mismatches.
c
      integer   l
      character, allocatable :: util(:)*2
      real,      allocatable :: work(:,:)
c
      allocate( util(idm*jdm+14) )
      allocate( work(idm,jdm)    )
c
      l = len_trim(flnm)
c
      if     (flnm(l-1:l).eq.'.a' .or. flnm(l-1:l).eq.'.b') then
c ---   HYCOM 2.0 array I/O archive file.
        call getdata(flnm,time,artype,.false.,
     &               initl,icegln,iexpt,yrflag,kkin,      work)
      else
c ---   HYCOM 1.0 pakked archive file.
        call getdatp(flnm,time(1),
     &               initl,icegln,iexpt,yrflag,kkin, util,work)
        time(2) = time(1)
        time(3) = time(1)
        artype  = 1
      endif
      deallocate( util, work )
      return
      end

      subroutine getdatp(flnm,time,initl,icegln,
     &                   iexpt,yrflag,kkin, util,work)
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
      character        flnm*(*)
      double precision time
      logical          initl,icegln
      integer          iexpt,yrflag,kkin
      character        util(idm*jdm+14)*2
      real             work(idm,jdm)
c
c --- read model fields and extract portion of global fields.
c --- HYCOM 1.0 pakked archive file.
c --- (no time-averaged fluxes in this version)
c
      character cline*80
      character text*8,preambl(5)*79
      character cvarin*6,c2*2
      real      hminb,hmaxb,thet
      integer   i,iversn,ios,j,k,l,lgth
c
      integer ni
      data    ni/14/
c
      if (initl) then
        open (unit=ni,file=flnm,form='formatted',
     &        status='old',action='read',iostat=ios)
        if (ios.ne.0) then
          write(lp,*)
          write(lp,*) "error in getdat - can't open ",
     &                trim(flnm)
          write(lp,*) 'iostat = ',ios
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
      endif
c
c --- which version of file format?
      read(ni,'(a2)') c2
      rewind(unit=ni)
      if (c2.ne.'::') then
        read (ni,'(a80/a80/a80/a80)') ctitle
        write(lp,'(a80/a80/a80/a80)') ctitle
        read (ni,*) iversn,cvarin
        write(lp,*) cvarin,' = ',iversn
        if (cvarin.ne.'iversn') then
          write(lp,*)
          write(lp,*) 'error in getdat - input ',cvarin,
     &                          ' but should be iversn'
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        read (ni,*) iexpt,cvarin
        write(lp,*) cvarin,' = ',iexpt
        if (cvarin.ne.'iexpt ') then
          write(lp,*)
          write(lp,*) 'error in getdat - input ',cvarin,
     &                          ' but should be iexpt '
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        read (ni,*) yrflag,cvarin
        write(lp,*) cvarin,' = ',yrflag
        if (cvarin.ne.'yrflag') then
          write(lp,*)
          write(lp,*) 'error in getdat - input ',cvarin,
     &                          ' but should be yrflag'
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
      endif
c
      do 14 k=1,kk
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
 100  format (' ::',a8,i8,f8.1,f8.3,2i5,i9,1x,14a2)
      if     (text(1:8).ne.'u-vel.  ') then
        write(lp,*)
        write(lp,*) 'error in getdat - layer ',k,
     &             ' does not exist (kk= ',kk,')'
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
      call unpakk(  work(idm-i,1),idm,i,j,util,lgth)
      call extrct_u(work,idm,jdm,iorign,jorign, 
     &              u(1,1,2*k),ii,jj)
      write(lp,'("unpakk ",a," into ",a,i3)') text(1:8),'u       ',k
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work(1,2),idm,i,j,util,lgth)
      call extrct_v(work,idm,jdm,iorign,jorign, 
     &              v(1,1,2*k),ii,jj)
      write(lp,'("unpakk ",a," into ",a,i3)') text(1:8),'v       ',k
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              dp(1,1,k),ii,jj)
      write(lp,'("unpakk ",a," into ",a,i3)') text(1:8),'dp      ',k
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              temp(1,1,2*k),ii,jj)
      write(lp,'("unpakk ",a," into ",a,i3)') text(1:8),'temp    ',k
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              saln(1,1,2*k),ii,jj)
      write(lp,'("unpakk ",a," into ",a,i3)') text(1:8),'saln    ',k
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              th3d(1,1,2*k),ii,jj)
      write(lp,'("unpakk ",a," into ",a,i3)') text(1:8),'th3d    ',k
c
      if(ntracr.eq.1) then
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              trcr(1,1,2*k,1),ii,jj)
      write(lp,'("unpakk ",a," into ",a,2i3)') text(1:8),'trcr    ',k,1
      end if
c
ccc      write(lp,'(a,i4)') 'shown below: density in layer',k
ccc      call zebra(th3d(1,1,2*k),ii,ii1,jj1)
c
      write(lp,'(a,f9.5)') 'finished reading data for layer',thet
      call flush(lp)
      theta(k)=thet
 14   continue
      kkin=kk
c
c --- now unpack barotropic velocity field
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work(idm-i,1),idm,i,j,util,lgth)
      call extrct_u(work,idm,jdm,iorign,jorign, 
     &              ubaro,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'ubaro   '
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work(1,2),idm,i,j,util,lgth)
      call extrct_v(work,idm,jdm,iorign,jorign, 
     &              vbaro,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'vbaro   '
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              montg,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'montg   '
* --- discard montg
*     write(lp,'("unpakk ",a," into ",a)') text(1:8),'work    '
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              srfht,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'srfht   '
c
      surflx(:,:) = 0.0
      salflx(:,:) = 0.0
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              dpbl,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'dpbl    '
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              dpmixl,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'dpmixl  '
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              tmix,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'tmix    '
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              smix,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'smix    '
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work,idm,i,j,util,lgth)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              thmix,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'thmix   '
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work(idm-i,1),idm,i,j,util,lgth)
      call extrct_u(work,idm,jdm,iorign,jorign, 
     &              umix,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'umix    '
c
      read (ni,'(3x,a,a)',end=6) text,cline
      read (cline,      *,end=6) nstep,time,thet,i,j,lgth
      read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
      write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,(util(l),l=1,14)
      call unpakk(  work(1,2),idm,i,j,util,lgth)
      call extrct_v(work,idm,jdm,iorign,jorign, 
     &              vmix,ii,jj)
      write(lp,'("unpakk ",a," into ",a)') text(1:8),'vmix    '
c
c --- is there ice?
      read (ni,'(3x,a,a)',iostat=ios) text,cline
      icegln = ios.eq.0 .and. text.eq.'covice  '
      if     (icegln) then
        read (cline,      *,end=6) nstep,time,thet,i,j,lgth
        read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
        write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,
     &                          (util(l),l=1,14)
        call unpakk(  work,idm,i,j,util,lgth)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                covice,ii,jj)
c
        read (ni,'(3x,a,a)',end=6) text,cline
        read (cline,      *,end=6) nstep,time,thet,i,j,lgth
        read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
        write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,
     &                          (util(l),l=1,14)
        call unpakk(  work,idm,i,j,util,lgth)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                thkice,ii,jj)
c
        read (ni,'(3x,a,a)',end=6) text,cline
        read (cline,      *,end=6) nstep,time,thet,i,j,lgth
        read (ni,  '(40a2)',end=6) (util(l),l=1,lgth)
        write(lp,100) text(1:8),nstep,time,thet,i,j,lgth,
     &                          (util(l),l=1,14)
        call unpakk(  work,idm,i,j,util,lgth)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                temice,ii,jj)
      else
        covice(:,:) = 0.0
        thkice(:,:) = 0.0
        temice(:,:) = 0.0
      endif
c
      if (initl) then
c ---   acquire basin depths and land/sea mask
        call getdepth(work)
      end if
c
*     write(lp,'(a)') 'shown below: sea surface height'
*     call zebra(srfht,ii,ii1,jj1)
*     call flush(lp)
c
      return
c
c --- unexpected end of file
 6    continue
      write (lp,*) '***** unexpected end of archive file *****'
      call flush(lp)
      call clsgks
      stop '(e-o-f)'
      end

      subroutine getdata(flnm,time,artype,lrange,initl,icegln,
     &                   iexpt,yrflag,kkin, work)
      use mod_plot  ! HYCOM plot array interface
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      character        flnm*(*)
      double precision time(3)
      logical          lrange,initl,icegln
      integer          artype,iexpt,yrflag,kkin
      real             work(idm,jdm)
c
c --- read model fields and extract portion of global fields.
c --- HYCOM 2.0 array I/O archive file.
c --- (no time-averaged fluxes in this version)
c
      character cline*80
      character preambl(5)*79
      character cvarin*6
      real      hminb,hmaxb,thet,thbase
      integer   i,idmtst,j,jdmtst,iversn,ios,k,l,layer,sigver
      integer   itr,ktr,ntr
      logical   lke,nodens
c
      integer ni
      data    ni/14/
c
      if (initl) then
        l = len_trim(flnm)
        call zaiopf(flnm(1:l-2)//'.a','old', ni)
        open (unit=ni,file=flnm(1:l-2)//'.b',form='formatted',
     &        status='old',action='read',iostat=ios)
        if (ios.ne.0) then
          write(lp,*)
          write(lp,*) "error in getdat - can't open ",
     &                flnm(1:l-2)//'.b'
          write(lp,*) 'iostat = ',ios
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
      endif
c
      read( ni,'(a80/a80/a80/a80)') ctitle
      write(lp,'(a80/a80/a80/a80)') ctitle
      read( ni,*) iversn,cvarin
      write(lp,*) cvarin,' = ',iversn
      if (cvarin.ne.'iversn') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                          ' but should be iversn'
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
      read( ni,*) iexpt,cvarin
      write(lp,*) cvarin,' = ',iexpt
      if (cvarin.eq.'jexpt ') then
c ---   skip jexpt, artype of 4 set below
        read( ni,*) iexpt,cvarin
        write(lp,*) cvarin,' = ',iexpt
      endif
      if (cvarin.ne.'iexpt ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be iexpt '
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
      read( ni,*) yrflag,cvarin
      write(lp,*) cvarin,' = ',yrflag
      if (cvarin.ne.'yrflag') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be yrflag'
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
      read( ni,*) idmtst,cvarin
      write(lp,*) cvarin,' = ',idmtst
      if (cvarin.ne.'idm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be idm   '
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
      read( ni,*) jdmtst,cvarin
      write(lp,*) cvarin,' = ',jdmtst
      if (cvarin.ne.'jdm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be jdm   '
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
c
      if (idmtst.ne.idm .or. jdmtst.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error in getdat - input idm,jdm',
     &                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',idm,   jdm,   '  (REGION.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
c
c --- artype==1 for normal archive files
c --- artype==2 for   mean archive files
c --- artype==3 for stddev archive files
c --- artype==4 for   diff archive files
c
      read( ni,'(a)') cline
      write(lp,'(a)') cline(1:len_trim(cline))
      if     (cline(25:28).eq.'mean') then
        artype = 2
      elseif (cline(25:28).eq.'std.') then
        artype = 3
      elseif (cline(25:28).eq.'diff') then
        artype = 4
      else
        artype = 1
      endif
      write(lp,'(a,i2)') 'artype =',artype
c
      if     (artype.eq.3 .or. artype.eq.4) then
c ---   always have thickness (deviation) as the 1st tracer, if any
        itr = min(2,ntracr+1)
      else
        itr = 1
      endif
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(1),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
      call extrct_p(work,idm,jdm,iorign,jorign,
     &              montg,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'montg   '
* --- discard montg
*     write(lp,'("input  ",a," into ",a)') cline(1:8),'work    '
c
c --- detect version 2.2 normal archive files
      nodens = layer.ne.0
      if     (nodens) then
        sigver = layer
        thbase = thet
      endif
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(2),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
      call extrct_p(work,idm,jdm,iorign,jorign,
     &              srfht,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'srfht   '
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
      if     (cline(1:8).ne.'steric  ') then
        steric(:,:) = 0.0
      else
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                steric,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'steric  '
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
      endif
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              surflx,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'surflx  '
c
      call time_hour(time)  !reset, assuming time is on the hour
      write(lp,*) 'time3 = ',time
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              salflx,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'salflx  '
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              dpbl,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'dpbl    '
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              dpmixl,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'dpmixl  '
c
      if     (.not. nodens) then
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                tmix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'tmix    '
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                smix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'smix    '
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                thmix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'thmix   '
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .true. ,lrange)
        call extrct_u(work,idm,jdm,iorign,jorign, 
     &                umix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'umix    '
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .true. ,lrange)
        call extrct_v(work,idm,jdm,iorign,jorign, 
     &                vmix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'vmix    '
      endif !.not. nodens
c
c --- is there ke?
      read (ni,'(a)',iostat=ios) cline
      write(lp,'(a)')            cline(1:len_trim(cline))
      lke = artype.gt.1 .and. ios.eq.0 .and. cline(1:8).eq.'kemix   '
      if     (lke) then  ! mean or std. or diff archive with ke
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .true. ,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                kemix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'kemix   '
        read (ni,'(a)',iostat=ios) cline
        write(lp,'(a)')            cline(1:len_trim(cline))
      endif
c
c --- is there ice?
*     read (ni,'(a)',iostat=ios) cline
*     write(lp,'(a)')            cline(1:len_trim(cline))
      icegln = ios.eq.0 .and. cline(1:8).eq.'covice  '
      if     (icegln) then
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                covice,ii,jj)
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                thkice,ii,jj)
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                temice,ii,jj)
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
      else
        covice(:,:) = 0.0
        thkice(:,:) = 0.0
        temice(:,:) = 0.0
      endif
c
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .true. ,lrange)
      call extrct_u(work,idm,jdm,iorign,jorign, 
     &              ubaro,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'ubaro   '
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .true. ,lrange)
      call extrct_v(work,idm,jdm,iorign,jorign, 
     &              vbaro,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'vbaro   '
c
      if     (lke) then  ! mean or std. archive with ke
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .true. ,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                kebaro,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'kebaro  '
      endif
c
      kkin=1
      do 14 k=1,kk
      if     (k.eq.2) then
c ---   already input at end of k=1 loop.
      else
        read (ni,'(a)',end=6)   cline
        write(lp,'(a)')         cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      if     (cline(1:8).ne.'u-vel.  ') then
        write(lp,*)
        write(lp,*) 'error in getdat - layer ',k,
     &             ' does not exist (kk= ',kk,')'
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
      call getfld(  work, ni, hminb,hmaxb, .true. ,lrange)
      call extrct_u(work,idm,jdm,iorign,jorign, 
     &              u(1,1,2*k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'u       ',k
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .true. ,lrange)
      call extrct_v(work,idm,jdm,iorign,jorign, 
     &              v(1,1,2*k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'v       ',k
c
      if     (lke) then  ! mean or std. archive with ke
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                ke(1,1,2*k),ii,jj)
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'ke      ',k
      endif
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .true., lrange)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              dp(1,1,k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dp      ',k
      if     (cline(1:8).eq.'mnthknss') then
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .true., lrange)
        if     (itr.eq.1) then
c ---     discard std. thickness.
          write(lp,'("skip   ",a)') cline(1:8)
        else
          call extrct_p(work,idm,jdm,iorign,jorign, 
     &                  trcr(1,1,2*k,1),ii,jj)
          write(lp,'("input  ",a," into ",a,2i3)')
     &      cline(1:8),'trcr    ',k,1
        endif
      endif
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              temp(1,1,2*k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'temp    ',k
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
      call extrct_p(work,idm,jdm,iorign,jorign, 
     &              saln(1,1,2*k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'saln    ',k
c
      if     (.not. nodens) then
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
        call extrct_p(work,idm,jdm,iorign,jorign, 
     &                th3d(1,1,2*k),ii,jj)
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'th3d    ',k
      else
        call th3d_p(temp(1,1,2*k),saln(1,1,2*k),
     &              th3d(1,1,2*k),ii,jj, sigver,thbase)
        write(lp,'("    ",a8,"calculate ",a,i3)') " ",'th3d    ',k
        if     (k.eq.1) then
           tmix(:,:) = temp(:,:,2)
           smix(:,:) = saln(:,:,2)
          thmix(:,:) = th3d(:,:,2)
           umix(:,:) =    u(:,:,2)
           vmix(:,:) =    v(:,:,2)
          write(lp,'("copy   ",a," into ",a)') 'temp.1  ','tmix    '
          write(lp,'("copy   ",a," into ",a)') 'saln.1  ','smix    '
          write(lp,'("copy   ",a," into ",a)') 'th3d.1  ','thmix   '
          write(lp,'("copy   ",a," into ",a)') '   u.1  ','umix    '
          write(lp,'("copy   ",a," into ",a)') '   v.1  ','vmix    '
        endif !k==1
      endif !.not.nodens:else
c
c --- tracers and q2,q2l and and visc/diff, 
c ---  may be more in archive than plotted.
c
      if     (k.eq.1) then
        do ktr= itr,999
          read (ni,'(a)',iostat=ios) cline
          write(lp,'(a)')            cline(1:len_trim(cline))
          if (ios.ne.0) then
            write(lp,'(a,f9.5)') 'finished reading data for layer',thet
            call flush(lp)
            theta(k)=thet
            goto 114  ! archive containing only 1 layer
          elseif (cline(1:8).ne.'tracer  ' .and.
     &            cline(1:8).ne.'q2      ' .and.
     &            cline(1:8).ne.'q2l     ' .and.
     &            cline(1:8).ne.'viscty  ' .and.
     &            cline(1:8).ne.'t-diff  ' .and.
     &            cline(1:8).ne.'s-diff  '      ) then
            exit !end of tracers and visc/diff
          else
            i = index(cline,'=')
            read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
            call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
            if     (ktr.le.ntracr) then
              call extrct_p(work,idm,jdm,iorign,jorign, 
     &                      trcr(1,1,2*k,ktr),ii,jj)
              if     (cline(1:8).eq.'viscty  ' .or.
     &                cline(1:8).eq.'t-diff  ' .or.
     &                cline(1:8).eq.'s-diff  '     ) then
                do j= 1,jj
                  do i= 1,ii
                    if     (trcr(i,j,2*k,ktr).lt.2.0**99) then
                      trcr(i,j,2*k,ktr) = 1.e4*trcr(i,j,2*k,ktr)  !cm^2/s
                    endif
                  enddo
                enddo
              endif !visc/diff
              if     (cline(1:8).eq.'q2      ' .or.
     &                cline(1:8).eq.'q2l     '     ) then
                do j= 1,jj
                  do i= 1,ii
                    if     (trcr(i,j,2*k,ktr).lt.2.0**99) then
                      trcr(i,j,2*k,ktr) = 1.e7*trcr(i,j,2*k,ktr)  
                    endif
                  enddo
                enddo
              endif !q2/q2l
              write(lp,'("input  ",a," into ",a,2i3)')
     &          cline(1:8),'trcr    ',k,ktr
            endif
          endif
        enddo !ktr
        ntr=ktr-1
        if     (ntracr.gt.ntr) then
          write(lp,*)
          write(lp,*) 'error in getdat - fewer tracers than requested'
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
      else !k.gt.1
        do ktr= itr,ntr
          read (ni,'(a)',end=6) cline
          write(lp,'(a)')       cline(1:len_trim(cline))
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
          call getfld(  work, ni, hminb,hmaxb, .false.,lrange)
          if     (ktr.le.ntracr) then
            call extrct_p(work,idm,jdm,iorign,jorign, 
     &                    trcr(1,1,2*k,ktr),ii,jj)
            if     (cline(1:8).ne.'tracer  ') then !visc/diff
              do j= 1,jj
                do i= 1,ii
                  if     (trcr(i,j,2*k,ktr).lt.2.0**99) then
                    trcr(i,j,2*k,ktr) = 1.e4*trcr(i,j,2*k,ktr)  !cm^2/s
                  endif
                enddo
              enddo
            endif !visc/diff
            write(lp,'("input  ",a," into ",a,2i3)')
     &        cline(1:8),'trcr    ',k,ktr
          endif
        enddo !ktr
      endif !tracers+visc/diff
c
ccc      write(lp,'(a,i4)') 'shown below: density in layer',k
ccc      call zebra(th3d(1,1,2*k),ii,ii1,jj1)
c
      write(lp,'(a,f9.5)') 'finished reading data for layer',thet
      call flush(lp)
      theta(k)=thet
 14   continue
      kkin=kk
114   continue
c
      call time_hour(time)  !reset, assuming time is on the hour
c
      if (initl) then
c ---   acquire basin depths and land/sea mask
        call getdepth(work)
      end if
c
*     write(lp,'(a)') 'shown below: sea surface height'
*     call zebra(srfht,ii,ii1,jj1)
*     call flush(lp)
c
      return
c
c --- unexpected end of file
 6    continue
      write (lp,*) '***** unexpected end of archive file *****'
      call flush(lp)
      call clsgks
      stop '(e-o-f)'
      end

      subroutine getfld(work, iunit, hminb,hmaxb, lzero,lrange)
      use mod_za ! HYCOM array I/O interface
      implicit none
c
c --- read a single array
c
      logical lzero,lrange
      integer iunit
      real    work(idm,jdm), hminb,hmaxb
c
      integer mask(1)  !dummy which is never accessed
      integer i,j
      real    hmina,hmaxa
c
      call zaiord(work,mask,.false., hmina,hmaxa, iunit)
c
      if     (lrange) then
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          call flush(lp)
          call clsgks
          stop
        endif
      endif
c
      if     (lzero) then
        do j= 1,jdm
          do i= 1,idm
            if     (work(i,j).gt.2.0**99) then
              work(i,j) = 0.0
            endif
          enddo
        enddo
      endif
      return
      end

      subroutine time_hour(time)
      implicit none
c
      double precision time(3)
c
c --- reset time to an exact hour if very close to an hour.
c
      integer k
      double precision day,hour,ihr
c
      do k= 1,3
        day  = int(time(k))
        hour = (time(k)-day)*24.d0
        ihr  = nint(hour)
        if     (abs(hour-ihr).le.0.15d0) then
          time(k) = day + ihr/24.d0
        endif
      enddo
      return
      end

      subroutine th3d_p(temp,saln,th3d,no,mo,sigver,thbase)
      implicit none
c
      integer no,mo,sigver
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the appropriate equation of state.
c
      if     (sigver.eq.1) then
        call th3d_p1(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.2) then
        call th3d_p2(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.3) then
        call th3d_p3(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.4) then
        call th3d_p4(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.5) then
        call th3d_p5(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.6) then
        call th3d_p6(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.7) then
        call th3d_p7(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.8) then
        call th3d_p8(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.46) then
        call th3d_p46(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.48) then
        call th3d_p48(temp,saln,th3d,no,mo,thbase)
      else  !unknown
        th3d(:,:) = 0.0
      endif
      return
      end
      subroutine th3d_p1(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA0_7term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p2(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA2_7term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p3(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA0_9term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p4(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA2_9term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p5(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA0_17term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p6(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA2_17term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p7(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA0_12term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p8(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA2_12term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p46(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA4_17term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p48(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA4_12term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
