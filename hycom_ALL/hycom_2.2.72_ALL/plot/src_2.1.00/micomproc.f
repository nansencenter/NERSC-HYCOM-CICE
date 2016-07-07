      program hycomproc
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- hycom/micom processor (MKS real-basin version)
c
c --- use environment variable CROSS_LABELS to provide the
c --- name of a file containing a 4-character label for each
c --- layer on cross-section plots.  Note that CROSS_LABELS="NONE"
c --- indicates no cross section plots are needed (i.e. no[ij]sec=0),
c --- and can speed up the processing in such cases.
c
      parameter (kout=80)
c
      real, allocatable, dimension (:,:,:,:) ::
     &   trtr
      real, allocatable, dimension (:,:,:) ::
     &   utr,vtr,wtr,ttr,str, w,puv,xyij, thtr,
     &   tr2d
      real, allocatable, dimension (:,:)   ::
     &   uflux,vflux,vort,strmf,
     &   util1,depth1,depthu,depthv,dpdx,dpdy, work,
     &   uv2d,vv2d,wv2d,th2d,tm2d,sl2d, p2d
      real, allocatable, dimension (:)     ::
     &   coord,pout,dpbl2d,dpml2d,xlonlat
      real ubi,ubmi,vbi,vbim1,strmft,strmfu
c
      common/conrng/ amn,amx
c
      character label*81,cmonth(12)*3,text*18,flnm*240,region*8
      character cline*80
      character crflnm*240,crlabl(99)*4
      character ctrc_title(99)*8,ctrc_units(99)*8,ctrc_format(99)*8,
     &          ctrc_name(99)*56
      logical plotuv,plotw,plotnv,plotth,plotem,plosal,plotbl,plotml
     &       ,smooth,mthin,mdens,initl,lhycom,icegln,lperiod,lsecanom
     &       ,baclin
      integer lrfi(6),lrfo(6)
c
      character        csec_name(105)*56
      real             qqsec(105),qcsec(105),qqtrc(99),qqctrc(99)
      double precision sumsec(105),samsec(105)
c
      integer          artype,iexpt,kkin,yrflag,yrflgi,kpalet,mxlflg
      double precision dt,dt0,time3(3),time,year
c
c --- 'lhycom' -- hycom (vs micom) input file
      data lhycom/.false./
c
      data lrfi/6*0/,lrfo/6*0/
      data cmonth/'Jan','Feb','Mar','Apr','May','Jun',
     &            'Jul','Aug','Sep','Oct','Nov','Dec'/
c
      data tenm/10./,onem/1./,tencm/.1/,onecm/.01/,onemm/.001/
c
      data initl /.true. /
      data thref/1.e-3/,spcifh/3990./
      data flag/-.03125/,nquad/0/,ncount/0/
      character blank*40
      data blank/'                                        '/
      common/perframe/nperfr,locbar
      real dudxdn,dudxup,dvdydn,dvdyup
c
c --- color options
c --- ipalet = 0      --  contour lines only, no color
c --- ipalet = 1      --  alternate pastel shading of contour intervals
c --- ipalet = 2      --  use canonical sst color palette  ( 64 intervals)
c --- ipalet = 3      --  use rainer's gaudy color palette (100 intervals)
c --- ipalet = 4      --  two-tone shades                  ( 64 intervals)
c --- ipalet = 5      --  NRL's 100 false color palette    (100 intervals)
c --- ipalet = 6      --  NRL's 100 inverted fc palette    (100 intervals)
c --- ipalet = 7      --  NRL's  20 false color palette    ( 20 intervals)
c --- ipalet = 8      --  NRL's  20 fc, 2xwhite palette    ( 20 intervals)
      common /colopt/ ipalet,nbase,ibase(2,10)
c
c --- number of contour intervals for each palette
      real cntrs(-2:8)
      data cntrs / 16.0, 16.0, 16.0,
     &             16.0, 64.0, 100.0, 64.0, 100.0, 100.0, 20.0, 20.0 /
c
      call xcspmd
      call zaiost
      lp=6
      film=onemm
c
c --- establish time limits for plotting
c
ccc      open (unit=9,file='limits',status='old',readonly)
ccc      read (9,*) yearlo,yearhi,days
ccc      close (unit=9)
c
c --- read model data
c
c ---   'flnm  ' = name of file containing the actual data
c ---   'region' = name of model region (e.g. ATLa2.00)
c ---   'iexpt ' = experiment number x10  (000=from archive file)
c ---   'yrflag' = days in year flag (-1=none,0=360,1=366J16,2=366J01,3=actual)
c ---   'ntracr' = number of tracers (to plot, optional with default 0)
c ---   one line per tracer: 8-letter title, 8-letter format, long name
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'kdm   ' = number of layers
        read (*,'(a)') flnm
        write (lp,'(2a)') 'input file: ',flnm(1:len_trim(flnm))
        call flush(lp)
        read (*,'(a)') region
        write (lp,'(2a)') '    region: ',region
        call flush(lp)
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflgi,'yrflag')
        call blkini2(i,j,  'ntracr','idm   ')  !read ntracr or idm
        if (j.eq.1) then
          ntracr = i
          do ktr= 1,ntracr
            read(*,'(a8,a8,a8,a)')
     &        ctrc_title(ktr),ctrc_format(ktr),ctrc_units(ktr),
     &        ctrc_name(ktr)
            write (lp,'(2x,i2,3a)') 
     &        ktr,' title  = "',trim(ctrc_title( ktr)),'"',
     &        ktr,' format = "',trim(ctrc_format(ktr)),'"',
     &        ktr,' units  = "',trim(ctrc_units( ktr)),'"',
     &        ktr,' l.name = "',trim(ctrc_name(  ktr)),'"'
            call flush(lp)
          enddo
          call blkini(ii,  'idm   ')
        else
          ntracr = 0
          ii     = i
        endif
        call blkini(jj,    'jdm   ')
        call blkini(kk,    'kdm   ')
        if     (ii.ne.idm .or. jj.ne.jdm) then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                           idm,jdm,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
c
c ---   'thbase' = reference density (sigma units)
        call blkinr(thbase,
     &             'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c ---   'nperfr' = number of horizontal plots per frame
c ---   'lalolb' = spacing of latitude/longitude labels (<0 use array spacing)
c ---   'lalogr' = spacing of latitude/longitude grid over land (<0 land+sea)
c ---      (abs(lalogr)>1000: spacing is (abs(lalogr)-1000)/100.0)
c ---   'loclab' = flag indicating the location of the contour label
c ---      (1=upper-right,2=lower-right,3=lower-left,4=upper-left)
c ---   'locbar' = flag indicating the location of the color bar
c ---      (vertical:   10=right, 11=ur,12=lr,13=ll,14=ul,15=cr,16=cl)
c ---      (horizontal: 20=bottom,21=ur,22=lr,23=ll,24=ul,25=ct,26=cb)
c ---   'kpalet' = palete (0=none,1=pastel,2=sst,3=gaudy,4=2tone,5=fc,6=ifc)
c ---              paletes 2-8 require input of the central contour
c ---   'smooth' = smooth fields before plotting
c ---   'mthin ' = mask thin layers from plots (0=F,1=T,2=T&mdens=T)
c ---   'baclin' = plot baroclinic velocity (0=total:DEFAULT,1=baroclinic)
c ---   'i_th  ' = draw only every i_th vector in every (i_th/2) row
        call blkini(nperfr,'nperfr')
        call blkini(lalolb,'lalolb')
        call blkini(lalogr,'lalogr')
        call blkini(loclab,'loclab')
        call blkini(locbar,'locbar')
        call blkini(kpalet,'kpalet')
        call blkinl(smooth,'smooth')
        call blkini(ithin, 'mthin ')
        mthin = ithin.ge.1  !mask thin layers
        mdens = ithin.ge.2  !mask non-iospycnal layers
        call blkini2(i,j,  'baclin','i_th  ')  !read baclin or i_th
        if (j.eq.1) then
          baclin = i.eq.1  !baroclinic, vs total, velocity
          call blkini(i_th,  'i_th  ')
        else
          baclin = .false. !total velocity
          i_th   = i
        endif
c
        if     (kpalet.lt.0 .or. kpalet.gt.8) then
          write(lp,*)
          write(lp,*) 'error - illegal kpalet'
          write(lp,*)
          stop
        endif
c
c ---   'iorign' = i-origin of plotted subregion
c ---   'jorign' = j-origin of plotted subregion
c ---   'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
c ---   'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
        call blkini(iorign,'iorign')
        call blkini(jorign,'jorign')
        call blkini(ii,    'idmp  ')
        call blkini(jj,    'jdmp  ')
        if     (ii.eq.0) then
          ii=idm
        endif
        if     (jj.eq.0) then
          jj=jdm
        endif
c ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
c ---   from the full history grid (dimensioned idm x jdm). 
c ---   The size of the subgrid is determined by ii,jj.
        write (lp,'(2(a,i5),9x,2(a,i5))') 'extracting i =',iorign,
     &    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-1
        call flush(lp)
        if     (iorign.lt.1 .or. jorign.lt.1) then
          write(lp,*)
          write(lp,*) 'error - [ij]orign must be positive'
          write(lp,*)
          stop
        endif
c
c --- array allocation
c
      call plot_alloc
c
      iout=max(ii,jj)
      iwrk=max(iout,ii)
      jwrk=max(kout,jj)
c
      allocate(   trtr(ii,jj,kout,ntracr) )
c
      allocate(    utr(ii,jj,kout) )
      allocate(    vtr(ii,jj,kout) )
      allocate(    wtr(ii,jj,kout) )
      allocate(    ttr(ii,jj,kout) )
      allocate(    str(ii,jj,kout) )
      allocate(   thtr(ii,jj,kout) )
      allocate(      w(ii,jj,kk*2) )
      allocate(    puv(ii,jj,kk*2) )
      allocate(   xyij(ii,jj,2)    )
      allocate(   tr2d(iout,kout,ntracr) )
c
      allocate(  strmf(ii,jj) )
      allocate(  uflux(ii,jj) )
      allocate(  vflux(ii,jj) )
      allocate(   vort(ii,jj) )
      allocate(  util1(ii,jj) )
      allocate( depth1(ii,jj) )
      allocate( depthu(ii,jj) )
      allocate( depthv(ii,jj) )
      allocate(   dpdx(ii,jj) )
      allocate(   dpdy(ii,jj) )
      allocate(   uv2d(iout,kout) )
      allocate(   vv2d(iout,kout) )
      allocate(   wv2d(iout,kout) )
      allocate(   th2d(iout,kout) )
      allocate(   tm2d(iout,kout) )
      allocate(   sl2d(iout,kout) )
      allocate(   work(iwrk,jwrk) )
      allocate(    p2d(iout,kk+1) )
c
      allocate(  dpbl2d(iout) )
      allocate(  dpml2d(iout) )
      allocate( xlonlat(iout) )
c
      allocate(    pout(kout) )
c
      allocate(   coord(2*kk) )
c
      do k=1,2*kk
        coord(k)=float(k)
      enddo
c
         p2d(:,:) = flag
      dpbl2d(:)   = flag
      dpml2d(:)   = flag
c
      dpthfil = 'regional.depth'
c
      if     (nperfr.eq.1 .or. nquad.lt.0) then
        siz=.46
        csn=-0.9
        csb=-1.2
      elseif (nperfr.eq.2) then
        siz=.43
        csn=-0.9
        csb=-1.2
      elseif (nperfr.eq.4) then
        siz=.35
        csn=-0.9*0.7
        csb=-1.2*0.7
      else
        siz=.25
        csn=-0.9*0.4
        csb=-1.2*0.4
      endif
      chrsiz=.009*float(ii-1)/
     &       (2.*siz*min(1.,float(ii-1)/float(jj-1)))
      write(lp,*) 'chrsiz = ',chrsiz,ii/chrsiz,jj/chrsiz
      if     (loclab.eq.1) then
c ---   contour stats in upper right of horizontal plot
        xlab0 = ii-1-10.0*chrsiz
        xlab1 = ii-1-10.0*chrsiz
        xlab2 = ii-1-10.0*chrsiz
        ylab0 = jj-1- 3.0*chrsiz
        ylab1 = jj-1- 6.0*chrsiz
        ylab2 = jj-1- 9.0*chrsiz
      elseif (loclab.eq.2) then
c ---   contour stats in lower right of horizontal plot
        xlab0 = ii-1-10.0*chrsiz
        xlab1 = ii-1-10.0*chrsiz
        xlab2 = ii-1-10.0*chrsiz
        ylab0 =    1 + 9.0*chrsiz
        ylab1 =    1 + 5.0*chrsiz
        ylab2 =    1 + 3.0*chrsiz
      elseif (loclab.eq.3) then
c ---   contour stats in lower left  of horizontal plot
        xlab0 =    1 +10.0*chrsiz
        xlab1 =    1 +10.0*chrsiz
        xlab2 =    1 +10.0*chrsiz
        ylab0 =    1 + 9.0*chrsiz
        ylab1 =    1 + 5.0*chrsiz
        ylab2 =    1 + 3.0*chrsiz
      elseif (loclab.eq.4) then
c ---   contour stats in upper left  of horizontal plot
        xlab0 =    1 +10.0*chrsiz
        xlab1 =    1 +10.0*chrsiz
        xlab2 =    1 +10.0*chrsiz
        ylab0 = jj-1- 3.0*chrsiz
        ylab1 = jj-1- 6.0*chrsiz
        ylab2 = jj-1- 9.0*chrsiz
      else
        write(lp,*)
        write(lp,*) 'error - unknown loclab = ',loclab
        write(lp,*)
        call flush(lp)
        stop
      endif
      write (lp,'(a,i3,a)') '... plotting',nperfr,' maps per frame'
      call flush(lp)
c
      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.
        enddo
      enddo
c
c --- read the archive file.
c
      if (lhycom) then
        call getdat(flnm,time3,artype,initl,icegln,
     &              iexpt,yrflag,kkin)              ! hycom input
        time = time3(3)
      else
        call getdtm(flnm,time,initl, thbase)        ! micom input
        yrflag = max(yrflgi,0)
        artype = 1
      endif
c
      if     (artype.eq.3) then
        smooth = .false.
        mthin  = .false.
      endif
c
      if     (yrflag.eq.0) then
        year  = 360.0d0
      elseif (yrflag.lt.3) then
        year  = 366.0d0
      else
        year  = 365.25d0
      endif
c
*     do k=2,2*kkin,2
*     write(51,9999) 39,34,k/2,100.0*v(39,34,k),100.0*vbaro(39,34)
*9999 format(3i5,1p,2e14.6)
*     end do
c
      write(lp,'(/a,2f8.2/a,2f8.2)') 
     &     'sub-domain longitude range = ',
     &    minval(plon(:,:)),maxval(plon(:,:)),
     &     'sub-domain latitude  range = ',
     &    minval(plat(:,:)),maxval(plat(:,:))
c
      lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
      if     (lperiod) then
        write(lp,'(/a/)') 'sub-domain assumed to be periodic'
      else
        write(lp,'(/a/)') 'sub-domain assumed to be non-periodic'
      endif
c
      if     (.not.lperiod .and. iorign+ii-1.gt.idm) then
        write(lp,*)
        write(lp,*) 'error - sub-domain not inside the full domain'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
      call bigrid(depths)
c
c --- check that bathymetry is consistent with this archive.
c
      ibadl = 0
      ibads = 0
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            if     (srfht(i,j).gt.2.0**99) then
              ibads = ibads + 1   ! topo sea, srfht land
*             if     (mod(ibads,100).eq.1) then
*               write(lp,*) 'topo sea, srfht land at i,j = ',i,j
*             endif
            endif
          else
            if     (srfht(i,j).lt.2.0**99) then
              ibadl = ibadl + 1   ! topo land, srfht sea
*             if     (mod(ibadl,100).eq.1) then
*               write(lp,*) 'topo land, srfht sea at i,j = ',i,j
*    &                      ,srfht(i,j)
*             endif
            endif
          endif
        enddo
      enddo
      if     (ibads.ne.0) then
        write(lp,*)
        write(lp,*) 'error - wrong bathymetry for this archive file'
        write(lp,*) 'number of topo sea  mismatches = ',ibads
        write(lp,*) 'number of topo land mismatches = ',ibadl
        write(lp,*)
        call flush(lp)
        stop
      endif
      if     (ibadl.ne.0 .and. lhycom) then
        write(lp,*)
*       write(lp,*) 'warning - wrong bathymetry for this archive file'
        write(lp,*) 'error - wrong bathymetry for this archive file'
        write(lp,*) 'number of topo sea  mismatches = ',ibads
        write(lp,*) 'number of topo land mismatches = ',ibadl
        write(lp,*)
        call flush(lp)
        stop
      endif
c
cdiag call prtmsk(ip,depths,work,ii,ii1,jj1,0.,1.,'water depth (m)')
c
c --- cover islands with thin film of water for stream functions
        do 23 j=1,jj1
        do 23 i=1,ii1
 23     depth1(i,j)=depths(i,j)
c
        call sbmerg(depth1,film)
c
        call bigrd1(depth1)
c
cdiag   call prtmsk(ip,depth1,work,ii,ii1,jj1,0.,1.,'water depth (m)')
c
      call opngks
c
      ipalet=kpalet
      nbase =0
      ibase =0
      call colors	!  define color table
c
      call gsclip(0)
c
      do 3 k=1,kkin
      do 3 j=1,jj
      do 3 i=1,ii
c
c --- convert baroclinic to total velocities by adding barotropic component
c --- note that mean archives contain total velocity
      if (k.eq.1) then
        if     (iu(i,j).eq.1) then
          if     (artype.eq.1) then
            umix(i,j)=umix(i,j)+ubaro(i,j)  !ignore baclin, always total
          end if
        else !iu(i,j).ne.1
          umix(i,j)=0.
        end if
        if     (iv(i,j).eq.1) then
          if     (artype.eq.1) then
            vmix(i,j)=vmix(i,j)+vbaro(i,j)  !ignore baclin, always total
          end if
        else !iv(i,j).ne.1
          vmix(i,j)=0.
        end if
      endif
      if     (iu(i,j).eq.1) then
        if     (artype.eq.1 .and. .not.baclin) then
          u(i,j,2*k)=u(i,j,2*k)+ubaro(i,j)  !total velocity
        elseif (artype.eq.2 .and.      baclin) then
          u(i,j,2*k)=u(i,j,2*k)-ubaro(i,j)  !baroclinic velocity
        end if
      else !iu(i,j).ne.1
        u(i,j,2*k)=0.
      end if
      if     (iv(i,j).eq.1) then
        if     (artype.eq.1 .and. .not.baclin) then
          v(i,j,2*k)=v(i,j,2*k)+vbaro(i,j)  !total velocity
        elseif (artype.eq.2 .and.      baclin) then
          v(i,j,2*k)=v(i,j,2*k)-vbaro(i,j)  !baroclinic velocity
        end if
      else !iv(i,j).ne.1
        v(i,j,2*k)=0.
      end if
c
c --- convert layer thickness to meters
      if (depths(i,j).gt.0.) then
        dp(i,j,k)=dp(i,j,k)/9806.
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
        th3d(i,j,2*k)=th3d(i,j,2*k)+thbase
      else
        saln(i,j,2*k)=flag
        temp(i,j,2*k)=flag
        th3d(i,j,2*k)=flag
        ke(  i,j,2*k)=flag
        dp(i,j,k)=flag
        p(i,j,k+1)=flag
        if (depths(i,j).eq.film) p(i,j,k+1)=film
        do ktr= 1,ntracr
          trcr(i,j,2*k,ktr)=flag
        enddo
      endif
 3    continue
c
c --- fluid vertical velocity
      do j= 1,jj-1
        do i= 1,ii-1
          if     (ip(i,j).eq.1) then
            dudxdn=
     &            (u(i+1,j  ,2)*scuy(i+1,j  )-u(i,j,2)*scuy(i,j))
     &            /(scpx(i,j)*scpy(i,j))
            dvdydn=
     &            (v(i  ,j+1,2)*scvx(i  ,j+1)-v(i,j,2)*scvx(i,j))
     &        /(scpx(i,j)*scpy(i,j))
            w(i,j,2)=0.5*dp(i,j,1)*(dudxdn+dvdydn)
          else
            w(i,j,2) = flag
          endif
        enddo
      enddo
      do k= 2,kkin
        do j= 1,jj-1
          do i= 1,ii-1
            if     (iu(i,j).eq.1) then
              dpdx(i,j)=
     &                (p(i,j,k)*scpy(i,j)-p(i-1,j  ,k)*scpy(i-1,j  ))
     &                /(scux(i,j)*scuy(i,j))
            endif
          enddo
        enddo
        do j= 1,jj-1
          do i= 1,ii-1
            if     (iv(i,j).eq.1) then
              dpdy(i,j)=
     &                (p(i,j,k)*scpx(i,j)-p(i  ,j-1,k)*scpx(i  ,j-1))
     &                /(scvx(i,j)*scvy(i,j))
            endif
          enddo
        enddo
        do j=1,jj
          if     (iu(2   ,j).eq.1) then
            dpdx(1 ,j)=dpdx(2   ,j)
          endif
          if     (iu(ii-1,j).eq.1) then
            dpdx(ii,j)=dpdx(ii-1,j)
          endif
        enddo
        do i=1,ii
          if     (iv(i,2   ).eq.1) then
            dpdy(i,1 )=dpdy(i,2   )
          endif
          if     (iv(i,jj-1).eq.1) then
            dpdy(i,jj)=dpdy(i,jj-1)
          endif
        enddo
        do j= 1,jj-1
          do i= 1,ii-1
            if     (ip(i,j).eq.1) then
              dudxup=
     &              (u(i+1,j  ,2*k-2)*scuy(i+1,j  )-
     &               u(i  ,j  ,2*k-2)*scuy(i  ,j  ))
     &              /(scpx(i,j)*scpy(i,j))
              dvdyup=
     &              (v(i  ,j+1,2*k-2)*scvx(i  ,j+1)-
     &               v(i  ,j  ,2*k-2)*scvx(i  ,j  ))
     &              /(scpx(i,j)*scpy(i,j))
              dudxdn=
     &              (u(i+1,j  ,2*k  )*scuy(i+1,j  )-
     &               u(i  ,j  ,2*k  )*scuy(i  ,j  ))
     &              /(scpx(i,j)*scpy(i,j))
              dvdydn=
     &              (v(i  ,j+1,2*k  )*scvx(i  ,j+1)-
     &               v(i  ,j  ,2*k  )*scvx(i  ,j  ))
     &              /(scpx(i,j)*scpy(i,j))
              w(i,j,2*k)=w(i,j,2*k-2)+0.5*(dp(i,j,k-1)*(dudxup+dvdyup)+
     &                                     dp(i,j,k  )*(dudxdn+dvdydn)-
     &                  (u(i  ,j  ,2*k)-u(i  ,j  ,2*k-2))*dpdx(i  ,j  )-
     &                  (u(i+1,j  ,2*k)-u(i+1,j  ,2*k-2))*dpdx(i+1,j  )-
     &                  (v(i  ,j  ,2*k)-v(i  ,j  ,2*k-2))*dpdy(i  ,j  )-
     &                  (v(i  ,j+1,2*k)-v(i  ,j+1,2*k-2))*dpdy(i  ,j+1))
            else
              w(i,j,2*k) = flag
            endif
          enddo
        enddo
        do i= 1,ii
          w(i ,jj,2*k) = flag
        enddo
        do j= 1,jj
          w(ii,j ,2*k) = flag
        enddo
      enddo
c
c --- w is noisy - smooth at least once.
c
      if     (smooth) then
        do k= 1,kkin
          call psmoo(w(1,1,2*k),work)
          call psmoo(w(1,1,2*k),work)
        enddo
      else
        do k= 1,kkin
          call psmoo(w(1,1,2*k),work)
        enddo
      endif
c
ccc      x=thrufl(107,209,122,212,'(Drake Passage)')
ccc      x=thrufl(41,199,44,201,'(Florida Straits)')
ccc      x=thrufl(63,76,69,94,'(Indonesia)')
c
      do 7 j=1,jj
      do 7 i=1,ii
      if (depths(i,j).gt.0.) then
        srfht( i,j)=srfht( i,j)/(thref*98.06)  ! cm
        dpbl(  i,j)=dpbl(  i,j)/9806.          ! m
        dpmixl(i,j)=dpmixl(i,j)/9806.          ! m
        thmix( i,j)=thmix( i,j)+thbase         ! SigmaT
        if     (artype.ne.3) then
          ttrend(i,j)=surflx(i,j)*thref*8.64E4
     &                    /spcifh/dpbl(i,j)      ! deg/day
          strend(i,j)=salflx(i,j)*thref*8.64E4
     &                           /dpbl(i,j)      ! psu/day
          emnp(  i,j)=salflx(i,j)*thref*8.64E7
     &                           /saln(i,j,2)    ! mm/day
        else
          ttrend(i,j)=0.0
          strend(i,j)=0.0
          emnp(  i,j)=salflx(i,j)*thref*8.64E7
     &                           /35.0           ! mm/day (approx)
        endif
        if     (covice(i,j).eq.0.0) then
          thkice(i,j)= 0.0
          temice(i,j)=-1.8
        endif
      else
        srfht( i,j)=flag
        surflx(i,j)=flag
        salflx(i,j)=flag
        ttrend(i,j)=flag
        strend(i,j)=flag
          emnp(i,j)=flag
        covice(i,j)=flag
        thkice(i,j)=flag
        temice(i,j)=flag
        dpbl(  i,j)=flag
        dpmixl(i,j)=flag
        kebaro(i,j)=flag
        tmix( i,j)=flag
        smix( i,j)=flag
        thmix(i,j)=flag
        kemix(i,j)=flag
      end if
      if (kkin.eq.1 .or. kkin.lt.kk) then
        if (depths(i,j).gt.0.) then
          p(i,j,kk+1)=depths(i,j)
        else
          p(i,j,kk+1)=flag
        endif
      endif
 7    continue
      qq=contur(srfht,ii,ii1,jj1)
      write(6,*) 'srfht = ',amn,amx,qq
c
      dpth=0.5*onecm
c
c --- put vertically averaged t,s values into massless layers
c --- only if not plotted horizontally (improves vertical plots)
c
      if (mthin) then
      do 70 k=2,kkin
c
      do 70 j=1,jj1
      do 70 i=1,ii1
c
      if (depths(i,j).gt.0.) then
        temp(i,j,2*k-1)=0.
        saln(i,j,2*k-1)=0.
        th3d(i,j,2*k-1)=0.
        do ktr= 1,ntracr
          trcr(i,j,2*k-1,ktr)=0.
        enddo
        pmid=.5*(p(i,j,k)+p(i,j,k+1))
        phi=pmid+dpth
        plo=pmid-dpth
c
        sum=0.
        do k1=1,kkin
          delp=max(0.,min(p(i,j,k1+1),phi)-max(p(i,j,k1),plo))
          sum=sum+delp
          temp(i,j,2*k-1)=temp(i,j,2*k-1)+temp(i,j,2*k1)*delp
          saln(i,j,2*k-1)=saln(i,j,2*k-1)+saln(i,j,2*k1)*delp
          th3d(i,j,2*k-1)=th3d(i,j,2*k-1)+th3d(i,j,2*k1)*delp
          do ktr= 1,ntracr
            trcr(i,j,2*k-1,ktr)=trcr(i,j,2*k-1,ktr)+
     &                          trcr(i,j,2*k1, ktr)*delp
          enddo
        enddo
c
        if     (sum.gt.0.0) then
          temp(i,j,2*k)=temp(i,j,2*k-1)/sum
          saln(i,j,2*k)=saln(i,j,2*k-1)/sum
          th3d(i,j,2*k)=th3d(i,j,2*k-1)/sum
          do ktr= 1,ntracr
            trcr(i,j,2*k,ktr)=trcr(i,j,2*k-1,ktr)/sum
          enddo
        else
          write(lp,'(/ a,2i5,a /)')
     &      'error - subregion point i,j =',i,j,
     &      ' is sea in the bathymetry but land in the archive'
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
      end if
 70   continue
      endif !mthin
c
      if (smooth) then
c
c --- smooth mass field variables
c
      call psmoo(temp(1,1,2),work)
      call psmoo(saln(1,1,2),work)
      call psmoo(th3d(1,1,2),work)
      call psmoo(tmix,work)
      call psmoo(smix,work)
      call psmoo(thmix,work)
      do ktr= 1,ntracr
        call psmoo(trcr(1,1,2,ktr),work)
      enddo
c
      do 38 k=2,kkin
c
      do 76 j=1,jj1
      do 76 i=1,ii1
      if (depths(i,j).gt.0.) then
        util1(i,j)=max(onemm,dp(i,j,k))
        temp(i,j,2*k-1)=temp(i,j,2*k)*util1(i,j)
        saln(i,j,2*k-1)=saln(i,j,2*k)*util1(i,j)
        th3d(i,j,2*k-1)=th3d(i,j,2*k)*util1(i,j)
        do ktr= 1,ntracr
          trcr(i,j,2*k-1,ktr)=trcr(i,j,2*k,ktr)*util1(i,j)
        enddo
      else
        temp(i,j,2*k-1)=flag
        saln(i,j,2*k-1)=flag
        th3d(i,j,2*k-1)=flag
        do ktr= 1,ntracr
          trcr(i,j,2*k-1,ktr)=flag
        enddo
      end if
 76   continue
c
      call psmoo(util1,work)
      call psmoo(temp(1,1,2*k-1),work)
      call psmoo(saln(1,1,2*k-1),work)
      call psmoo(th3d(1,1,2*k-1),work)
      do ktr= 1,ntracr
        call psmoo(trcr(1,1,2*k-1,ktr),work)
      enddo
c
      do 38 j=1,jj1
      do 38 i=1,ii1
      if (depths(i,j).gt.0.) then
        temp(i,j,2*k)=temp(i,j,2*k-1)/util1(i,j)
        saln(i,j,2*k)=saln(i,j,2*k-1)/util1(i,j)
        th3d(i,j,2*k)=th3d(i,j,2*k-1)/util1(i,j)
        do ktr= 1,ntracr
          trcr(i,j,2*k,ktr)=trcr(i,j,2*k-1,ktr)/util1(i,j)
        enddo
      end if
 38   continue
c
c --- smooth velocity and layer thickness fields
c
      do 30 k=1,kkin
c
      do 31 j=1,jj1
      do 31 i=2,ii1
      if(k.eq.1) umix(i,j)=umix(i,j)*(dpmixl(i,j)+dpmixl(i-1,j))
 31   uflux(i,j)=u(i,j,2*k)*max(onecm,dp(i,j,k)+dp(i-1,j,k))
c
      do 32 j=2,jj1
      do 32 i=1,ii1
      if(k.eq.1) vmix(i,j)=vmix(i,j)*(dpmixl(i,j)+dpmixl(i,j-1))
 32   vflux(i,j)=v(i,j,2*k)*max(onecm,dp(i,j,k)+dp(i,j-1,k))
c
      if(k.eq.1) then
        call usmoo(umix,work)
        call vsmoo(vmix,work)
      end if
      call usmoo(uflux,work)
      call vsmoo(vflux,work)
      call psmoo(dp(1,1,k),work)
c --- (warning: smoothed -dp- field unsuitable for deriving interface depths)
c
      do 33 j=1,jj1
      do 33 i=2,ii1
      if(k.eq.1) umix(i,j)=umix(i,j)/(dpmixl(i,j)+dpmixl(i-1,j))
 33   u(i,j,2*k)=uflux(i,j)/max(onecm,dp(i,j,k)+dp(i-1,j,k))
c
      do 34 j=2,jj1
      do 34 i=1,ii1
      if(k.eq.1) vmix(i,j)=vmix(i,j)/(dpmixl(i,j)+dpmixl(i,j-1))
 34   v(i,j,2*k)=vflux(i,j)/max(onecm,dp(i,j,k)+dp(i,j-1,k))
c
c --- now smooth layer interfaces and find corresponding -dp- field
      if (k.lt.kkin) call psmo1(p(1,1,k+1),work,p(1,1,kk+1))
c --- now smooth boundary layer thickness and mixed layer base
      if (k.eq.1) then
        call psmo1(dpbl,work,p(1,1,kk+1))
        call psmo1(dpmixl,work,p(1,1,kk+1))
      end if
c
      do 35 j=1,jj1
      do 35 i=1,ii1
      if (depths(i,j).gt.0.) dp(i,j,k)=p(i,j,k+1)-p(i,j,k)
 35   continue
c
 30   continue
c
      end if			!  smooth = .true.
c
c --- put vertically averaged u,v values into massless layers
c --- only if not plotted horizontally (improves vertical plots)
c
      if (mthin) then
      do 74 k=2,kkin
c
      do 72 j=1,jj1
      do 72 i=2,ii1
      if (min(depths(i,j),depths(i-1,j)).gt.0.) then
c
        u(i,j,2*k-1)=0.
        pmid=.25*(p(i,j,k)+p(i-1,j,k)+p(i,j,k+1)+p(i-1,j,k+1))
        plo=pmid-dpth
        phi=pmid+dpth
c
        sum=0.
        do 73 k1=1,kkin
        delp=max(0.,min(.5*(p(i,j,k1+1)+p(i-1,j,k1+1)),phi)
     &             -max(.5*(p(i,j,k1  )+p(i-1,j,k1  )),plo))
        sum=sum+delp
 73     u(i,j,2*k-1)=u(i,j,2*k-1)+u(i,j,2*k1)*delp
c
        u(i,j,2*k)=u(i,j,2*k-1)/sum
      endif
 72   continue
c
      do 74 j=2,jj1
      do 74 i=1,ii1
      if (min(depths(i,j),depths(i,j-1)).gt.0.) then
c
        v(i,j,2*k-1)=0.
        pmid=.25*(p(i,j,k)+p(i,j-1,k)+p(i,j,k+1)+p(i,j-1,k+1))
        plo=pmid-dpth
        phi=pmid+dpth
c
        sum=0.
        do 75 k1=1,kkin
        delp=max(0.,min(.5*(p(i,j,k1+1)+p(i,j-1,k1+1)),phi)
     &             -max(.5*(p(i,j,k1  )+p(i,j-1,k1  )),plo))
        sum=sum+delp
 75     v(i,j,2*k-1)=v(i,j,2*k-1)+v(i,j,2*k1)*delp
c
        v(i,j,2*k)=v(i,j,2*k-1)/sum
      endif
 74   continue
      endif !mthin
c
      do 97 k=1,kkin
      do 97 j=1,jj1
      do 97 i=1,ii1
        ke(i,j,2*k-1)=  ke(i,j,2*k)
         w(i,j,2*k-1)=   w(i,j,2*k)
      temp(i,j,2*k-1)=temp(i,j,2*k)
      saln(i,j,2*k-1)=saln(i,j,2*k)
      th3d(i,j,2*k-1)=th3d(i,j,2*k)
      do ktr= 1,ntracr
        trcr(i,j,2*k-1,ktr)=trcr(i,j,2*k,ktr)
      enddo
 97   continue
c
      if     (artype.eq.1) then
        if     (yrflgi.lt.0) then  ! use model day
          write (label(51:72),'(a,f8.2,2x)') '  model day:',time
        else  ! use model date
          call fordate(time,yrflag, iyear,month,iday,ihour)
          if     (yrflag.lt.3) then
            write (label(51:72),112) time/year,cmonth(month),iday
          else
*           write (label(51:72),113) cmonth(month),iday,iyear
            write (label(51:72),123) cmonth(month),iday,iyear,ihour
          endif
        endif
      else  ! mean or sdev archive
        write(lp,*) 'time3 = ',time3
        dt = 0.5*(time3(2)-time3(1))/(nstep-1)
        if     (yrflag.eq.0) then
          dt0 = 15.0
        elseif (yrflag.eq.1) then
          dt0 = 15.25
        elseif (yrflag.eq.2) then
          dt0 = 0.0
        else
          dt0 = 0.0
        endif
        if     (artype.eq.2) then
          if     (yrflag.ne.3) then
            write(label(51:72),114)
     &        ' mean: ',(time3(1)-dt+dt0)/year,
     &                  (time3(2)+dt+dt0)/year
          else
            write(label(51:72),114)
     &        ' mean: ',(time3(1)-dt+dt0)/year+1901.0,
     &                  (time3(2)+dt+dt0)/year+1901.0
          endif
        else
          if     (yrflag.ne.3) then
            write(label(51:72),114) 
     &        ' sdev: ',(time3(1)-dt+dt0)/year,
     &                  (time3(2)+dt+dt0)/year
          else
            write(label(51:72),114) 
     &        ' sdev: ',(time3(1)-dt+dt0)/year+1901.0,
     &                  (time3(2)+dt+dt0)/year+1901.0
          endif
        endif
      endif
      if (lhycom) then
        write (label(73:81),115) iexpt/10,mod(iexpt,10),'H'
      else
        write (label(73:81),115) iexpt/10,mod(iexpt,10),'M'
      endif
 112  format ('  year',f7.2,' (',a3,i3.2,')')
*113  format ('  date: ',a3,i3.2,',',i5,'  ')
 123  format ('   ',a3,i3.2,',',i5,'  ',i2.2,'Z  ')
 114  format (a7,f7.2,'-',f7.2)
 115  format (' [',i2.2,'.',i1.1,a1,']')
      write (lp,'('' now plotting:'',a30)') label(51:81)
      call flush(lp)
c
c --- -----------------------
c --- plot non-layered fields
c --- -----------------------
c
      k=1
c
      nquad=-1				!  next plot full frame
c
c --- 'botqq ' = bathymetry contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'botqq ','("blkinr: ",a6," =",f11.4," m")')
      if (qqin.ge.0.0) then
      
      label(33:50)=' bathymetry       '
c
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=p(i,j,kk+1)
        enddo
      enddo
c
      qq=max(10.,contur(util1,ii,ii1,jj1))
      if (qqin.ne.0.0) qq=qqin
      ipalet=0
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,.2*qq,.2*qq,.2*qq,
     &         495,blank,.5,nquad,.false.,lalolb,lalogr)
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," m")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.1,'' m'')') qq
        call pcloqu(xlab1,ylab1,text(1:10),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
c
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- -------------------
c --- plot surface fluxes
c --- -------------------
c
c --- 'flxqq ' = surf. heat  flux contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'flxqq ','("blkinr: ",a6," =",f11.4," cm")')
      if (qqin.ge.0.0) then
      qq=min(50.,max( 2.,contur(surflx,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      write(6,*) 'surflx= ',amn,amx,qq
      label(33:50)=' surf. heat flux  '
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," cm")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(surflx,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.1,'' w/m^2'')') qq
        call pcloqu(xlab1,ylab1,text(1:14),csn,0.,0)
        write (text,'(f7.1,'' to'',f8.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:18),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c --- 'empqq ' = surf. evap-pcip contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'empqq ','("blkinr: ",a6," =",
     &                 f11.4," mm/day")')
      if (qqin.ge.0.0) then
      qq=min(500.,max( 2.,contur(  emnp,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      write(6,*) '  emnp= ',amn,amx,qq
      label(33:50)=' surf. evap-precip'
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",
     &                 f11.4," mm/day")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(  emnp,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.1,'' mm/day'')') qq
        call pcloqu(xlab1,ylab1,text(1:15),csn,0.,0)
        write (text,'(f7.1,'' to'',f8.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:18),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c --- 'ttrqq ' = surf. temp trend contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'ttrqq ','("blkinr: ",a6," =",f11.4," deg/day")')
      if (qqin.ge.0.0) then
      qq=min(50.,max( 2.,contur(ttrend,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      write(6,*) 'ttrend= ',amn,amx,qq
      label(33:50)=' surf. temp. trend'
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",
     &                              f11.4," deg/day")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(ttrend,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.3,'' deg/day'')') qq
        call pcloqu(xlab1,ylab1,text(1:16),csn,0.,0)
        write (text,'(f7.3,'' to'',f8.3)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:18),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c --- 'strqq ' = surf. saln trend contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'strqq ','("blkinr: ",a6," =",f11.4," psu/day")')
      if (qqin.ge.0.0) then
      qq=min(50.,max( 2.,contur(strend,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      write(6,*) 'strend= ',amn,amx,qq
      label(33:50)=' surf. saln. trend'
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",
     &                              f11.4," psu/day")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(strend,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.3,'' psu/day'')') qq
        call pcloqu(xlab1,ylab1,text(1:16),csn,0.,0)
        write (text,'(f7.3,'' to'',f8.3)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:18),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
      if (nquad.gt.0) then
        call fram(ncount)
        nquad=0
      end if
c
c --- ---------------
c --- plot ice fields
c --- ---------------
c
*     if     (icegln) then
c
c --- 'icvqq ' = ice coverage contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'icvqq ','("blkinr: ",a6," =",f11.4," ")')
      if (qqin.ge.0.0) then
        qq=min(1.0,max(0.01,contur(covice,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin
        write(6,*) 'covice = ',amn,amx,qq
        label(33:50)='     ice coverage '
        ipalet=kpalet
        if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---     'center' = central contoured value
          call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," ")')
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(covice,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.2)') qq
        call pcloqu(xlab1,ylab1,text(1:8), csn,0.,0)
        write (text,'(f5.3,'' to'',f6.3)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:14),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
      endif
c
c --- 'ithqq ' = ice thickness contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'ithqq ','("blkinr: ",a6," =",f11.4," m")')
      if (qqin.ge.0.0) then
        qq=min(10.0,max(0.01,contur(thkice,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin
        write(6,*) 'thkice = ',amn,amx,qq
        label(33:50)='    ice thickness '
        ipalet=kpalet
        if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---     'center' = central contoured value
          call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," m")')
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(thkice,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.2,'' m'')') qq
        call pcloqu(xlab1,ylab1,text(1:10),csn,0.,0)
        write (text,'(f5.2,'' to'',f6.2)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:14),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
      endif
c
c --- 'ictqq ' = ice temperature contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'ictqq ','("blkinr: ",a6," =",f11.4," deg")')
      if (qqin.ge.0.0) then
        qq=min(5.,max(.1,contur(temice,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin
        label(33:50)='   ice temperture '
        ipalet=kpalet
        if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---     'center' = central contoured value
          call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," deg")')
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(temice,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f6.3,'' deg'')') qq
        call pcloqu(xlab1,ylab1,text(1:13),csn,0.,0)
        write (text,'(f6.2,'' to'',f6.2)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:15),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
      endif
c
      if (nquad.gt.0) then
        call fram(ncount)
        nquad=0
      end if
c
*     endif  ! icegln
c
c --- -------------------
c --- plot surface fields
c --- -------------------
c
c --- 'sshqq ' = sea surf. height contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'sshqq ','("blkinr: ",a6," =",f11.4," cm")')
      if (qqin.ge.0.0) then
      qq=min(50.,max( 2.,contur(srfht,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      write(6,*) 'srfht = ',amn,amx,qq
      label(33:50)=' sea surf. height '
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," cm")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(srfht,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm'')') qq
        call pcloqu(xlab1,ylab1,text(1:10),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
      if     (artype.gt.1) then  ! mean or std. archive
c --- 'bkeqq ' = baro. kinetic energy contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'bkeqq ','("blkinr: ",a6," =",
     &                 f11.4," m^2/s^2")')
      if (qqin.ge.0.0) then
      qq=contur(kebaro,ii,ii1,jj1)
      if (qqin.ne.0.0) qq=qqin
      write(6,*) 'kebaro= ',amn,amx,qq
      label(33:50)='    baro.k.e./mass'
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",
     &                 f11.4," m^2/s^2")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(kebaro,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f6.4,'' m^2/s^2'')') qq
        call pcloqu(xlab1,ylab1,text(1:17),csn,0.,0)
        write (text,'(f6.3,'' to'',f7.3)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
      endif  ! mean or std. archive
c
c --- ------------------------
c --- plot barotropic velocity
c --- ------------------------
c --- 'bspdqq' = barotropic speed contour interval  (<0 no plot; 0 from field)
c --- 'center' = central contoured value (if needed)
c --- 'bthrsh' = barotropic velocity plot threshold (0 no plot; <0 contour int)
      call blkinr2(qq,i,
     &             'bspdqq','("blkinr: ",a6," =",f11.4," cm/s")',
     &             'bsfqq ','("blkinr: ",a6," =",f11.4," sv"  )')
      if     (i.eq.1) then !'bspdqq'
        spdqq  = qq
        if (kpalet.ge.2 .and. spdqq.gt.0.0) then
c ---     'center' = central contoured value
          call blkinr(qqcspd,
     &                'center','("blkinr: ",a6," =",f11.4," cm/s")')
        endif
c ---   also have 'vthrsh' and 'bsfqq '
        call blkinr(thresh,
     &              'bthrsh','("blkinr: ",a6," =",f11.4," cm/s")')
        call blkinr(qqstrm,
     &              'bsfqq ','("blkinr: ",a6," =",f11.4," sv")')
      else !'bsfqq '
        spdqq  = -1.0 !no barotropic speed plot
        thresh =  0.0 !no barotropic velocity plots
        qqstrm =  qq  !'bthrsh'
      endif !i.eq.1:else
c
      if     (spdqq.ge.0.0) then  ! signals separate speed contour plot
c
c ---   speed.
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              util1(i,j)=50.0*sqrt( (ubaro(i,j)+ubaro(i+1,j))**2 +
     &                              (vbaro(i,j)+vbaro(i,j+1))**2  )
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (spdqq.gt.0.0) qq=spdqq
        write(6,*) 'bspeed = ',amn,amx,qq
        label(33:50)='    baro. speed   '
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc=qqcspd
          qqmn = qqc-0.5*qq*cntrs(kpalet)
          qqmx = qqc+0.5*qq*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
      endif  !spdqq.ge.0.0
c       
      if     (thresh.lt.0.0) then  ! signals u-v-speed contour plots
        qqin = -thresh
c       
c ---   u-velocity
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii) then
              util1(i,j)=50.0*(ubaro(i,j)+ubaro(i+1,j))
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'u-bvel = ',amn,amx,qq
        label(33:50)='    baro. u-vel.  '
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc = 0.0
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
c       
c ---   v-velocity
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. j.lt.jj) then
              util1(i,j)=50.0*(vbaro(i,j)+vbaro(i,j+1))
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'v-bvel = ',amn,amx,qq
        label(33:50)='    baro. v-vel.  '
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc = 0.0
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
c       
c ---   speed.
        qqin=0.5*qqin  ! because speed is positive
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              util1(i,j)=50.0*sqrt( (ubaro(i,j)+ubaro(i+1,j))**2 +
     &                              (vbaro(i,j)+vbaro(i,j+1))**2  )
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'speed  = ',amn,amx,qq
        label(33:50)='    baro. speed   '
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqmn = 0.0
          qqmx = qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
      endif  !thresh.lt.0.0
c       
      if (thresh.gt.0.0) then
      label(33:50)='    baro. velocity'
c --- call horplt to initiate next frame but don't draw any isolines
      do j=1,jwrk
        do i=1,iwrk
          work(i,j)=0.
        enddo
      enddo
      ipalet=0
      call horplt(work,plon,plat,ii,jj,ii1,jj1,1.,1.,2.,
     &            0,label(33:81),0.,nquad,.true.,lalolb,lalogr)
ccc     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
c    
c --- velocities in the range (0...thresh) are indicated by length of vector.
c --- velocities > thresh are indicated by barbs (1  per 'thresh' increment).
c --- draw only every i_th vector in every (i_th/2) row
      do i=2,ii1,max(1,i_th/2)
        do j=2+mod(i,i_th),jj1,i_th
          qu=100.0*(ubaro(i,j)+ubaro(i+1,j))*float(i_th)/(4.*thresh)
          qv=100.0*(vbaro(i,j)+vbaro(i,j+1))*float(i_th)/(4.*thresh)
          call arrow1(float(i)-qu,float(j)-qv
     &               ,float(i)+qu,float(j)+qv,float(i_th))
        enddo !j
      enddo !i
      call legend1(xlab1,ylab1,thresh,float(i_th))
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- 'bsfqq ' = baro. strmfn. contour int (<0 no plot; 0 from field)
******call blkinr(qqin,'bsfqq ','("blkinr: ",a6," =",f11.4," sv")')
      qqin = qqstrm
      if (qqin.ge.0.0) then
      label(33:50)='barotropic strmf. '
c
      do j= 1,jj
        jm1=max(1,j-1)
        do i= 1,ii
          im1=max(1,i-1)
c
          util1(i,j)=1.
c
          if (iu(i,j).eq.1) then
            ubaro(i,j)=ubaro(i,j)*min(p(i,j,kk+1),p(im1,j,kk+1))
          else
            ubaro(i,j)=0.0
          endif
c
          if (iv(i,j).eq.1) then
            vbaro(i,j)=vbaro(i,j)*min(p(i,j,kk+1),p(i,jm1,kk+1))
          else
            vbaro(i,j)=0.0
          endif
c
          if (iq(i,j).eq.1) then
            vort(i,j)=
     &       (vbaro(i,j)*scvy(i,j)-vbaro(im1,j)*scvy(im1,j)
     &       -ubaro(i,j)*scux(i,j)+ubaro(i,jm1)*scux(i  ,jm1))*1.e-6
          else
            vort(i,j)=0.0
          endif
        enddo 
      enddo 
c
c --- calculate first estimate for strmf values
      ubi=0.0
      do j=1,jj-1
        ubim1=ubi
        ubi  =ubaro(1,j)*scux(1,j)
        if     (j.eq.1) then
          strmfu=0.0
        else
          strmfu=strmfu-0.5*(ubi+ubim1)
        endif
        strmft=strmfu
        if     (iq(1,j).eq.1) then
          strmf(1,j)=strmft*1.e-6
        else
          strmf(1,j)=0.
        endif
        vbi=0.0
        do i=2,ii
          vbim1 =vbi
          vbi   =vbaro(i,j)*scvy(i,j)
          strmft=strmft+0.5*(vbi+vbim1)
          if     (iq(i,j).eq.1) then
            strmf(i,j)=strmft*1.e-6
          else
            strmf(i,j)=0.
          endif
        enddo
      enddo
      do i=1,ii
        strmf(i,jj)=0.
      enddo
c
      call poisnd(ii,jj,vort,strmf,util1,util1,work)
      call qsmoo(strmf,work)
c
c --- interpolate from q-grid to p-grid.
      do j=1,jj
        do i=1,ii
          if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
            s1=0.0
            s2=0.0
            if     (iq(i,  j)  .eq.1) then
              s1 = s1 + 1.0
              s2 = s2 + strmf(i,  j)
            endif
            if     (iq(i+1,j)  .eq.1) then
              s1 = s1 + 1.0
              s2 = s2 + strmf(i+1,j)
            endif
            if     (iq(i,  j+1).eq.1) then
              s1 = s1 + 1.0
              s2 = s2 + strmf(i,  j+1)
            endif
            if     (iq(i+1,j+1).eq.1) then
              s1 = s1 + 1.0
              s2 = s2 + strmf(i+1,j+1)
            endif
            if     (s1.ne.0.0) then
              util1(i,j)=s2/s1
            else
              util1(i,j)=flag
            endif
          else
            util1(i,j)=flag
          endif
        enddo
      enddo
c
ccc   call zebra(util1,ii,ii,jj)
c
      qq=min(1.,max(50.,contur(util1,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," sv")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = -15.5*qq
        qqmx =  15.5*qq
      endif
      call horplt(util1,plon,plat,ii,jj,ii,jj,qqmn,qqmx,qq,
     &            -495,label(33:81),0.,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' sv'')') qq
        call pcloqu(xlab1,ylab1,text(1:10),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
      if (nquad.gt.0) then
        call fram(ncount)
        nquad=0
      end if
c
c --- -----------------------
c --- plot mixed layer fields
c --- -----------------------
c
      k=1
c
      nquad=-1				!  next plot full frame
c
c --- 'mthrsh' = mix lay velocity plot threshold (0 no plot; <0 contour int)
      call blkinr(thresh,'mthrsh','("blkinr: ",a6," =",f11.4," cm")')
c
      if     (thresh.lt.0.0) then  ! signals u-v-speed contour plots
        qqin = -thresh
c
c ---   mixed layer u-velocity
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii) then
              util1(i,j)=50.0*(umix(i,j)+umix(i+1,j))
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'u-mix  = ',amn,amx,qq
        label(33:50)='mix.l. u-velocity '
        ipalet=kpalet
        if (kpalet.ge.2 .and. qqin.gt.0.0) then
          qqc = 0.0
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
c
c ---   mixed layer v-velocity
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. j.lt.jj) then
              util1(i,j)=50.0*(vmix(i,j)+vmix(i,j+1))
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'v-mix  = ',amn,amx,qq
        label(33:50)='mix.l. v-velocity '
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc = 0.0
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
c
c ---   mixed layer speed.
        qqin=0.5*qqin  ! because speed is positive
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              util1(i,j)=50.0*sqrt( (umix(i,j)+umix(i+1,j))**2 +
     &                              (vmix(i,j)+vmix(i,j+1))**2  )
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'speed  = ',amn,amx,qq
        label(33:50)='mixed-layer speed '
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqmn = 0.0
          qqmx = qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
      endif  !thresh.lt.0.0
c
      if (thresh.gt.0.0) then
      label(33:50)='mix.lyr. velocity '
c --- call horplt to initiate next frame but don't draw any isolines
      do 66 j=1,jwrk
      do 66 i=1,iwrk
 66    work(i,j)=0.
      ipalet=0
      call horplt(work,plon,plat,ii,jj,ii1,jj1,1.,1.,2.,
     &            0,label(33:81),0.,nquad,.true.,lalolb,lalogr)
ccc     .            0,label(33:81),.5,nquad,.true.)
c
c --- velocities in the range (0...thresh) are indicated by length of vector.
c --- velocities > thresh are indicated by barbs (1  per 'thresh' increment).
c --- draw only every i_th vector in every (i_th/2) row
      do 648 i=2,ii1,max(1,i_th/2)
      do 648 j=2+mod(i,i_th),jj1,i_th
      qu=100.0*(umix(i,j)+umix(i+1,j))*float(i_th)/(4.*thresh)
      qv=100.0*(vmix(i,j)+vmix(i,j+1))*float(i_th)/(4.*thresh)
      if (min(depths(i,j  ),depths(i-1,j  ),
     &        depths(i,j-1),depths(i-1,j-1),qu*qu+qv*qv).gt.film) then
        if     (mdens) then
          if (p(i,j,k)+onecm.le.p(i,j,k+1) .and.
     &        abs(th3d(i,j,2*k)-theta(k)).le.0.002) then
            call arrow1(float(i)-qu,float(j)-qv
     &                 ,float(i)+qu,float(j)+qv,float(i_th))
          endif
        elseif (mthin) then
          if (p(i,j,k)+onecm.le.p(i,j,k+1)) then
            call arrow1(float(i)-qu,float(j)-qv
     &                 ,float(i)+qu,float(j)+qv,float(i_th))
          endif
        else
          call arrow1(float(i)-qu,float(j)-qv
     &               ,float(i)+qu,float(j)+qv,float(i_th))
        endif
      endif
 648  continue
      call legend1(xlab1,ylab1,thresh,float(i_th))
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- 'bltqq ' = bnd. lay. thick. contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'bltqq ','("blkinr: ",a6," =",f11.4," m")')
      if (qqin.ge.0.0) then
      qq=max(10.,contur(dpbl,ii,ii1,jj1))
      if (qqin.ne.0.0) qq=qqin
      label(33:50)='bnd.layr.thickness'
      ipalet=0
      call horplt(dpbl,plon,plat,ii,jj,ii1,jj1,.2*qq,.2*qq,.2*qq,
     &         495,blank,.5,nquad,.false.,lalolb,lalogr)
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," m")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(dpbl,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.1,'' m'')') qq
        call pcloqu(xlab1,ylab1,text(1:10),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- 'mltqq ' = mix. lay. thick. contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'mltqq ','("blkinr: ",a6," =",f11.4," m")')
      if (qqin.ge.0.0) then
      qq=max(10.,contur(dpmixl,ii,ii1,jj1))
      if (qqin.ne.0.0) qq=qqin
      label(33:50)='mix.layr.thickness'
      ipalet=0
      call horplt(dpmixl,plon,plat,ii,jj,ii1,jj1,.2*qq,.2*qq,.2*qq,
     &         495,blank,.5,nquad,.false.,lalolb,lalogr)
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," m")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(dpmixl,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.1,'' m'')') qq
        call pcloqu(xlab1,ylab1,text(1:10),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- 'sstqq ' = mix. lay. temp.  contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'sstqq ','("blkinr: ",a6," =",f11.4," deg")')
      if (qqin.ge.0.0) then
      qq=min(5.,max(.1,contur(tmix,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      label(33:50)='mix.layr.temp     '
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," deg")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(tmix,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f6.3,'' deg'')') qq
        call pcloqu(xlab1,ylab1,text(1:13),csn,0.,0)
        write (text,'(f5.2,'' to'',f6.2)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:14),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- 'sssqq ' = mix. lay. saln.  contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'sssqq ','("blkinr: ",a6," =",f11.4," psu")')
      if (qqin.ge.0.0) then
      qq=min(2.,max(.02,contur(smix,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      label(33:50)='mix.layr.saln     '
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," psu")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(smix,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.3,'' psu'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.2,'' to'',f7.2)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- 'ssdqq ' = mix. lay. dens.  contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'ssdqq ','("blkinr: ",a6," =",f11.4," sig")')
      if (qqin.ge.0.0) then
      qq=min(.5,max(.1,contur(thmix,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      label(33:50)='mix.layr.dens     '
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," sig")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(thmix,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.3,'' sig'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f5.2,'' to'',f6.2)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:14),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
      if     (artype.gt.1) then  ! mean or std. archive
c --- 'mkeqq ' = m.l. kinetic energy contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'mkeqq ','("blkinr: ",a6," =",
     &                 f11.4," m^2/s^2")')
      if (qqin.ge.0.0) then
      qq=contur(kemix,ii,ii1,jj1)
      if (qqin.ne.0.0) qq=qqin
      write(6,*) 'kemix = ',amn,amx,qq
      label(33:50)='    mixl.k.e./mass'
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",
     &                 f11.4," m^2/s^2")')
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(kemix,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f6.4,'' m^2/s^2'')') qq
        call pcloqu(xlab1,ylab1,text(1:17),csn,0.,0)
        write (text,'(f6.3,'' to'',f7.3)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
      endif  ! mean or std. archive
c
      if (nquad.gt.0) then
        call fram(ncount)
        nquad=0
      end if
c
c --- --------------------
c --- plot selected layers
c --- --------------------
c
      do !layer loop
c --- 'kf    ' = first plot layer (=0 end layer plots; <0 label with layer #)
c --- 'kl    ' = last  plot layer
      call blkini(kin, 'kf    ')
      kf = abs(kin)
      if     (kf.eq.0) then
        exit
      endif
      call blkini(kl, 'kl    ')
      kl = abs(kl)
      if     (kl.gt.kkin) then
        write(lp,'(a)') 'error - kl larger than kdm'
        exit
      elseif (kl.lt.kf) then
        write(lp,'(a)') 'error - kl smaller than kf'
        exit
      endif
c
      do k= kf,kl  !k-loop
c
c --- -------------------
c --- plot layer velocity
c --- -------------------
c --- 'spdqq ' = layer k speed contour int (<0 no plot; 0 from field)
c --- 'vthrsh' = layer k velocity plot threshold (0 no plot; <0 contour int)
      if    (k.eq.kf) then
        call blkinr2(qq,i,
     &               'spdqq ','("blkinr: ",a6," =",f11.4," cm/s")',
     &               'vthrsh','("blkinr: ",a6," =",f11.4," cm/s")')
        if     (i.eq.1) then !''spdqq '
          spdqq  = qq
          if (kpalet.ge.2 .and. spdqq.gt.0.0) then
c ---       'center' = central contoured value
            call blkinr(qqcspd,
     &                  'center','("blkinr: ",a6," =",f11.4," cm/s")')
          endif
c ---     also have 'vthrsh'
          call blkinr(thresh,
     &                'vthrsh','("blkinr: ",a6," =",f11.4," cm/s")')
        else !'vthrsh'
          spdqq  = -1.0 !no separate speed plot
          thresh = qq
        endif !i.eq.1:else
      endif !k.eq.kf
c
      if     (spdqq.ge.0.0) then  ! signals separate speed contour plot
c
c ---   speed.
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              util1(i,j)=50.0*sqrt( (u(i,j,2*k)+u(i+1,j,2*k))**2 +
     &                              (v(i,j,2*k)+v(i,j+1,2*k))**2  )
              if (mthin .and. p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                util1(i,j)=flag
              elseif (mdens .and.
     &                abs(th3d(i,j,2*k)-theta(k)).gt.0.002) then
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (spdqq.gt.0.0) qq=spdqq
        write(6,*) 'speed  = ',amn,amx,qq
        if     (baclin) then
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' bcl.spd.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  bcl.spd.'')') k
          endif
        else !total velocity
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' speed   '')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  speed   '')') k
          endif
        endif
        label(33:50)=text
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc=qqcspd
          qqmn = qqc-0.5*qq*cntrs(kpalet)
          qqmx = qqc+0.5*qq*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
      endif  !spdqq.ge.0.0
c
      if     (thresh.lt.0.0) then  ! signals u-v-speed contour plots
        qqin = -thresh
c
c ---   u-velocity
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii) then
              util1(i,j)=50.0*(u(i,j,2*k)+u(i+1,j,2*k))
              if (mthin .and. p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                util1(i,j)=flag
              elseif (mdens .and.
     &                abs(th3d(i,j,2*k)-theta(k)).gt.0.002) then
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'u-vel  = ',amn,amx,qq
        if     (baclin) then
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' u-b.vel.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  u-b.vel.'')') k
          endif
        else !total velocity
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' u-veloc.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  u-veloc.'')') k
          endif
        endif
        label(33:50)=text
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc = 0.0
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
c
c ---   v-velocity
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. j.lt.jj) then
              util1(i,j)=50.0*(v(i,j,2*k)+v(i,j+1,2*k))
              if (mthin .and. p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                util1(i,j)=flag
              elseif (mdens .and.
     &                abs(th3d(i,j,2*k)-theta(k)).gt.0.002) then
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'v-vel  = ',amn,amx,qq
        if     (baclin) then
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' v-b.vel.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  v-b.vel.'')') k
          endif
        else !total velocity
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' v-veloc.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  v-veloc.'')') k
          endif
        endif
        label(33:50)=text
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc = 0.0
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
c
c ---   speed.
        qqin=0.5*qqin  ! because speed is positive
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              util1(i,j)=50.0*sqrt( (u(i,j,2*k)+u(i+1,j,2*k))**2 +
     &                              (v(i,j,2*k)+v(i,j+1,2*k))**2  )
              if (mthin .and. p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                util1(i,j)=flag
              elseif (mdens .and.
     &                abs(th3d(i,j,2*k)-theta(k)).gt.0.002) then
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'speed  = ',amn,amx,qq
        if     (baclin) then
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' bcl.spd.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  bcl.spd.'')') k
          endif
        else !total velocity
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' speed   '')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  speed   '')') k
          endif
        endif
        label(33:50)=text
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqmn = 0.0
          qqmx = qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' cm/s'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
      endif  !thresh.lt.0.0
c
      if (thresh.gt.0.0) then
      if     (baclin) then
        if (kin.gt.0) then
          write (text,'(''sig='',f5.2,   '' bcl.vel.'')') theta(k)
        else
          write (text,'(''layer='',i2.2,''  bcl.vel.'')') k
        endif
      else !total velocity
        if (kin.gt.0) then
          write (text,'(''sig='',f5.2,   '' velocity'')') theta(k)
        else
          write (text,'(''layer='',i2.2,''  velocity'')') k
        endif
      endif
      label(33:50)=text
c --- call horplt to initiate next frame but don't draw any isolines
      do 6 j=1,jwrk
      do 6 i=1,iwrk
 6    work(i,j)=0.
      ipalet=0
      call horplt(work,plon,plat,ii,jj,ii1,jj1,1.,1.,2.,
     &            0,label(33:81),0.,nquad,.true.,lalolb,lalogr)
ccc     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
c
c --- velocities in the range (0...thresh) are indicated by length of vector.
c --- velocities > thresh are indicated by barbs (1  per 'thresh' increment).
c --- draw only every i_th vector in every (i_th/2) row
      do 48 i=2,ii1,max(1,i_th/2)
      do 48 j=2+mod(i,i_th),jj1,i_th
      qu=100.0*(u(i,j,2*k)+u(i+1,j,2*k))*float(i_th)/(4.*thresh)
      qv=100.0*(v(i,j,2*k)+v(i,j+1,2*k))*float(i_th)/(4.*thresh)
      if (min(depths(i,j  ),depths(i-1,j  ),
     &        depths(i,j-1),depths(i-1,j-1),qu*qu+qv*qv).gt.film) then
        if (mthin) then
          if (p(i,j,k)+onecm.le.p(i,j,k+1)) then
            call arrow1(float(i)-qu,float(j)-qv
     &                 ,float(i)+qu,float(j)+qv,float(i_th))
          endif
        else
          call arrow1(float(i)-qu,float(j)-qv
     &               ,float(i)+qu,float(j)+qv,float(i_th))
        endif
      endif
  48  continue
      call legend1(xlab1,ylab1,thresh,float(i_th))
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- ----------------------------
c --- plot fluid vertical velocity
c --- ----------------------------
      if     (k.eq.kf) then
c --- 'wvelqq' = lay.  k w-vel. contour int (<0 no plot; 0 from field)
c --- 'infqq ' = lay.  k thick. contour int (<0 no plot; 0 from field)
      call blkinr2(qq,i,'wvelqq','("blkinr: ",a6," =",f11.4," cm/day")',
     &                  'infqq ','("blkinr: ",a6," =",f11.4," m")'    )
        if     (i.eq.1) then
          qqwvel = max( qq, -1.0 )
        else
          qqwvel = -99.0  ! signal that infqq has been input
          qqinf = qq
        endif
      endif
      qqin=qqwvel
      if (qqin.ge.0.0) then
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   '' w.veloc '')') theta(k)
      else
        write (text,'(''layer='',i2.2,''  w.veloc '')') k
      endif
      label(33:50)=text
c
      do j=1,jj
        do i=1,ii
          if (w(i,j,2*k).ne.flag) then
            util1(i,j)=8640000.0*w(i,j,2*k)
          else
            util1(i,j)=flag
          endif
        enddo
      enddo
c
      qq=max(10.,contur(util1,ii,ii1,jj1))
      if (qqin.ne.0.0) qq=qqin
      ipalet=0
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,.2*qq,.2*qq,.2*qq,
     &         495,blank,.5,nquad,.false.,lalolb,lalogr)
   
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqciwf,'center',
     &                     '("blkinr: ",a6," =",f11.4," cm/day")')
        endif
        qqc=qqciwf
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.2,'' cm/day'')') qq
        call pcloqu(xlab1,ylab1,text(1:15),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
c
      if (nquad.eq.0) call fram(ncount)
      endif !qqin.ge.0.0
c
c --- --------------------
c --- plot interface depth
c --- --------------------
c
c --- 'infqq ' = lay.  k thick. contour int (<0 no plot; 0 from field)
      if     (k.eq.kf .and. qqwvel.ne.-99.0) then
      call blkinr(qqinf,'infqq ','("blkinr: ",a6," =",f11.4," m")')
      endif
      qqin=qqinf
      if (qqin.ge.0.0) then
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   ''  i.depth'')') theta(k)
      else
        write (text,'(''layer='',i2.2,''   i.depth'')') k
      endif
      label(33:50)=text
c
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=p(i,j,k+1)
        enddo
      enddo
c
      qq=max(10.,contur(util1,ii,ii1,jj1))
      if (qqin.ne.0.0) qq=qqin
      ipalet=0
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,.2*qq,.2*qq,.2*qq,
     &         495,blank,.5,nquad,.false.,lalolb,lalogr)
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqcinf,'center','("blkinr: ",a6," =",f11.4," m")')
        endif
        qqc=qqcinf
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.1,'' m'')') qq
        call pcloqu(xlab1,ylab1,text(1:10),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
c
      if (nquad.eq.0) call fram(ncount)
      endif !qqin.ge.0.0
c
c --- --------------------
c --- plot layer thickness
c --- --------------------
c
c --- 'thkqq ' = lay.  k thick. contour int (<0 no plot; 0 from field)
      if     (k.eq.kf) then
      call blkinr(qqthk,'thkqq ','("blkinr: ",a6," =",f11.4," m")')
      endif
      qqin=qqthk
      if (qqin.ge.0.0) then
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   ''  thknss '')') theta(k)
      else
        write (text,'(''layer='',i2.2,''   thknss '')') k
      endif
      label(33:50)=text
c
      do 39 j=1,jj1
      do 39 i=1,ii1
      util1(i,j)=dp(i,j,k)
      if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
        util1(i,j)=flag
      endif
 39   continue
c
      qq=max(10.,contur(util1,ii,ii1,jj1))
      if (qqin.ne.0.0) qq=qqin
      ipalet=0
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,.2*qq,.2*qq,.2*qq,
     &         495,blank,.5,nquad,.false.,lalolb,lalogr)
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqcthk,'center','("blkinr: ",a6," =",f11.4," m")')
        endif
        qqc=qqcthk
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.1,'' m'')') qq
        call pcloqu(xlab1,ylab1,text(1:10),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
c
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- ----------------
c --- plot temperature
c --- ----------------
c
c --- 'temqq ' = layer k temp contour int (<0 no plot; 0 from field)
      if     (k.eq.kf) then
      call blkinr(qqtem,'temqq ','("blkinr: ",a6," =",f11.4," deg")')
      endif
      qqin=qqtem
      if (qqin.ge.0.0) then
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   ''  temp   '')') theta(k)
      else
        write (text,'(''layer='',i2.2,''   temp   '')') k
      endif
      label(33:50)=text
c
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=temp(i,j,2*k)
          if     (mdens) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1) .or.
     &          abs(th3d(i,j,2*k)-theta(k)).gt.0.002) then
              util1(i,j)=flag
            endif
          elseif (mthin) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
              util1(i,j)=flag
            endif
          elseif (artype.ne.3) then
            if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
              util1(i,j)=flag
            endif
          endif
        enddo
      enddo
      qq=min(.5,max(.1,contur(util1,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqctem,'center','("blkinr: ",a6," =",f11.4," deg")')
        endif
        qqc=qqctem
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f6.3,'' deg'')') qq
        call pcloqu(xlab1,ylab1,text(1:13),csn,0.,0)
        write (text,'(f5.2,'' to'',f6.2)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:14),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- -------------
c --- plot salinity
c --- -------------
c
c --- 'salqq ' = lay.  k saln. contour int (<0 no plot; 0 from field)
      if     (k.eq.kf) then
      call blkinr(qqsal,'salqq ','("blkinr: ",a6," =",f11.4," psu")')
      endif
      qqin=qqsal
      if (qqin.ge.0.0) then
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   '' salinity'')') theta(k)
      else
        write (text,'(''layer='',i2.2,''  salinity'')') k
      endif
      label(33:50)=text
c
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=saln(i,j,2*k)
          if     (mdens) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1) .or.
     &          abs(th3d(i,j,2*k)-theta(k)).gt.0.002) then
              util1(i,j)=flag
            endif
          elseif (mthin) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
              util1(i,j)=flag
            endif
          elseif (artype.ne.3) then
            if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
              util1(i,j)=flag
            endif
          endif
        enddo
      enddo
c
      qq=qqin
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqcsal,'center','("blkinr: ",a6," =",f11.4," psu")')
        endif
        qqc=qqcsal
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.3,'' psu'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        qq=contur(util1,ii,ii1,jj1)
        write (text,'(f5.2,'' to'',f6.2)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:14),csn,0.,0)
        if (nquad.eq.0) call fram(ncount)
      endif
c
c --- -------------
c --- plot density
c --- -------------
c
c --- 'tthqq ' = layer k density contour int (<0 no plot; 0 from field)
      if     (k.eq.kf) then
      call blkinr(qqtth,'tthqq ','("blkinr: ",a6," =",f11.4," sig")')
      endif
      qqin=qqtth
      if (qqin.ge.0.0) then
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=th3d(i,j,2*k)
          if     (mdens) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1) .or.
     &          abs(th3d(i,j,2*k)-theta(k)).gt.0.002) then
              util1(i,j)=flag
            endif
          elseif (mthin) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
              util1(i,j)=flag
            endif
          elseif (artype.ne.3) then
            if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
              util1(i,j)=flag
            endif
          endif
        enddo
      enddo
c
      qq=min(.5,max(.01,contur(util1,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   '' density '')') theta(k)
      else
        write (text,'(''layer='',i2.2,''  density '')') k
      endif
      label(33:50)=text
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqctth,'center',
     &                     '("blkinr: ",a6," =",f11.4,"kg/m^3")')
        endif
        qqc=qqctth
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f5.3,'' sig'')') qq
        call pcloqu(xlab1,ylab1,text(1:12),csn,0.,0)
        if     (kin.gt.0) then
          write (text,'(f6.3,'' to'',f6.3)') amn,amx
        else
          write (text,'(f6.2,'' to'',f6.2)') amn,amx
        endif
        call pcloqu(xlab2,ylab2,text(1:15),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- ------------
c --- plot tracers
c --- ------------
c
      do ktr= 1,ntracr
c --- 'trcqq ' = layer k tracer contour int (<0 no plot; 0 from field)
      if     (k.eq.kf) then
      write(cline,'(a,a,a)')
     &  '("blkinr: ",a6," =",f11.4,"', trim(ctrc_units(ktr)), '")'
      call blkinr(qqtrc(ktr),'trcqq ',cline)
      endif
      qqin=qqtrc(ktr)
      if (qqin.ge.0.0) then
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,  1x,a8)') theta(k),ctrc_title(ktr)
      else
        write (text,'(''layer='',i2.2,2x,a8)')       k, ctrc_title(ktr)
      endif
      label(33:50)=text
c
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=trcr(i,j,2*k,ktr)
          if     (mdens) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1) .or.
     &          abs(th3d(i,j,2*k)-theta(k)).gt.0.002) then
              util1(i,j)=flag
            endif
          elseif (mthin) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
              util1(i,j)=flag
            endif
          elseif (artype.ne.3) then
            if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
              util1(i,j)=flag
            endif
          endif
        enddo
      enddo
      qq=min(.5,max(.1,contur(util1,ii,ii1,jj1)))
      if (qqin.ne.0.0) qq=qqin
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        write(cline,'(a,a,a)')
     &    '("blkinr: ",a6," =",f11.4,"', trim(ctrc_units(ktr)), '")'
        call blkinr(qqctrc(ktr),'center',cline)
        endif
        qqc=qqctrc(ktr)
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write(cline,'(5a)')
     &    '("ci ",',ctrc_format(ktr),',1x,"',trim(ctrc_units(ktr)),'")'
        write(lp,'(a)') trim(cline)
        write (text,cline) qq
        call pcloqu(xlab1,ylab1,trim(text),csn,0.,0)
        write(cline,'(5a)')
     &    '(',ctrc_format(ktr),'," to",',ctrc_format(ktr),')'
        write(lp,'(a)') trim(cline)
        write (text,cline) amn,amx
        call pcloqu(xlab2,ylab2,trim(text),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
      enddo !ktr
c
c --- --------------------------
c --- plot layer kinetic energy
c --- --------------------------
c
      if     (artype.gt.1) then  ! mean or std. archive
c --- 'keqq  ' = kinetic energy contour int (<0 no plot; 0 from field)
      if     (k.eq.kf) then
      call blkinr(qqke,'keqq  ','("blkinr: ",a6," =",
     &                 f11.4," m^2/s^2")')
      endif
      qqin=qqke
      if (qqin.ge.0.0) then
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=ke(i,j,2*k)
          if     (mdens) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1) .or.
     &          abs(th3d(i,j,2*k)-theta(k)).gt.0.002) then
              util1(i,j)=flag
            endif
          elseif (mthin) then
            if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
              util1(i,j)=flag
            endif
          elseif (artype.ne.3) then
            if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
              util1(i,j)=flag
            endif
          endif
        enddo
      enddo
c
      qq=contur(util1,ii,ii1,jj1)
      if (qqin.ne.0.0) qq=qqin
      write(6,*) 'kemix = ',amn,amx,qq
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   '' ke/mass '')') theta(k)
      else
        write (text,'(''layer='',i2.2,''  ke/mass '')') k
      endif
      label(33:50)=text
c
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqcke,'center','("blkinr: ",a6," =",
     &                 f11.4," m^2/s^2")')
        endif
        qqc=qqcke
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),.5,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f6.4,'' m^2/s^2'')') qq
        call pcloqu(xlab1,ylab1,text(1:17),csn,0.,0)
        write (text,'(f6.3,'' to'',f7.3)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
      endif  ! mean or std. archive
c
c --- --------------------------
c --- plot layer stream function
c --- --------------------------
c
c --- 'sfnqq ' = layer k strmfn. contour int (<0 no plot; 0 from field)
      if     (k.eq.kf) then
      call blkinr(qqsfn,'sfnqq ','("blkinr: ",a6," =",f11.4," sv")')
      endif
      qqin=qqsfn
      if (qqin.ge.0.0) then
      write (text,'(''layer 1     strmf '')')
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   ''   strmf '')') theta(k)
      else
        write (text,'(''layer='',i2.2,''    strmf '')') k
      endif
      label(33:50)=text
c
      do j=1,jj
        jm1=max(1,j-1)
        do i= 1,ii
          im1=max(1,i-1)
c
          util1(i,j)=1.
c
          if (iu(i,j).eq.1) then
            ubaro(i,j)=u(i,j,2*k)*.5*(dp(i,j,k)+dp(im1,j,k))
          else
            ubaro(i,j)=0.0
          endif
c
          if (iv(i,j).eq.1) then
            vbaro(i,j)=v(i,j,2*k)*.5*(dp(i,j,k)+dp(i,jm1,k))
          else
            vbaro(i,j)=0.0
          endif
c
          if (iq(i,j).eq.1) then
            vort(i,j)=
     &       (vbaro(i,j)*scvy(i,j)-vbaro(im1,j)*scvy(im1,j)
     &       -ubaro(i,j)*scux(i,j)+ubaro(i,jm1)*scux(i  ,jm1))*1.e-6
          else
            vort(i,j)=0.0
          endif
        enddo 
      enddo 
c
c --- calculate first estimate for strmf values
      ubi=0.0
      do j=1,jj-1
        ubim1=ubi
        ubi  =ubaro(1,j)*scux(1,j)
        if     (j.eq.1) then
          strmfu=0.0
        else
          strmfu=strmfu-0.5*(ubi+ubim1)
        endif
        strmft=strmfu
        if     (iq(1,j).eq.1) then
          strmf(1,j)=strmft*1.e-6
        else
          strmf(1,j)=0.
        endif
        vbi=0.0
        do i=2,ii
          vbim1 =vbi
          vbi   =vbaro(i,j)*scvy(i,j)
          strmft=strmft+0.5*(vbi+vbim1)
          if     (iq(i,j).eq.1) then
            strmf(i,j)=strmft*1.e-6
          else
            strmf(i,j)=0.
          endif
        enddo
      enddo
      do i=1,ii
        strmf(i,jj)=0.
      enddo
c
      call poisnd(ii,jj,vort,strmf,util1,util1,work)
      call qsmoo(strmf,work)
c
c --- interpolate from q-grid to p-grid.
      do j=1,jj
        do i=1,ii
          if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
            s1=0.0
            s2=0.0
            if     (iq(i,  j)  .eq.1) then
              s1 = s1 + 1.0
              s2 = s2 + strmf(i,  j)
            endif
            if     (iq(i+1,j)  .eq.1) then
              s1 = s1 + 1.0
              s2 = s2 + strmf(i+1,j)
            endif
            if     (iq(i,  j+1).eq.1) then
              s1 = s1 + 1.0
              s2 = s2 + strmf(i,  j+1)
            endif
            if     (iq(i+1,j+1).eq.1) then
              s1 = s1 + 1.0
              s2 = s2 + strmf(i+1,j+1)
            endif
            if     (s1.ne.0.0) then
              util1(i,j)=s2/s1
            else
              util1(i,j)=flag
            endif
          else
            util1(i,j)=flag
          endif
        enddo
      enddo
ccc   call zebra(util1,ii,ii,jj)
ccc   write (*,'('' shown above: layer''i3'' stream function'')') k
c
      qq=max(.25,contur(util1,ii,ii,jj))
      if (qqin.ne.0.0) qq=qqin
      ipalet=kpalet
      if (kpalet.ge.2 .and. qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqcsfn,'center','("blkinr: ",a6," =",f11.4," ")')
        endif
        qqc=qqcsfn
        qqmn = qqc-0.5*qqin*cntrs(kpalet)
        qqmx = qqc+0.5*qqin*cntrs(kpalet)
      else
        qqmn = -15.5*qq
        qqmx =  15.5*qq
      endif
      call horplt(util1,plon,plat,ii,jj,ii,jj,qqmn,qqmx,qq,
     &            -495,label(33:81),0.,nquad,.true.,lalolb,lalogr)
        call pcloqu(xlab0,ylab0,region,    csb,0.,0)
        write (text,'(''ci '',f4.1,'' sv'')') qq
        call pcloqu(xlab1,ylab1,text(1:10),csn,0.,0)
        write (text,'(f6.1,'' to'',f7.1)') amn,amx
        call pcloqu(xlab2,ylab2,text(1:16),csn,0.,0)
      if (nquad.eq.0) call fram(ncount)
      endif
c
      if (nquad.gt.0) then
        call fram(ncount)
        nquad=0
      end if
      enddo  !k-loop
      enddo  !layer loop
c
      if (nquad.gt.0) call fram(ncount)
      nquad=0
c
c --- use environment variable CROSS_LABELS to provide the
c --- name of a file containing a 4-character label for each
c --- layer on cross-section plots.  Note that CROSS_LABELS="NONE"
c --- indicates no cross section plots are needed (i.e. no[ij]sec=0),
c --- and can speed up the processing in such cases.
c
      crflnm = ' '
      call getenv('CROSS_LABEL',crflnm)
c
      if (kkin.eq.1 .or. crflnm.eq."NONE") then
        call clsgks
        stop '(normal: no cross-sections)'
      endif
c
c --- ----------------------------
c --- c r o s s    s e c t i o n s
c --- ----------------------------
c
c --- 'topsec' = cross section plot top (OPTIONAL, default 0.0)
c --- 'depth ' = cross section plot depth
c --- 'tthovr' = isopycnal overlay contour int (=0.0 overlay interfaces)
c --- 'vstep ' = velocity contours (1.0 stairstep, to 0.0 gently curved)
c --- 'velqq ' = uvel contour int (<0 no vel  plot; 0 from field)
c --- 'velqq ' = vvel contour int (<0 no vel  plot; 0 from field) (OPTIONAL)
c --- 'center' = central contoured value (ignored if kpalet<2)
c --- 'wvelqq' = wvel contour int (<0 no vel  plot; 0 from field) (OPTIONAL)
c --- 'center' = central contoured value (ignored if kpalet<2)
c --- 'temqq ' = temp contour int (<0 no temp plot; 0 from field)
c --- 'center' = central contoured value (ignored if kpalet<2)
c --- 'salqq ' = saln contour int (<0 no saln plot; 0 from field)
c --- 'center' = central contoured value (ignored if kpalet<2)
c --- 'tthqq ' = dens contour int (<0 no dens plot; 0 from field)
c --- 'center' = central contoured value (ignored if kpalet<2)
c --- 'mxlflg' = plot mixed layer (0=no,1-2=mxl,3-4=bl,5-6=mxl&bl,even=smooth)
c --- 'kpalet' = palete (0=none,1=pastel,2:7=sst,gaudy,2tone,fc,ifc,fc20)
c ---            negative to also line contour the field (no interfaces)
c --- 
c --- if one 'velqq ' is provided it is for the normal velocity
c
c --- one section plot for each non-negative velqq,wvelqq,temqq,salqq,tthqq
c
      if     (crflnm.eq.' ') then
        do k= 1,kk
          if     (k.lt.10) then
            write(crlabl(k),'(i1)') k
          else
            write(crlabl(k),'(i2)') k
          endif
        enddo
      else
        open(unit=99,file=crflnm,status='old',form='formatted')
        do k= 1,kk
          read(99,'(a)') crlabl(k)
        enddo
        close(99)
      endif
c
      if     (baclin) then
        csec_name(1) = 'u-bcl.vel.'
        csec_name(2) = 'v-bcl.vel.'
      else !total velocity
        csec_name(1) = 'u-velocity'
        csec_name(2) = 'v-velocity'
      endif
      csec_name(3) = 'w-velocity'
      csec_name(4) = 'temperature'
      csec_name(5) = 'salinity'
      csec_name(6) = 'density'
c
      call blkinr2(qq,i, 'topsec','("blkinr: ",a6," =",f11.4," m")',
     &                   'depth ','("blkinr: ",a6," =",f11.4," m")')
      if     (i.eq.1) then !topsec
      topsec=nint(qq)
      call blkinr(depth, 'depth ','("blkinr: ",a6," =",f11.4," m")')
      depth=nint(depth-topsec)  !total depth input, but use depth from topsec
      else !depth
      topsec=0.0
      depth =nint(qq)
      endif
      call blkinr2(qq,i, 'tthovr','("blkinr: ",a6," =",f11.4," sig")',
     &                   'vstep ','("blkinr: ",a6," =",f11.4," ")'   )
      if     (i.eq.1) then !tthovr
      tthovr=qq
      call blkinr(vstep, 'vstep ','("blkinr: ",a6," =",f11.4," ")')
      else !vstep
      tthovr=0.0
      vstep =qq
      endif
      call blkinr(qqsec(1),
     &            'velqq ','("blkinr: ",a6," =",f11.4," cm/s")')
      call blkinr2(qq,i,
     &             'velqq ','("blkinr: ",a6," =",f11.4," cm/s")',
     &             'center','("blkinr: ",a6," =",f11.4," cm/s")')
      if     (i.eq.2) then !'center'
        plotnv    = .true.  ! plot normal velocity only
        qqsec(2)  = -1.0
        qcsec(1)  = qq
        qcsec(2)  = qq
      else !velqq
        plotnv    = .false. ! plot u and/or v velocity
        qqsec(2)  = qq
      call blkinr(qcsec(1),
     &            'center','("blkinr: ",a6," =",f11.4," cm/s")')
        qcsec(2)  = qcsec(1)
      endif
      call blkinr2(qq,i,
     &             'wvelqq','("blkinr: ",a6," =",f11.4," cm/s")',
     &             'temqq ','("blkinr: ",a6," =",f11.4," deg")' )
      if     (i.eq.1) then !'wvelqq'
        qqsec(3) = qq
      call blkinr(qcsec(3),
     &            'center','("blkinr: ",a6," =",f11.4," cm/s")')
      call blkinr(qqsec(4),
     &            'temqq ','("blkinr: ",a6," =",f11.4," deg")' )
      else !'temqq '
        qqsec(3) = -1.0 !turn off 'wvelqq'
        qcsec(3) =  0.0
        qqsec(4) =  qq
      endif
      call blkinr(qcsec(4),
     &            'center','("blkinr: ",a6," =",f11.4," deg")')
      call blkinr(qqsec(5),
     &            'salqq ','("blkinr: ",a6," =",f11.4," psu")')
      call blkinr(qcsec(5),
     &            'center','("blkinr: ",a6," =",f11.4," psu")')
      call blkinr(qqsec(6),
     &            'tthqq ','("blkinr: ",a6," =",f11.4," sig")')
      call blkinr(qcsec(6),
     &            'center','("blkinr: ",a6," =",f11.4," sig")')
      do ktr= 1,ntracr
        write(cline,'(a,a,a)')
     &    '("blkinr: ",a6," =",f11.4,"', trim(ctrc_units(ktr)), '")'
        call blkinr(qqsec(6+ktr),'trcqq ',cline)
        call blkinr(qcsec(6+ktr),'center',cline)
        csec_name(6+ktr) = ctrc_name(ktr)
      enddo
      call blkini(mxlflg,'mxlflg')
      call blkini(kpalet,'kpalet')
c
      ipalet=abs(kpalet)
      call colors	!  redefine color table
c
      plotuv = qqsec(1).ge.0.0 .or. qqsec(2).ge.0.0
      plotw  = qqsec(3).ge.0.0
      plotem = qqsec(4).ge.0.0
      plosal = qqsec(5).ge.0.0
      plotth = qqsec(6).ge.0.0
      if     (mxlflg.eq.0) then
        plotbl = .false.
        plotml = .false.
      elseif (mod(mxlflg,2).eq.1) then
        plotbl = mxlflg.eq.3 .or. mxlflg.eq.5
        plotml = mxlflg.eq.1 .or. mxlflg.eq.5
      else  ! smooth
        plotbl = mxlflg.eq.4 .or. mxlflg.eq.6
        plotml = mxlflg.eq.2 .or. mxlflg.eq.6
        call psmo1(dpbl,  work,p(1,1,kk+1))
        call psmo1(dpmixl,work,p(1,1,kk+1))
      endif
c
c --- qcsec (i.e. 'center') == 0.0 is a signal to plot anomalies.
      lsecanom = .false.
      do kp= 4,6+ntracr
        if     (qqsec(kp).ge.0.0 .and. qcsec(kp).eq.0.0) then
          lsecanom = .true.
          csec_name(kp) = trim(csec_name(kp)) // ".anom"
        endif
      enddo
c
      do k=1,kout
        pout(k)=topsec+float(k-1)/float(kout-1)*depth
      enddo
c
c --- define pressure at mass points (2 per layer) for use in grdtrns
c
      dpth=.01*depth
c
      do 26 j=1,jj1
      do 26 i=1,ii1
c
c --- 'vstep' controls the appearance of velocity contours
c ---         vstep = 1.0 produces stairstep type contour lines
c ---                    with discontinuities at the layer interfaces
c ---         vstep = 0.0 produces gently curved contour lines
      do 25 k=1,kk
      ptop=p(i,j,k  )
      pbot=p(i,j,k+1)
      pmid=.51*ptop+.49*pbot
      pedg=min(ptop+dpth,pmid)
      puv(i,j,2*k-1)=(1.-vstep)*pmid+vstep*pedg
      pmid=.51*pbot+.49*ptop
      pedg=max(pbot-dpth,pmid)
 25   puv(i,j,2*k  )=(1.-vstep)*pmid+vstep*pedg
c
c --- put uppermost depth point at sea surface and lowest point at bottom
      puv(i,j,1)=0.
      puv(i,j,2*kk)=p(i,j,kk+1)
c
      xyij(i,j,1)=float(i)
 26   xyij(i,j,2)=float(j)
c
      lrfi(1)=ii1
      lrfi(2)=jj1
c
      if (plotuv) then
c
c --- interpolate -u,v- to -p- points for use in grdtrns (2 per layer)
c
      do 15 k=1,kk
c
      do 2 j=1,jj1
      do 2 i=1,ii1
      u(i,j,2*k-1)=.5*(u(i,j,2*k)+u(i+1,j,2*k))
 2    v(i,j,2*k-1)=.5*(v(i,j,2*k)+v(i,j+1,2*k))
c
      do 15 j=1,jj1
      do 15 i=1,ii1
      u(i,j,2*k)=u(i,j,2*k-1)
 15   v(i,j,2*k)=v(i,j,2*k-1)
c
      do 16 k=1,kout
      do 16 j=1,jj
      do 16 i=1,ii
      utr(i,j,k)=flag
      vtr(i,j,k)=flag
 16   wtr(i,j,k)=flag
c
      call grdtrns(0,puv,u,dum,ii,jj,2*kk,coord,lrfi,utr,utr,ii,jj,
     &             kout,pout,lrfo,xyij)
      call grdtrns(0,puv,v,dum,ii,jj,2*kk,coord,lrfi,vtr,vtr,ii,jj,
     &             kout,pout,lrfo,xyij)
ccc   print 101,((utr(i,i,k),i=1,25),k=1,kout)
ccc   print 101,((vtr(i,i,k),i=1,25),k=1,kout)
 101  format (/(1x,25f5.1))
c
      endif
c
      if (plotw) then
      ttr(:,:,:)=flag
      call grdtrns(0,puv,w,dum,ii,jj,2*kk,coord,lrfi,wtr,wtr,ii,jj,
     &             kout,pout,lrfo,xyij)
ccc   print 101,((wtr(i,i,k),i=1,25),k=1,kout)
      endif
c
      if (plotem) then
      ttr(:,:,:)=flag
      call grdtrns(0,puv,temp,dum,ii,jj,2*kk,coord,lrfi,ttr,ttr,
     &             ii,jj,kout,pout,lrfo,xyij)
ccc   print 102,((int(100.*ttr(i,i,k)),i=ii/2-7,ii/2+7),k=1,kout)
 102  format (//(1x,15i5))
      endif
c
      if (plosal) then
      str(:,:,:)=flag
      call grdtrns(0,puv,saln,dum,ii,jj,2*kk,coord,lrfi,str,str,
     &             ii,jj,kout,pout,lrfo,xyij)
ccc   print 102,((int(100.*str(i,i,k)),i=ii/2-7,ii/2+7),k=1,kout)
      endif
c
      if (plotth.or.tthovr.ne.0.0) then
      thtr(:,:,:)=flag
      call grdtrns(0,puv,th3d,dum,ii,jj,2*kk,coord,lrfi,thtr,thtr,
     &             ii,jj,kout,pout,lrfo,xyij)
ccc   print 102,((int(100.*thtr(i,i,k)),i=ii/2-7,ii/2+7),k=1,kout)
      endif
c
      do ktr= 1,ntracr
        if (qqsec(6+ktr).ge.0.0) then
          write(lp,*) 'trcr = ',trcr(6,11,1:2*kk,ktr)
          trtr(:,:,:,ktr)=flag
          call grdtrns(0,puv,trcr(1,1,1,ktr),dum,ii,jj,2*kk,
     &                 coord,lrfi,trtr(1,1,1,ktr),trtr(1,1,1,ktr),
     &                 ii,jj,kout,pout,lrfo,xyij)
          write(lp,*) 'trtr = ',trtr(6,11,1:kout,ktr)
ccc       print 102,((int(100.*trtr(i,i,k,ktr)),
ccc  &                i=ii/2-7,ii/2+7),k=1,kout)
        endif
      enddo
c
c --- 'noisec' = number of i cross sections
      call blkini(noisec,'noisec')
      nsplot=0
      do 4 ni= 1,noisec
c --- 'isec  ' = i cross section location, or
c --- 'i1st  ' = i cross section location at i=1
      call blkini2(i,j,  'isec  ','i1st  ')  !read isec or i1st
      if     (j.eq.1) then !isec
        if     (i.lt.1 .or. i.gt.ii1) then
          write(lp,*)
          write(lp,*) 'error - isec must be between 1 and',ii1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        i1st = i
        iend = i
c
        alon = mod(plon(i,jj/2)+1080.0,360.0) !guess mid-point is result
        ij   = 0
        do j=1,jj1
          if     (ip(i,j).eq.1) then
            blon = mod(plon(i,j)+1080.0,360.0)
            if     (blon-alon.gt. 180.0) then
              blon = blon - 360.0
            elseif (blon-alon.lt.-180.0) then
              blon = blon + 360.0
            endif
            ij = ij + 1
            xlonlat(ij) = blon
*           write(lp,'(a,i5,2f10.2)') 
*    &        'j,alon,blon = ',j,alon,blon
          endif
        enddo
        if     (ij.eq.0) then
          write (lp,'('' skip section at  i ='',i4,'' (all land)'')') i
          goto 4
        endif
        call ssort(xlonlat,xlonlat,ij,1)  ! sort xlonlat to get median
        xlonlat0 = mod(xlonlat((ij+1)/2)+1260.0,360.0)-180.0
        xlonlat1 = xlonlat0
        write (lp,'(" next section at  i,lon=",i5,f10.2)') i,xlonlat0
        call flush(lp)
      else  !i1st (diagonal line)
        i1st = i
c ---   'iend  ' = i cross section location at j=jdmp-1
        call blkini(iend,'iend  ')
        if     (i1st.lt.1 .or. i1st.gt.ii1) then
          write(lp,*)
          write(lp,*) 'error - i1st must be between 1 and',ii1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        if     (iend.lt.1 .or. iend.gt.ii1) then
          write(lp,*)
          write(lp,*) 'error - iend must be between 1 and',ii1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
c
        ij = 0
        do j=1,jj1
          i = i1st + nint((iend-i1st)*float(j-1)/float(jj1-1))
          if     (ip(i,j).eq.1) then
            ij = ij + 1
          endif
          if     (j.eq.jj1 .or. mod(j,max(3,jj1/5)).eq.1) then
            write(lp,'(a,2i5)') '    diag: i,j =',i,j
          endif
        enddo
        if     (ij.eq.0) then
          write (lp,'('' skip section at  i ='',i4,'' (all land)'')') i
          goto 4
        endif
        xlonlat0 = mod(plon(i1st,  1)+1260.0,360.0)-180.0
        xlonlat1 = mod(plon(iend,jj1)+1260.0,360.0)-180.0
        write (lp,'(" next section at  i,lon=",i5,f10.2,
     &              " to i,lon=",i5,f10.2)')
     &    i1st,xlonlat0,iend,xlonlat1
        call flush(lp)
      endif  !vertical or diagonal line
c
      do j=1,jj1
        i = i1st + nint((iend-i1st)*float(j-1)/float(jj1-1))
        if     (ip(i,j).eq.1) then
          dpbl2d(j) =   dpbl(i,j) - topsec
          dpml2d(j) = dpmixl(i,j) - topsec
          do k=1,kk+1
            p2d(j,k)=p(i,j,k) - topsec
          enddo !k
        else
          dpbl2d(j) = flag
          dpml2d(j) = flag
          do k=1,kk+1
            p2d(j,k)=flag
          enddo !k
        endif !ip
        sumsec(4:6+ntracr) = 0.d0
        samsec(4:6+ntracr) = 0.d0
        do k=1,kout
          if     (ip(i,j).eq.1) then
            if     (utr(i,j,kout+1-k).ne.flag) then
              uv2d(j,k)= utr(i,j,kout+1-k)*100.0
            else
              uv2d(j,k)= flag
            endif
            if     (vtr(i,j,kout+1-k).ne.flag) then
              vv2d(j,k)= vtr(i,j,kout+1-k)*100.0
            else
              vv2d(j,k)= flag
            endif
            if     (wtr(i,j,kout+1-k).ne.flag) then
              wv2d(j,k)= wtr(i,j,kout+1-k)*8640000.0
            else
              wv2d(j,k)= flag
            endif
            tm2d(j,k)= ttr(i,j,kout+1-k)
            sl2d(j,k)= str(i,j,kout+1-k)
            th2d(j,k)=thtr(i,j,kout+1-k)
            if     (tm2d(j,k).ne.flag) then
              sumsec(4) = sumsec(4)+tm2d(j,k)
              samsec(4) = samsec(4)+1.d0
            endif
            if     (sl2d(j,k).ne.flag) then
              sumsec(5) = sumsec(5)+sl2d(j,k)
              samsec(5) = samsec(5)+1.d0
            endif
            if     (th2d(j,k).ne.flag) then
              sumsec(6) = sumsec(6)+th2d(j,k)
              samsec(6) = samsec(6)+1.d0
            endif
            do ktr= 1,ntracr
              tr2d(j,k,ktr)=trtr(i,j,kout+1-k,ktr)
              if     (tr2d(j,k,ktr).ne.flag) then
                sumsec(6+ktr) = sumsec(6+ktr)+tr2d(j,k,ktr)
                samsec(6+ktr) = samsec(6+ktr)+1.d0
              endif
*             if     (j.eq.11) then
*               write(lp,*) 'k,tr2d = ',k,tr2d(j,k,ktr)
*             endif
            enddo
          else
            tm2d(j,k)=flag
            sl2d(j,k)=flag
            uv2d(j,k)=flag
            vv2d(j,k)=flag
            wv2d(j,k)=flag
            th2d(j,k)=flag
            do ktr= 1,ntracr
              tr2d(j,k,ktr)=flag
            enddo
          endif !ip
        enddo !k
        if     (lsecanom) then  !subtract the mean from at least one section
          do kp=4,6+ntracr
            sumsec(kp) = sumsec(kp)/max(samsec(kp),1.d0)
          enddo
          do k=1,kout
            if     (ip(i,j).eq.1) then
              if     (qcsec(4).eq.0.0 .and. tm2d(j,k).ne.flag) then
                tm2d(j,k)= tm2d(j,k) - sumsec(4)
              endif
              if     (qcsec(5).eq.0.0 .and. sl2d(j,k).ne.flag) then
                sl2d(j,k)= sl2d(j,k) - sumsec(5)
              endif
              if     (qcsec(6).eq.0.0 .and. th2d(j,k).ne.flag) then
                th2d(j,k)= th2d(j,k) - sumsec(6)
              endif
              do ktr= 1,ntracr
                if     (qcsec(6+ktr).eq.0.0 .and.
     &                  tr2d(j,k,ktr).ne.flag) then
                  tr2d(j,k,ktr)=tr2d(j,k,ktr) - sumsec(6+ktr)
                endif
              enddo !ktr
            endif !ip.eq.1
          enddo !k
        endif !lsecanom
      enddo !j
      i = (i1st+iend)/2
c
ccc      do iter=1,2
ccc      if (plotem) call filtr1(tm2d,work,iout,jj1,kout)
ccc      if (plosal) call filtr1(sl2d,work,iout,jj1,kout)
ccc      if (plotth) call filtr1(th2d,work,iout,jj1,kout)
ccc      end do
c
      nsplot=0
      do kp= 1,6+ntracr
      if     (qqsec(kp).ge.0.0) then
        qqin = qqsec(kp)
        if (abs(kpalet).ge.2 .and. qqin.gt.0.0) then
          qqlin = qcsec(kp)-0.5*qqin*cntrs(abs(kpalet))
          qqhin = qcsec(kp)+0.5*qqin*cntrs(abs(kpalet))
        else
          qqlin = 0.0
          qqhin = 0.0
        endif
        nsplot=nsplot+1
      else
        cycle
      endif
c
c --- draw sections along lines i = const. from ja to jb
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 1: split section in two
ccc      do 4 ja=1,jj/2,jj/2-1
ccc      jb=ja+jj/2-1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 2: display section in one piece
      ja=1
      jb=jj1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if     (kp.eq.1) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),uv2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.2) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),vv2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.3) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),wv2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.4) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),tm2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.5) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),sl2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.6) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),th2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      else
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),tr2d(1,1,kp-6),iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      endif
c
      call getset(xc,xd,yc,yd,xa,xb,ya,yb,dum)
      call set(xc,xd,yc,yd,float(ja),float(jb),ya,yb,dum)
      chrsiz=.008 * (yb-ya)/(yd-yc)
      if     (lalolb.gt.0) then
c ---   add latitude labels
        alalolb=lalolb
        do j=1,jj1
          i = i1st + nint((iend-i1st)*float(j-1)/float(jj1-1))
          alat = sign(alalolb,plat(i,j))*nint(abs(plat(i,j))/alalolb)
          if     (abs(plat(i,j+1)-plat(i,j)).gt.1.e-4) then
            q=(alat-plat(i,j))/(plat(i,j+1)-plat(i,j))
            if     (q.ge.0.0 .and. q.lt.1.0) then  !found a needed latitude
              x   = j+q
              lat = nint(alat)
              if     (lat.eq.0.0) then
                write (label(1:3),'(   a3)') ' EQ'
              elseif (lat.gt.0.0) then
                write (label(1:3),'(i2,a1)')  lat,'N'
              else
                write (label(1:3),'(i2,a1)') -lat,'S'
              endif
              call pcloqu(x,ya+2.*chrsiz,label(1:3),csn,0.,0.)
*             if (ni.eq.1) then
*               write(lp,'(a,a,f9.3)') 'lat,x =  ',label(1:3),x
*             endif
            endif !q
          endif !plat
        enddo !j
      elseif (lalolb.lt.0) then
c ---   add array index labels
        do j= ja-lalolb-mod(ja,-lalolb),jb-mod(jb,-lalolb),-lalolb
          x=j
          write (label(1:4),'(i4)') j
          call pcloqu(x,ya+2.*chrsiz,label(1:4),csn,0.,0.)
          if (ni.eq.1) then
            write(lp,'(a,a,f9.3)') 'label,x =  ',label(1:3),x
          endif
        enddo !i
      endif !lalolb
      if     (region.ne.'        ') then
        x=ja+0.08*(jb-ja)
        y=ya+5.0*chrsiz
        call pcloqu(x,y,region, csb,0.,0.)
*       write(lp,'(a,a,f9.3)') region,',x = ',x
      endif
      if (nquad.eq.0) call fram(ncount)       !2 plots per page
      enddo !kp
      if (nquad.gt.0 .and. nsplot.ge.3) then  !3 or more plots per section
        call fram(ncount)  !start each section on a new page
        nquad=0
      endif
  4   continue !ni
c
c --- 'nojsec' = number of j sections
      call blkini(nojsec,'nojsec')
c --- use a single frame if there is exactly one i and one j section plot.
      if (nsplot.gt.1 .or. noisec.ne.1 .or. nojsec.ne.1) then
        if (nquad.gt.0) call fram(ncount)
        nquad=0
      endif
      if     (plotnv) then
        qqsec(2) = qqsec(1)  ! normal is v-velocity
        qqsec(1) = -1.0
      endif
      do 5 nj= 1,nojsec
c --- 'jsec  ' = j cross section location, or
c --- 'j1st  ' = j cross section location at i=1
      call blkini2(j,i,  'jsec  ','j1st  ')  !read jsec or j1st
      if     (i.eq.1) then  !jsec
        if     (j.lt.1 .or. j.gt.jj1) then
          write(lp,*)
          write(lp,*) 'error - jsec must be between 1 and',jj1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        j1st = j
        jend = j
c
        ij    = 0
        do i=1,ii1
          if     (ip(i,j).eq.1) then
            ij = ij + 1
            xlonlat(ij) = plat(i,j)
          endif
        enddo !i
        if     (ij.eq.0) then
          write (lp,'('' skip section at  j='',i4,'' (all land)'')') j
          goto 5
        endif
        call ssort(xlonlat,xlonlat,ij,1)  ! sort xlonlat to get median
        xlonlat0 = xlonlat((ij+1)/2)
        xlonlat1 = xlonlat0
        write (lp,'(" next section at  j,lat=",i5,f10.2)') j,xlonlat0
      else  !j1st (diagonal line)
        j1st = j
c ---   'jend  ' = j cross section location at i=idmp-1
        call blkini(jend,'jend  ')
        if     (j1st.lt.1 .or. j1st.gt.jj1) then
          write(lp,*)
          write(lp,*) 'error - j1st must be between 1 and',jj1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        if     (jend.lt.1 .or. jend.gt.jj1) then
          write(lp,*)
          write(lp,*) 'error - jend must be between 1 and',jj1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
c
        ij = 0
        do i=1,ii1
          j = j1st + nint((jend-j1st)*float(i-1)/float(ii1-1))
          if     (ip(i,j).eq.1) then
            ij = ij + 1
          endif
          if     (i.eq.ii1 .or. mod(i,max(3,ii1/5)).eq.1) then
            write(lp,'(a,2i5)') '    diag: i,j =',i,j
          endif
        enddo !i
        if     (ij.eq.0) then
          write (lp,'('' skip section at  j ='',i4,'' (all land)'')') j
          goto 5
        endif
        xlonlat0 = plat(1,  j1st)
        xlonlat1 = plat(ii1,jend)
        write (lp,'(" next section at  j,lat=",i5,f10.2,
     &              " to j,lat=",i5,f10.2)')
     &    j1st,xlonlat0,jend,xlonlat1
        call flush(lp)
      endif  !vertical or diagonal line
c
      do i=1,ii1
        j = j1st + nint((jend-j1st)*float(i-1)/float(ii1-1))
        if     (ip(i,j).eq.1) then
          dpbl2d(i) =   dpbl(i,j) - topsec
          dpml2d(i) = dpmixl(i,j) - topsec
          do k=1,kk+1
            p2d(i,k)=p(i,j,k) - topsec
          enddo !k
        else
          dpbl2d(i) = flag
          dpml2d(i) = flag
          do k=1,kk+1
            p2d(i,k)=flag
          enddo !k
        endif !ip
        sumsec(4:6+ntracr) = 0.d0
        samsec(4:6+ntracr) = 0.d0
        do k=1,kout
          if     (ip(i,j).eq.1) then
            if     (utr(i,j,kout+1-k).ne.flag) then
              uv2d(i,k)= utr(i,j,kout+1-k)*100.0
            else
              uv2d(i,k)= flag
            endif
            if     (vtr(i,j,kout+1-k).ne.flag) then
              vv2d(i,k)= vtr(i,j,kout+1-k)*100.0
            else
              vv2d(i,k)= flag
            endif
            if     (wtr(i,j,kout+1-k).ne.flag) then
              wv2d(i,k)= wtr(i,j,kout+1-k)*8640000.0
            else
              wv2d(i,k)= flag
            endif
            tm2d(i,k)= ttr(i,j,kout+1-k)
            sl2d(i,k)= str(i,j,kout+1-k)
            th2d(i,k)=thtr(i,j,kout+1-k)
            if     (tm2d(i,k).ne.flag) then
              sumsec(4) = sumsec(4)+tm2d(i,k)
              samsec(4) = samsec(4)+1.d0
            endif
            if     (sl2d(i,k).ne.flag) then
              sumsec(5) = sumsec(5)+sl2d(i,k)
              samsec(5) = samsec(5)+1.d0
            endif
            if     (th2d(i,k).ne.flag) then
              sumsec(6) = sumsec(6)+th2d(i,k)
              samsec(6) = samsec(6)+1.d0
            endif
            do ktr= 1,ntracr
              tr2d(i,k,ktr)=trtr(i,j,kout+1-k,ktr)
              if     (tr2d(j,k,ktr).ne.flag) then
                sumsec(6+ktr) = sumsec(6+ktr)+tr2d(j,k,ktr)
                samsec(6+ktr) = samsec(6+ktr)+1.d0
              endif
*             if     (i.eq.6) then
*               write(lp,*) 'k,tr2d = ',k,tr2d(i,k,ktr)
*             endif
            enddo
          else
            tm2d(i,k)=flag
            sl2d(i,k)=flag
            uv2d(i,k)=flag
            vv2d(i,k)=flag
            wv2d(i,k)=flag
            th2d(i,k)=flag
            do ktr= 1,ntracr
              tr2d(i,k,ktr)=flag
            enddo
          endif !ip
        enddo !k
        if     (lsecanom) then  !subtract the mean from at least one section
          do kp=4,6+ntracr
            sumsec(kp) = sumsec(kp)/max(samsec(kp),1.d0)
          enddo
          do k=1,kout
            if     (ip(i,j).eq.1) then
              if     (qcsec(4).eq.0.0 .and. tm2d(i,k).ne.flag) then
                tm2d(i,k)= tm2d(i,k) - sumsec(4)
              endif
              if     (qcsec(5).eq.0.0 .and. sl2d(i,k).ne.flag) then
                sl2d(i,k)= sl2d(i,k) - sumsec(5)
              endif
              if     (qcsec(6).eq.0.0 .and. th2d(i,k).ne.flag) then
                th2d(i,k)= th2d(i,k) - sumsec(6)
              endif
              do ktr= 1,ntracr
                if     (qcsec(6+ktr).eq.0.0 .and.
     &                  tr2d(i,k,ktr).ne.flag) then
                  tr2d(i,k,ktr)=tr2d(i,k,ktr) - sumsec(6+ktr)
                endif
              enddo !ktr
            endif !ip.eq.1
          enddo !k
        endif !lsecanom
      enddo !i
      j = (j1st+jend)/2
c
ccc      do iter=1,2
ccc      if (plotem) call filtr1(tm2d,work,iout,ii1,kout)
ccc      if (plosal) call filtr1(sl2d,work,iout,ii1,kout)
ccc      if (plotth) call filtr1(th2d,work,iout,ii1,kout)
ccc      end do
c
      nsplot=0
      do kp= 1,6+ntracr
      if     (qqsec(kp).ge.0.0) then
        qqin = qqsec(kp)
        if (abs(kpalet).ge.2 .and. qqin.gt.0.0) then
          qqlin = qcsec(kp)-0.5*qqin*cntrs(abs(kpalet))
          qqhin = qcsec(kp)+0.5*qqin*cntrs(abs(kpalet))
        else
          qqlin = 0.0
          qqhin = 0.0
        endif
        nsplot=nsplot+1
      else
        cycle
      endif
c
c --- draw sections along lines j = const. from ia to ib
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 1: split section in two
ccc      do 5 ia=1,ii/2,ii/2-1
ccc      ib=ia+ii/2-1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 2: display section in one piece
      ia=1
      ib=ii1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if     (kp.eq.1) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),uv2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.2) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),vv2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.3) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),wv2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.4) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),tm2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.5) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),sl2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.6) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),th2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      else
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),tr2d(1,1,kp-6),iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      endif
c
      call getset(xc,xd,yc,yd,xa,xb,ya,yb,dum)
      call set(xc,xd,yc,yd,float(ia),float(ib),ya,yb,dum)
      chrsiz=.008 * (yb-ya)/(yd-yc)
c
      if     (lalolb.gt.0) then
c ---   add longitude labels
        alalolb=lalolb
        do i=1,ii1
          j = j1st + nint((jend-j1st)*float(i-1)/float(ii1-1))
          alon  = sign(alalolb,plon(i,j))*nint(abs(plon(i,j))/alalolb)
          plon1 = plon(i+1,j)
          q     = plon1 - plon(i,j)
          if     (q.ge. 720.0) then
            plon1 = plon1 -720.0
          elseif (q.ge. 360.0) then
            plon1 = plon1 -360.0
          elseif (q.le.-720.0) then
            plon1 = plon1 + 720.0
          elseif (q.le.-360.0) then
            plon1 = plon1 + 360.0
          endif
          if     (abs(plon1-plon(i,j)).gt.1.e-4) then
            q=(alon-plon(i,j))/(plon1-plon(i,j))
            if     (q.ge.0.0 .and. q.lt.1.0) then  !found a needed longitude
              x     = i+q
              lonew = mod(nint(alon)+1260,360)-180
              if     (lonew.ge.0.0) then
                write (label(1:4),'(i3,a1)')  lonew,'E'
              else
                write (label(1:4),'(i3,a1)') -lonew,'W'
              endif
              if      (label(2:2).eq.' ') then
                call pcloqu(x,ya+2.*chrsiz,label(3:4),csn,0.,0.)
              else if (label(1:1).eq.' ') then
                call pcloqu(x,ya+2.*chrsiz,label(2:4),csn,0.,0.)
              else
                call pcloqu(x,ya+2.*chrsiz,label(1:4),csn,0.,0.)
              endif
*             if (nj.eq.1) then
*               write(lp,'(a,a,f9.3)') 'lon,x = ',label(1:4),x
*             endif
            endif !q
          endif !plon
        enddo !i
      elseif (lalolb.lt.0) then
c ---   add array index labels
        do i= ia-lalolb-mod(ia,-lalolb),ib-mod(ib,-lalolb),-lalolb
          x=i
          write (label(1:4),'(i4)')  i
          if      (label(3:3).eq.' ') then
            call pcloqu(x,ya+2.*chrsiz,label(4:4),csn,0.,0.)
          else if (label(2:2).eq.' ') then
            call pcloqu(x,ya+2.*chrsiz,label(3:4),csn,0.,0.)
          else if (label(1:1).eq.' ') then
            call pcloqu(x,ya+2.*chrsiz,label(2:4),csn,0.,0.)
          else
            call pcloqu(x,ya+2.*chrsiz,label(1:4),csn,0.,0.)
          endif
          if (nj.eq.1) then
            write(lp,'(a,a,f9.1)') 'label,x = ',label(1:4),x
          endif
        enddo !i
      endif !lalolb
      if     (region.ne.'        ') then
        if     (loclab.eq.3 .or. loclab.eq.4) then
          x=ia+0.08*(ib-ia)  ! lower-left
        else
          x=ib-0.08*(ib-ia)  ! lower-right
        endif
        y=ya+5.0*chrsiz
        call pcloqu(x,y,region, csb,0.,0.)
        if (nj.eq.1) write(lp,'(a,a,f9.3)') region,',x = ',x
      endif
      if (nquad.eq.0) call fram(ncount)       !2 plots per page
      enddo !kp
      if (nquad.gt.0 .and. nsplot.ge.3) then  !3 or more plots per section
        call fram(ncount)  !start each section on a new page
        nquad=0
      endif
  5   continue !nj
      if (nquad.gt.0) call fram(ncount)
      nquad=0
c
      call clsgks
      stop '(normal)'
      end
