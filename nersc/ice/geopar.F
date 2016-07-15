      subroutine geopar
      use mod_xc         ! HYCOM communication interface
      !KAL 20160715 use mod_cb_arrays  ! HYCOM saved arrays
      use mod_za         ! HYCOM I/O interface
      use mod_common_ice, only : 
     &   plon,plat,
     &   ulon,ulat,
     &   vlon,vlat,
     &   qlon,qlat,
     &   scpx,scpy,
     &   scux,scuy,
     &   scvx,scvy,
     &   scqx,scqy,
     &   corio, depths,
     &   scu2,  scv2,
     &   scp2,  scq2,
     &   pang,pasp,
     &   flnmgrd,
     &   flnmdep,
     &   scu2,scv2,
     &   scp2,scq2,
     &   scp2i,scq2i,
     &   scuxi,scvyi,
     &   aspux, aspuy,
     &   aspvx, aspvy,
     &   scqx,scqy, 
     &   mapflg,itest,jtest, epsil, hugel,
     &   util1,util2,util3
#if defined(USE_CCSM3)
      use ccsm3_grid, only : ANGLET
#endif
c
c --- set up model parameters related to geography
c
c --- hycom version 2.1
c
      implicit none
c
      real      dp0kf,dpm,dpms,ds0kf,dsm,dsms
      real      hmina,hminb,hmaxa,hmaxb
      real*8    sum_ip,sum_is,sum_isa
      integer   i,ios,j,k,ktr,l,nishlf
      character preambl(5)*79,cline*80
#if defined(USE_CCSM3)
      real      plinei(itdm),plinej(jtdm)
      save      plinei,plinej
#endif
c
      real       aspmax
      parameter (aspmax=2.0)  ! maximum grid aspect ratio for diffusion
*     parameter (aspmax=1.0)  ! ignore  grid aspect ratio in  diffusion
c
c --- read grid location,spacing,coriolis arrays
c
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        write (lp,'(3a)') ' reading grid file from ',
     &                         trim(flnmgrd),'.[ab]'
        open (unit=uoff+9,file=trim(flnmgrd)//'.b',
     &        status='old')
      endif
      call xcsync(flush_lp)
      call zagetc(cline,ios, uoff+9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      read(cline,*) i
c
      call zagetc(cline,ios, uoff+9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      read (cline,*) j
c
      if     (i.ne.itdm .or. j.ne.jtdm) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)')
     &    'error - wrong array size in grid file'
        endif
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      call zagetc(cline,ios, uoff+9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      if     (mnproc.eq.1) then
      write (lp,'(a)') trim(cline)
      endif
      read (cline,*) mapflg
c
      call zaiopf(trim(flnmgrd)//'.a','old', 9)
c
      do k= 1,15
        call zagetc(cline,ios, uoff+9)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)')
     &        'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
          endif !1st tile
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
        i = index(cline,'=')
        read (cline(i+1:),*) hminb,hmaxb
        if     (mnproc.eq.1) then
        write (lp,'(a)') trim(cline)
        endif
        call xcsync(flush_lp)
c
        if     (k.eq.1) then
          call zaiord(plon, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.2) then
          call zaiord(plat, ip,.false., hmina,hmaxa, 9)
          do i= 1,2  !skip qlon,qlat
            call zagetc(cline,ios, uoff+9)
            if     (ios.ne.0) then
              if     (mnproc.eq.1) then
                write(lp,'(/ a,i4,i9 /)')
     &            'geopar: I/O error from zagetc, iunit,ios = ',
     &            uoff+9,ios
              endif !1st tile
              call xcstop('(geopar)')
                     stop '(geopar)'
            endif
            call zaiosk(9)
          enddo
        elseif (k.eq.3) then
          call zaiord(ulon, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.4) then
          call zaiord(ulat, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.5) then
          call zaiord(vlon, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.6) then
          call zaiord(vlat, ip,.false., hmina,hmaxa, 9)
          call zagetc(cline,ios, uoff+9)
          if     (ios.ne.0) then
            if     (mnproc.eq.1) then
              write(lp,'(/ a,i4,i9 /)')
     &          'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
            endif !1st tile
            call xcstop('(geopar)')
                   stop '(geopar)'
          endif
#if defined(USE_CCSM3)
c         pang in ANGLET
          i = index(cline,'=')
          read (cline(i+1:),*) hminb,hmaxb
          if     (mnproc.eq.1) then
          write (lp,'(a)') trim(cline)
          endif
          call xcsync(flush_lp)
          call zaiord(ANGLET, ip,.false., hmina,hmaxa, 9)
#elif defined(ENABLE_ATM)
c         pang
          i = index(cline,'=')
          read (cline(i+1:),*) hminb,hmaxb
          if     (mnproc.eq.1) then
          write (lp,'(a)') trim(cline)
          endif
          call xcsync(flush_lp)
          call zaiord(pang, ip,.false., hmina,hmaxa, 9)
#else
c         skip pang
          call zaiosk(9)
#endif
        elseif (k.eq.7) then
          call zaiord(scpx, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.8) then
          call zaiord(scpy, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.9) then
          call zaiord(scqx, iq,.false., hmina,hmaxa, 9)
        elseif (k.eq.10) then
          call zaiord(scqy, iq,.false., hmina,hmaxa, 9)
        elseif (k.eq.11) then
          call zaiord(scux, iu,.false., hmina,hmaxa, 9)
        elseif (k.eq.12) then
          call zaiord(scuy, iu,.false., hmina,hmaxa, 9)
        elseif (k.eq.13) then
          call zaiord(scvx, iv,.false., hmina,hmaxa, 9)
        elseif (k.eq.14) then
          call zaiord(scvy, iv,.false., hmina,hmaxa, 9)
        else
          call zaiord(corio,iq,.false., hmina,hmaxa, 9)
        endif
c
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          endif
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
      enddo
c
      call zaiocl(9)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        close(unit=uoff+9)
      endif
c
      if (itest.gt.0 .and. jtest.gt.0) then
        i=itest
        j=jtest
        write (lp,'(/ a,2i5,a,f8.3,a,f12.9,2f10.2/)')
     &   ' i,j=',i+i0,j+j0,
     &   ' plat=',plat(i,j),
     &   ' corio,scux,vy=',corio(i,j),scux(i,j),scvy(i,j)
      endif
      call xcsync(flush_lp)
#if defined(USE_CCSM3)
c --- printout similar to ccsm ice model
      call xclget(plinei,itdm, plon, 1,1, +1, 0, 1)
      call xclget(plinej,jtdm, plat, 1,1,  0,+1, 1)
      if     (mnproc.eq.1) then
        write (lp,*)
        write (lp,'(a,4f9.3,a,4f9.3)')
     &    '(domain) plon(:,1): ',
     &    plinei(1:4),' ...', plinei(itdm-2:itdm)
        write (lp,'(a,4f9.3,a,4f9.3)')
     &    '(domain) plat(1,:): ',
     &    plinej(1:4),' ...', plinej(jtdm-2:jtdm)
        write (lp,*)
      endif
      call xcsync(flush_lp)
#endif
c
c --- read basin depth array
c
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        write (lp,'(3a)') ' reading bathymetry file from ',
     &                         trim(flnmdep),'.[ab]'
        open (unit=uoff+9,file=trim(flnmdep)//'.b',
     &        status='old')
        read (     uoff+9,'(a79)')  preambl
      endif
      call xcsync(flush_lp)
      call zagetc(cline,ios, uoff+9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        close(unit=uoff+9)
        write (lp,'(/(1x,a))') preambl,cline
      endif
c
      call zaiopf(trim(flnmdep)//'.a','old', 9)
      call zaiord(depths,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        endif
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.0.5*hugel) then
            depths(i,j) = 0.0
          endif
        enddo
      enddo
c
c --- determine do-loop limits for u,v,p,q points, and update halo for depths
      call bigrid(depths, mapflg, util1,util2,util3)
ccc      call prtmsk(ip,depths,util1,idm,ii,jj,0.0,1.0,
ccc     &     'bottom depth (m)')
c
c     now safe to apply halo to arrays.
c
      vland = 1.0
      call xctilr(plon,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(plat,  1,1, nbdy,nbdy, halo_ps)
#if defined(USE_CCSM3)
      call xctilr(ANGLET,1,1, nbdy,nbdy, halo_ps)
#endif
#if defined(ENABLE_ATM)
      call xctilr(pang,  1,1, nbdy,nbdy, halo_ps)
#endif
      call xctilr(scpx,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(scpy,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(ulon,  1,1, nbdy,nbdy, halo_us)
      call xctilr(ulat,  1,1, nbdy,nbdy, halo_us)
      call xctilr(scux,  1,1, nbdy,nbdy, halo_us)
      call xctilr(scuy,  1,1, nbdy,nbdy, halo_us)
      call xctilr(vlon,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(vlat,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(scvx,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(scvy,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(corio, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scqx,  1,1, nbdy,nbdy, halo_qs)
      call xctilr(scqy,  1,1, nbdy,nbdy, halo_qs)
      vland = 0.0
c
c --- area of grid cells (length x width) at u,v,p,q points resp.
c
******!$OMP PARALLEL DO PRIVATE(j,i)
******!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          scu2(i,j)=scux(i,j)*scuy(i,j)
          scv2(i,j)=scvx(i,j)*scvy(i,j)
          scp2(i,j)=scpx(i,j)*scpy(i,j)
          scq2(i,j)=scqx(i,j)*scqy(i,j)
c
          scuxi(i,j)=1.0/max(scux(i,j),epsil)
          scvyi(i,j)=1.0/max(scvy(i,j),epsil)
          scp2i(i,j)=1.0/max(scp2(i,j),epsil)
          scq2i(i,j)=1.0/max(scq2(i,j),epsil)
c
c ---     largest grid spacing (within limits) used in all diffusion
c ---     coefficients: min(max(sc?x,sc?y),sc?x*aspmax,sc?y*aspmax)
          aspux(i,j)=min(max(scux(i,j),scuy(i,j)),
     &                   min(scux(i,j),scuy(i,j))*aspmax)
     &               /max(scux(i,j),epsil)
          aspuy(i,j)=min(max(scux(i,j),scuy(i,j)),
     &                   min(scux(i,j),scuy(i,j))*aspmax)
     &               /max(scuy(i,j),epsil)
          aspvx(i,j)=min(max(scvx(i,j),scvy(i,j)),
     &                   min(scvx(i,j),scvy(i,j))*aspmax)
     &               /max(scvx(i,j),epsil)
          aspvy(i,j)=min(max(scvx(i,j),scvy(i,j)),
     &                   min(scvx(i,j),scvy(i,j))*aspmax)
     &               /max(scvy(i,j),epsil)
c
          util1(i,j)=depths(i,j)*scp2(i,j)
        enddo
      enddo

CKAL
CKALc
CKALc --- read ice shelf depth array
CKALc
CKAL      if     (ishelf.eq.0) then
CKAL        ishlf(:,:) = ip(:,:)  !no ice shelf
CKAL      else
CKAL        if     (mnproc.eq.1) then  ! .b file from 1st tile only
CKAL          write (lp,'(3a)') ' reading ice shelf file from ',
CKAL     &                           trim(flnmshlf),'.[ab]'
CKAL          open (unit=uoff+9,file=trim(flnmshlf)//'.b',
CKAL     &          status='old')
CKAL          read (     uoff+9,'(a79)')  preambl
CKAL        endif
CKAL        call xcsync(flush_lp)
CKAL        call zagetc(cline,ios, uoff+9)
CKAL        if     (ios.ne.0) then
CKAL          if     (mnproc.eq.1) then
CKAL            write(lp,'(/ a,i4,i9 /)')
CKAL     &        'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
CKAL          endif !1st tile
CKAL          call xcstop('(geopar)')
CKAL                 stop '(geopar)'
CKAL        endif
CKAL        i = index(cline,'=')
CKAL        read (cline(i+1:),*)   hminb,hmaxb
CKAL        if     (mnproc.eq.1) then  ! .b file from 1st tile only
CKAL          close(unit=uoff+9)
CKAL          write (lp,'(/(1x,a))') preambl,cline
CKAL        endif
CKALc
CKAL        call zaiopf(trim(flnmshlf)//'.a','old', 9)
CKAL        call zaiord(util3,ip,.false., hmina,hmaxa, 9)
CKAL        call zaiocl(9)
CKALc
CKAL        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
CKAL     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
CKAL          if     (mnproc.eq.1) then
CKAL          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
CKAL     &      'error - .a and .b files not consistent:',
CKAL     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
CKAL     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
CKAL          endif
CKAL          call xcstop('(geopar)')
CKAL                 stop '(geopar)'
CKAL        endif
CKALc
CKAL!$OMP   PARALLEL DO PRIVATE(j,i)
CKAL!$OMP&           SCHEDULE(STATIC,jblk)
CKAL        do j= 1,jj
CKAL          do i= 1,ii
CKAL            if     (ip(i,j).eq.0) then
CKAL              util3(i,j) = 0.0  !land
CKAL            elseif (util3(i,j).gt.0.5*hugel) then
CKAL              util3(i,j) = 0.0  !ice shelf over ocean
CKAL            elseif (util3(i,j).le.0.0) then
CKAL              util3(i,j) = 0.0  !ice shelf over ocean
CKAL            else
CKAL              util3(i,j) = 1.0  !open ocean
CKAL            endif
CKAL          enddo
CKAL        enddo
CKAL        call xctilr(util3,1,1, nbdy,nbdy, halo_ps)
CKAL        ishlf(:,:) = 0  !for jj:jdm and ii:idm
CKAL!$OMP   PARALLEL DO PRIVATE(j,i)
CKAL!$OMP&           SCHEDULE(STATIC,jblk)
CKAL        do j= 1-nbdy,jj+nbdy
CKAL          do i= 1-nbdy,ii+nbdy
CKAL            ishlf(i,j) = util3(i,j)
CKAL            util2(i,j) = ip(i,j)
CKAL          enddo
CKAL        enddo
CKALc
CKAL        call xcsum(sum_is,  util3,ip)
CKAL        call xcsum(sum_ip,  util2,ip)
CKAL        call xcsum(sum_isa, scp2, ishlf)
CKAL        call xcsum(area,    scp2, ip)
CKAL        nishlf = nint(sum_ip) - nint(sum_is)
CKAL        if     (mnproc.eq.1) then
CKAL        write (lp,'(/a,i9,f10.2)')
CKAL     &         ' number of ice shelf points and area (10^6 km^2):',
CKAL     &         nishlf,(area-sum_isa)*1.d-12
CKAL        endif
CKAL        call xcsync(flush_lp)
CKAL      endif !ishelf
CKALc
CKALc --- In arctic (tripole) domain, top row of mass points is redundent,
CKALc ---  so always use ipa, based on ishlf, for mass sums
CKAL#if defined(ARCTIC)
CKAL      ipa(:,:) = ishlf(:,:)
CKAL      if     (jj+j0.eq.jtdm) then
CKALc ---   mask top row of mass points
CKAL        ipa(:,jj:jj+nbdy) = 0
CKAL      endif
CKAL#else
CKALc --- Not a tripole domain, so ipa=ishlf
CKAL      ipa(:,:) = ishlf(:,:)
CKAL#endif
CKALc
CKAL      call xcsum(avgbot, util1,ipa)
CKAL      call xcsum(area,   scp2, ipa)
CKAL      avgbot=avgbot/area
CKAL      if     (mnproc.eq.1) then
CKAL      write (lp,'(/a,f9.1,f10.2)')
CKAL     &       ' mean basin depth (m) and area (10^6 km^2):',
CKAL     &       avgbot,area*1.e-12
CKAL      endif
CKAL      call xcsync(flush_lp)
CKALc
CKALc --- calculate dp0k and ds0k?
CKAL      if     (dp00.lt.0.0) then
CKALc ---   dp0k and ds0k already input
CKAL        dp00 =onem*dp0k(1)
CKAL        dp00x=onem*dp0k(kk-1)
CKAL        dp00i=onem*dp00i
CKAL        dpms = 0.0
CKAL        do k=1,kk
CKAL          dpm     = dp0k(k)
CKAL          dpms    = dpms + dpm
CKAL          dp0k(k) = dp0k(k)*onem
CKAL          if     (mnproc.eq.1) then
CKAL          write(lp,135) k,dp0k(k)*qonem,dpm,dpms
CKAL          endif
CKAL          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
CKAL            write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
CKAL          endif
CKAL          call xcsync(flush_lp)
CKAL        enddo !k
CKAL        dsms = 0.0
CKAL        do k=1,nsigma
CKAL          dsm     = ds0k(k)
CKAL          dsms    = dsms + dsm
CKAL          ds0k(k) = ds0k(k)*onem
CKAL          if     (mnproc.eq.1) then
CKAL          write(lp,130) k,ds0k(k)*qonem,dsm,dsms
CKAL          endif
CKAL          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
CKAL            write(6,*) 'geopar: ds0k   = ',ds0k(k),k,mnproc
CKAL          endif
CKAL          call xcsync(flush_lp)
CKAL        enddo !k
CKAL        if     (mnproc.eq.1) then
CKAL        write(lp,*)
CKAL        endif
CKAL      else
CKALc ---   calculate dp0k and ds0k
CKALc
CKALc ---   logorithmic k-dependence of dp0 (deep z's)
CKAL        dp00 =onem*dp00
CKAL        dp00x=onem*dp00x
CKAL        dp00i=onem*dp00i
CKAL        if     (isopyc) then
CKAL          dp0k(1)=thkmin*onem
CKAL        else
CKAL          dp0k(1)=dp00
CKAL        endif
CKAL        dpm  = dp0k(1)*qonem
CKAL        dpms = dpm
CKAL        if     (mnproc.eq.1) then
CKAL        write(lp,*)
CKAL        write(lp,135) 1,dp0k(1)*qonem,dpm,dpms
CKAL        endif
CKAL 135    format('dp0k(',i2,') =',f7.2,' m',
CKAL     &            '    thkns =',f7.2,' m',
CKAL     &            '    depth =',f8.2,' m')
CKAL        call xcsync(flush_lp)
CKALc
CKAL        dp0kf=1.0
CKAL        do k=2,kk
CKAL          dp0kf=dp0kf*dp00f
CKAL          if     (k.le.nhybrd) then
CKAL            if     (dp00f.ge.1.0) then
CKAL              dp0k(k)=min(dp00*dp0kf,dp00x)
CKAL            else
CKAL              dp0k(k)=max(dp00*dp0kf,dp00x)
CKAL            endif
CKAL          else
CKAL            dp0k(k)=0.0
CKAL          endif
CKAL          dpm  = dp0k(k)*qonem
CKAL          dpms = dpms + dpm
CKAL          if     (mnproc.eq.1) then
CKAL          write(lp,135) k,dp0k(k)*qonem,dpm,dpms
CKAL          endif
CKAL          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
CKAL            write(6,*) 'geopar: dp0kf  = ',dp0kf,    mnproc
CKAL            write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
CKAL          endif
CKAL          call xcsync(flush_lp)
CKAL        enddo !k
CKALc
CKALc ---   logorithmic k-dependence of ds0 (shallow z-s)
CKAL        ds00 =onem*ds00
CKAL        ds00x=onem*ds00x
CKAL        if     (isopyc) then
CKAL          ds0k(1)=thkmin*onem
CKAL        else
CKAL          ds0k(1)=ds00
CKAL        endif
CKAL        dsm  = ds0k(1)*qonem
CKAL        dsms = dsm
CKAL        if     (mnproc.eq.1) then
CKAL        write(lp,*)
CKAL        write(lp,130) 1,ds0k(1)*qonem,dsm,dsms
CKAL        endif
CKAL 130    format('ds0k(',i2,') =',f7.2,' m',
CKAL     &            '    thkns =',f7.2,' m',
CKAL     &            '    depth =',f8.2,' m')
CKAL        call xcsync(flush_lp)
CKALc
CKAL        ds0kf=1.0
CKAL        do k=2,nsigma
CKAL          ds0kf=ds0kf*ds00f
CKAL          if     (ds00f.ge.1.0) then
CKAL            ds0k(k)=min(ds00*ds0kf,ds00x)
CKAL          else
CKAL            ds0k(k)=max(ds00*ds0kf,ds00x)
CKAL          endif
CKAL          dsm  = ds0k(k)*qonem
CKAL          dsms = dsms + dsm
CKAL          if     (mnproc.eq.1) then
CKAL          write(lp,130) k,ds0k(k)*qonem,dsm,dsms
CKAL          endif
CKAL          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
CKAL            write(6,*) 'geopar: ds0kf  = ',ds0kf,    mnproc
CKAL            write(6,*) 'geopar: ds0k   = ',ds0k(k),k,mnproc
CKAL          endif
CKAL          call xcsync(flush_lp)
CKAL        enddo !k
CKAL        if     (mnproc.eq.1) then
CKAL        write(lp,*)
CKAL        endif
CKAL      endif !input:calculate dp0k,ds0k
CKALc
CKALc --- start and stop depths for terrain following coordinate
CKAL      if     (nsigma.eq.0) then
CKAL        dpns    = dp0k(1)
CKAL        dsns    = 0.0
CKAL        ds0k(1) = dp0k(1)
CKAL        do k= 2,kk
CKAL          ds0k(k)=0.0
CKAL        enddo !k
CKAL      else
CKAL        dpns = 0.0
CKAL        dsns = 0.0
CKAL        do k=1,nsigma
CKAL          dpns = dpns + dp0k(k)
CKAL          dsns = dsns + ds0k(k)
CKAL        enddo !k
CKAL        do k= nsigma+1,kk
CKAL          ds0k(k)=0.0
CKAL        enddo !k
CKAL      endif !nsigma
CKAL      dpns = dpns*qonem  !depths is in m
CKAL      dsns = dsns*qonem  !depths is in m
CKALc
CKAL      if     (mnproc.eq.1) then
CKAL      write(lp,131) nsigma,dpns,dsns
CKAL      endif
CKAL 131  format('nsigma = ',i2,
CKAL     &       '    deep    =',f8.2,' m',
CKAL     &       '    shallow =',f8.2,' m' )
CKAL      call flush(lp)
CKALc
CKALc --- initialize thermobaric reference state arrays.
CKALc
CKAL      if     (kapref.eq.-1) then
CKAL        if     (mnproc.eq.1) then  ! .b file from 1st tile only
CKAL          write (lp,'(3a)') ' reading thermobaric reference file from ',
CKAL     &                           trim(flnmforw), 'tbaric.[ab]'
CKAL          open (unit=uoff+9,file=trim(flnmforw)//'tbaric.b',
CKAL     &          status='old')
CKAL          read (     uoff+9,'(a79)')  preambl
CKAL        endif
CKAL        call xcsync(flush_lp)
CKAL        call zagetc(cline,ios, uoff+9)
CKAL        if     (ios.ne.0) then
CKAL          if     (mnproc.eq.1) then
CKAL            write(lp,'(/ a,i4,i9 /)')
CKAL     &        'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
CKAL          endif !1st tile
CKAL          call xcstop('(geopar)')
CKAL                 stop '(geopar)'
CKAL        endif
CKAL        i = index(cline,'=')
CKAL        read (cline(i+1:),*)   hminb,hmaxb
CKAL        if     (mnproc.eq.1) then  ! .b file from 1st tile only
CKAL          close(unit=uoff+9)
CKAL          write (lp,'(/(1x,a))') preambl,cline
CKAL        endif
CKALc
CKALc ---   input field is between 1.0 and 3.0 and indicates the
CKALc ---   relative strength of the two nearest reference states,
CKALc ---     e.g. 1.7 is 70% ref2 and 30% ref1
CKALc ---     and  2.3 is 70% ref2 and 30% ref3.
CKALc
CKAL        call zaiopf(trim(flnmforw)//'tbaric.a','old', 9)
CKAL        call zaiord(util1,ip,.false., hmina,hmaxa, 9)
CKAL        call zaiocl(9)
CKALc
CKAL        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
CKAL     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
CKAL          if     (mnproc.eq.1) then
CKAL          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
CKAL     &      'error - .a and .b files not consistent:',
CKAL     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
CKAL     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
CKAL          endif
CKAL          call xcstop('(geopar)')
CKAL                 stop '(geopar)'
CKAL        endif
CKALc
CKAL        do j= 1,jj
CKAL          do i= 1,ii
CKAL            if     (ip(i,j).eq.0) then
CKAL              util1(i,j) = 1.0 !land
CKAL            endif
CKAL          enddo
CKAL        enddo
CKALc
CKAL        vland = 1.0
CKAL        call xctilr(util1,  1,1, nbdy,nbdy, halo_ps)
CKAL        vland = 0.0
CKALc
CKALc       kapi is the 2nd reference state (1st is always 2)
CKALc       skap is the scale factor (0.0-1.0) for the 1st reference state
CKALc
CKALc       assumes that reference states 1 and 3 are never next to each other.
CKALc
CKAL        do j= 1,jj
CKAL          do i= 1,ii
CKAL            if     (max(util1(i,  j),
CKAL     &                  util1(i-1,j),
CKAL     &                  util1(i+1,j),
CKAL     &                  util1(i,  j-1),
CKAL     &                  util1(i,  j+1) ).gt.2.0) then
CKAL              util2(i,j) = 3.0              !kapi
CKAL               skap(i,j) = 3.0 - util1(i,j)
CKAL            else
CKAL              util2(i,j) = 1.0              !kapi
CKAL               skap(i,j) = util1(i,j) - 1.0
CKAL            endif
CKAL          enddo
CKAL        enddo
CKAL        vland = 1.0
CKAL        call xctilr(util2, 1,1, nbdy,nbdy, halo_ps)
CKAL        call xctilr(skap,  1,1, nbdy,nbdy, halo_ps)
CKAL        vland = 0.0
CKALc
CKAL        kapi(:,:) = util2(:,:)
CKAL      else
CKAL        skap(:,:) = 1.0     !for diagnostics only
CKAL        kapi(:,:) = kapref  !for diagnostics only
CKAL      endif !kapref.eq.-1:else
CKALc
CKALc --- initialize some arrays
CKALc --- set depthu,dpu,utotn,pgfx,depthv,dpv,vtotn,pgfy to zero everywhere,
CKALc --- so that they can be used at "lateral neighbors" of u and v points.
CKALc --- similarly for pbot,dp at neighbors of q points.
CKALc
CKAL      disp_count=0
CKALc
CKAL!$OMP PARALLEL DO PRIVATE(j,i,k,ktr)
CKAL!$OMP&         SCHEDULE(STATIC,jblk)
CKAL      do j=1-nbdy,jj+nbdy
CKAL        do i=1-nbdy,ii+nbdy
CKAL          p(     i,j,1)=0.0
CKAL          pu(    i,j,1)=0.0
CKAL          pv(    i,j,1)=0.0
CKAL          utotn( i,j)=0.0
CKAL          vtotn( i,j)=0.0
CKAL          pgfx(  i,j)=0.0
CKAL          pgfy(  i,j)=0.0
CKAL          gradx( i,j)=0.0
CKAL          grady( i,j)=0.0
CKAL          depthu(i,j)=0.0
CKAL          depthv(i,j)=0.0
CKAL          pbot(  i,j)=0.0
CKALc
CKAL          displd_mn(i,j)=0.0
CKAL          dispqd_mn(i,j)=0.0
CKAL          tidepg_mn(i,j)=0.0
CKALc
CKAL          psikk( i,j,1)=0.0
CKAL          psikk( i,j,2)=0.0
CKAL          thkk(  i,j,1)=0.0
CKAL          thkk(  i,j,2)=0.0
CKALc
CKAL          ubavg( i,j,1)=hugel
CKAL          ubavg( i,j,2)=hugel
CKAL          ubavg( i,j,3)=hugel
CKAL          vbavg( i,j,1)=hugel
CKAL          vbavg( i,j,2)=hugel
CKAL          vbavg( i,j,3)=hugel
CKAL          utotm( i,j)=hugel
CKAL          vtotm( i,j)=hugel
CKAL          uflux( i,j)=hugel
CKAL          vflux( i,j)=hugel
CKAL          uflux1(i,j)=hugel
CKAL          vflux1(i,j)=hugel
CKAL          uflux2(i,j)=hugel
CKAL          vflux2(i,j)=hugel
CKAL          uflux3(i,j)=hugel
CKAL          vflux3(i,j)=hugel
CKAL          uja(   i,j)=hugel
CKAL          ujb(   i,j)=hugel
CKAL          via(   i,j)=hugel
CKAL          vib(   i,j)=hugel
CKAL          do k=1,kk
CKAL            dp( i,j,k,1)=0.0
CKAL            dp( i,j,k,2)=0.0
CKAL            dpu(i,j,k,1)=0.0
CKAL            dpu(i,j,k,2)=0.0
CKAL            dpv(i,j,k,1)=0.0
CKAL            dpv(i,j,k,2)=0.0
CKALc
CKAL            u(  i,j,k,1)=hugel
CKAL            u(  i,j,k,2)=hugel
CKAL            v(  i,j,k,1)=hugel
CKAL            v(  i,j,k,2)=hugel
CKALc
CKAL            uflx(  i,j,k)=hugel
CKAL            vflx(  i,j,k)=hugel
CKALc
CKAL            dpav(  i,j,k)=0.0
CKAL            uflxav(i,j,k)=0.0
CKAL            vflxav(i,j,k)=0.0
CKAL            diaflx(i,j,k)=0.0
CKALc
CKAL            do ktr= 1,ntracr
CKAL              tracer(i,j,k,1,ktr)=0.0
CKAL              tracer(i,j,k,2,ktr)=0.0
CKAL            enddo
CKAL          enddo
CKAL        enddo
CKAL      enddo
CKAL!$OMP END PARALLEL DO
CKALc
CKAL!$OMP PARALLEL DO PRIVATE(j,l,i,k)
CKAL!$OMP&         SCHEDULE(STATIC,jblk)
CKAL      do j=1,jj
CKAL        do l=1,isp(j) !ok
CKAL          do i=max(1,ifp(j,l)),min(ii,ilp(j,l)+1)
CKAL            ubavg(i,j,1)=0.0
CKAL            ubavg(i,j,2)=0.0
CKAL            ubavg(i,j,3)=0.0
CKAL            utotm (i,j)=0.0
CKAL            uflux (i,j)=0.0
CKAL            uflux2(i,j)=0.0
CKAL            uflux3(i,j)=0.0
CKAL            uja(i,j)=0.0
CKAL            ujb(i,j)=0.0
CKALc
CKAL            do k=1,kk
CKAL              uflx(i,j,k)=0.0
CKAL              u(i,j,k,1)=0.0
CKAL              u(i,j,k,2)=0.0
CKAL            enddo
CKAL          enddo
CKAL        enddo
CKAL      enddo
CKALc
CKAL      call xctilr(ubavg,    1,   3, nbdy,nbdy, halo_us)  ! note scalar
CKAL      call xctilr(utotm,    1,   1, nbdy,nbdy, halo_us)  ! note scalar
CKAL      call xctilr(uflux,    1,   1, nbdy,nbdy, halo_us)  ! note scalar
CKAL      call xctilr(uflux2,   1,   1, nbdy,nbdy, halo_us)  ! note scalar
CKAL      call xctilr(uflux3,   1,   1, nbdy,nbdy, halo_us)  ! note scalar
CKAL      call xctilr(uja,      1,   1, nbdy,nbdy, halo_us)
CKAL      call xctilr(ujb,      1,   1, nbdy,nbdy, halo_us)
CKAL      call xctilr(uflx,     1,  kk, nbdy,nbdy, halo_us)  ! note scalar
CKAL      call xctilr(u,        1,2*kk, nbdy,nbdy, halo_us)  ! note scalar
CKALc
CKAL!$OMP PARALLEL DO PRIVATE(i,l,j,k)
CKAL!$OMP&         SCHEDULE(STATIC)
CKAL      do i=1,ii
CKAL        do l=1,jsp(i) !ok
CKAL          do j=max(1,jfp(i,l)),min(jj,jlp(i,l)+1)
CKAL            vbavg(i,j,1)=0.0
CKAL            vbavg(i,j,2)=0.0
CKAL            vbavg(i,j,3)=0.0
CKAL            vtotm (i,j)=0.0
CKAL            vflux (i,j)=0.0
CKAL            vflux2(i,j)=0.0
CKAL            vflux3(i,j)=0.0
CKAL            via(i,j)=0.0
CKAL            vib(i,j)=0.0
CKALc
CKAL            do k=1,kk
CKAL              vflx(i,j,k)=0.0
CKAL              v(i,j,k,1)=0.0
CKAL              v(i,j,k,2)=0.0
CKAL            enddo
CKAL          enddo
CKAL        enddo
CKAL      enddo
CKALc
CKAL      call xctilr(vbavg,    1,   3, nbdy,nbdy, halo_vs)  ! note scalar
CKAL      call xctilr(vtotm,    1,   1, nbdy,nbdy, halo_vs)  ! note scalar
CKAL      call xctilr(vflux,    1,   1, nbdy,nbdy, halo_vs)  ! note scalar
CKAL      call xctilr(vflux2,   1,   1, nbdy,nbdy, halo_vs)  ! note scalar
CKAL      call xctilr(vflux3,   1,   1, nbdy,nbdy, halo_vs)  ! note scalar
CKAL      call xctilr(via,      1,   1, nbdy,nbdy, halo_vs)  ! note scalar
CKAL      call xctilr(vib,      1,   1, nbdy,nbdy, halo_vs)  ! note scalar
CKAL      call xctilr(vflx,     1,  kk, nbdy,nbdy, halo_vs)  ! note scalar
CKAL      call xctilr(v,        1,2*kk, nbdy,nbdy, halo_vs)  ! note scalar
CKALc
      return
      end
c
c
c> Revision history:
c>
c> May  1997 - extended list of variables set to 'hugel' on land
c> Oct. 1999 - added code that defines the vertical distribution of dp0
c>             used in hybgen
c> Jan. 2000 - added mapflg logic for different projections
c> Feb. 2000 - added dp00f for logorithmic z-level spacing
c> Mar. 2000 - added dp00s for sigma-spacing in shallow water
c> May  2000 - conversion to SI units (still wrong corio)
c> Feb. 2001 - removed rotated grid option
c> Jan. 2002 - more flexible Z-sigma-Z vertical configuration
c> Jan. 2002 - all grids now via array input
c> Sep. 2004 - define kapi and skap for thermobaricity
c> Oct. 2008 - dp0k and ds0k can now be input, see blkdat.F
c> Mar. 2012 - replaced dssk with dpns and dsns
c> Apr. 2014 - added ishlf
c> Apr. 2014 - added ipa
c> Feb. 2015 - added pang for coupled cases
