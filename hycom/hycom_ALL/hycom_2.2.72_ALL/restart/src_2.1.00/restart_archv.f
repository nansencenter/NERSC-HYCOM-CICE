      module mod_restart_archv
      contains

      subroutine archv_in(nstepx, dtimex)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
c
      integer nstepx
      real*8  dtimex
c
c     read in an archive file on unit 11.
c
      integer   i1,k
      character cline*80
c
      write(lp,*) 'flnmrsa = ',flnmrsa(1:len_trim(flnmrsa))
      write(lp,*) 'file=',flnmrsa(1:len_trim(flnmrsa))//'.b'
      call flush(lp)
c
      open (unit=11,file=flnmrsa(1:len_trim(flnmrsa))//'.b',
     &      status='old',action='read',form='formatted')
      call zaiopf(flnmrsa(1:len_trim(flnmrsa))//'.a','old', 11)
      do k= 1,10
        read( 11,'(a)') cline
        if     (mnproc.eq.1) then
        write(lp,'(a)') cline(1:len_trim(cline))
        call flush(lp)
        endif !1st tile
      enddo
c
      if     (cline(24:28).ne.'model') then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') 'error - input file is not a HYCOM archv'
        endif !1st tile
        call xcstop('(archv_in)')
               stop '(archv_in)'
      endif
c
      i1=1-nbdy
c
c --- srfhgt=montg1+thref*pbavg
      call archv_in2d(pbavg( i1,i1,1), ip, 'montg1  ')
      call archv_in2d(pbavg( i1,i1,2), ip, 'srfhgt  ')
      pbavg(:,:,1) = (pbavg(:,:,2) - pbavg(:,:,1))*1.e3
c
      call zaiosk(11)
      read( 11,'(a)') cline              ! 'surflx  '
      read(cline(11:33),*) nstepx, dtimex
c
      call archv_sk2d(                     'salflx  ')
      call archv_sk2d(                     'bl_dpth ')
c
      call archv_in2d(dpmixl(i1,i1,1), ip, 'mix_dpth')
c
      call archv_sk2d(                     'tmix    ')
      call archv_sk2d(                     'smix    ')
      call archv_sk2d(                     'thmix   ')
      call archv_sk2d(                     'umix    ')
      call archv_sk2d(                     'vmix    ')
      if (icegln) then
        call archv_in2d(covice(i1,i1),   ip, 'covice  ')
        call archv_in2d(thkice(i1,i1),   ip, 'thkice  ')
        call archv_in2d(temice(i1,i1),   ip, 'temice  ')
      else
        covice(:,:) = 0.0
        thkice(:,:) = 0.0
        temice(:,:) = 0.0
      endif
c
      call archv_in2d(ubavg( i1,i1,1), iu, 'u_btrop ')
      call archv_in2d(vbavg( i1,i1,1), iv, 'v_btrop ')
c
      do k= 1,kdm
        call archv_in2d(u(   i1,i1,k,1), iu, 'u-vel.  ')
        call archv_in2d(v(   i1,i1,k,1), iu, 'v-vel.  ')
        call archv_in2d(dp(  i1,i1,k,1), ip, 'thknss  ')
        call archv_in2d(temp(i1,i1,k,1), ip, 'temp    ')
        call archv_in2d(saln(i1,i1,k,1), ip, 'salin   ')
        call archv_in2d(th3d(i1,i1,k,1), ip, 'density ')
      enddo
c
      close (unit=11)
      call zaiocl(11)
      return
      end subroutine archv_in

      subroutine archv_in2d(field, mask, cfield)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & field
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & mask
      character cfield*8
c
c --- read a single archive 2-d array field from unit 11.
c
      integer   i,j,layer,nnstep
      real      hmina,hminb,hmaxa,hmaxb,thet,ttime
      character cline*80
c
      real       hspval
      parameter (hspval=0.5*2.0**100)
c
*     if     (mnproc.eq.1) then
*     write(lp,'(a,i3,2x,a)') 'archv_in2d - cfield = ',cfield
*     call flush(lp)
*     endif !1st tile
      call zaiord(field, mask,.false., hmina,hmaxa, 11)
c
      do j= 1-nbdy,jdm+nbdy
        do i= 1-nbdy,idm+nbdy
          if     (field(i,j).gt.hspval) then
            field(i,j) = 0.0
          endif
        enddo
      enddo
c
      read ( 11,'(a)')  cline
*     if     (mnproc.eq.1) then
*     write (lp,'(a)')  cline(1:len_trim(cline))
*     endif !1st tile
      if     (cline(1:8).ne.cfield) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,a /)') cline(1:len_trim(cline)),
     &         'error in archv_in2d - expected ',cfield
        endif !1st tile
        call xcstop('(archv_in2d)')
               stop '(archv_in2d)'
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nnstep,ttime,layer,thet,hminb,hmaxb
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,i3 / a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    'iunit = ',11,
     &    cline,
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        endif !1st tile
        call xcstop('(archv_in2d)')
               stop '(archv_in2d)'
      endif
c
      return
      end subroutine archv_in2d

      subroutine archv_sk2d(cfield)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      character cfield*8
c
c --- skip a single archive 2-d array field from unit 11.
c
      character cline*80
c
*     if     (mnproc.eq.1) then
*     write(lp,'(a,i3,2x,a)') 'archv_sk2d - cfield = ',cfield
*     call flush(lp)
*     endif !1st tile
      call zaiosk(11)
c
      read ( 11,'(a)')  cline
*     if     (mnproc.eq.1) then
*     write (lp,'(a)')  cline(1:len_trim(cline))
*     endif !1st tile
      if     (cline(1:8).ne.cfield) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,a /)') cline(1:len_trim(cline)),
     &         'error in archv_sk2d - expected ',cfield
        endif !1st tile
        call xcstop('(archv_sk2d)')
               stop '(archv_sk2d)'
      endif
c
      return
      end subroutine archv_sk2d

      subroutine restart_in
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
c
c     read in a small part of a a restart file on unit 11.
c
      character cline*80
c
      open (unit=11,file=flnmrsi(1:len_trim(flnmrsi))//'.b',
     &      status='old',action='read',form='formatted')
      call zaiopf(flnmrsi(1:len_trim(flnmrsi))//'.a','old', 11)
      read( 11,'(a)') cline
      if     (mnproc.eq.1) then
      write(lp,'(a)') cline(1:len_trim(cline))
      endif !1st tile
      if     (cline(1:9).ne.'RESTART:') then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') 'error - input file is not a HYCOM restart'
        endif !1st tile
        call xcstop('(restart_in)')
               stop '(restart_in)'
      endif
      read( 11,'(a)') cline
      if     (mnproc.eq.1) then
      write(lp,'(a)') cline(1:len_trim(cline))
      call flush(lp)
      endif !1st tile
c
      call restart_sk3d(2*kdm, 'u       ')
      call restart_sk3d(2*kdm, 'v       ')
      call restart_sk3d(2*kdm, 'dp      ')
      call restart_sk3d(2*kdm, 'temp    ')
      call restart_sk3d(2*kdm, 'saln    ')
      call restart_sk3d(2*kdm, 'th3d    ')
c
      call restart_sk3d(3,     'ubavg   ')
      call restart_sk3d(3,     'vbavg   ')
      call restart_sk3d(3,     'pbavg   ')
c
      call restart_in3d(pbot,  1, ip, 'pbot    ')
      call restart_in3d(psikk, 1, ip, 'psikk   ')
      call restart_in3d(thkk,  1, ip, 'thkk    ')
c
      close (unit=11)
      call zaiocl(11)
      return
      end subroutine restart_in

      subroutine restart_in3d(field,l, mask, cfield)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      integer   l
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l) ::
     & field
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & mask
      character cfield*8
c
c --- read a single restart 3-d array field from unit 11.
c
      integer   i,layer,level,k
      real      hmina(l),hminb,hmaxa(l),hmaxb
      character cline*80
c
*     if     (mnproc.eq.1) then
*     write(lp,'(a,i3,2x,a)') 'restart_in3d - l,cfield = ',l,cfield
*     call flush(lp)
*     endif !1st tile
      call zaiord3(field,l, mask,.false., hmina,hmaxa, 11)
c
      do k= 1,l
        read ( 11,'(a)')  cline
*       if     (mnproc.eq.1) then
*       write (lp,'(a)')  cline(1:len_trim(cline))
*       endif !1st tile
        if     (cline(1:8).ne.cfield) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cline(1:len_trim(cline)),
     &           'error in restart_in3d - expected ',cfield
          endif !1st tile
          call xcstop('(restart_in3d)')
                 stop '(restart_in3d)'
        endif
        i = index(cline,'=')
        read (cline(i+1:),*) layer,level,hminb,hmaxb
        if     (abs(hmina(k)-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa(k)-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,3i3 / a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      'iunit,k,l = ',11,k,l,
     &      cline,
     &      '.a,.b min = ',hmina(k),hminb,hmina(k)-hminb,
     &      '.a,.b max = ',hmaxa(k),hmaxb,hmaxa(k)-hmaxb
          endif !1st tile
          call xcstop('(restart_in3d)')
                 stop '(restart_in3d)'
        endif
      enddo
c
      return
      end subroutine restart_in3d

      subroutine restart_sk3d(l, cfield)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      integer   l
      character cfield*8
c
c --- skip a single restart 3-d array field from unit 11.
c
      integer   k
      character cline*80
c
*     if     (mnproc.eq.1) then
*     write(lp,'(a,i3,2x,a)') 'restart_sk3d - l,cfield = ',l,cfield
*     call flush(lp)
*     endif !1st tile
c
      do k= 1,l
        call zaiosk(11)
        read ( 11,'(a)')  cline
*       if     (mnproc.eq.1) then
*       write (lp,'(a)')  cline(1:len_trim(cline))
*       endif !1st tile
        if     (cline(1:8).ne.cfield) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cline(1:len_trim(cline)),
     &           'error in restart_sk3d - expected ',cfield
          endif !1st tile
          call xcstop('(restart_in3d)')
                 stop '(restart_in3d)'
        endif
      enddo
c
      return
      end subroutine restart_sk3d

      subroutine restart_out(nstepx, dtimex)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
c
      integer nstepx
      real*8  dtimex
c
c     write out in a restart file on unit 12
c
      logical   lopen
      integer   i,j,k,l
      real      xmin(2*kdm),xmax(2*kdm)
      character cline*80
c
      call zaiopf(flnmrso(1:len_trim(flnmrso))//'.a','new', 12)
      if     (mnproc.eq.1) then
        open (unit=12,file=flnmrso(1:len_trim(flnmrso))//'.b',
     &        status='new',action='write',form='formatted')
        write(lp,'(a)') ' creating a new restart file'
      endif !1st tile
c
      if     (mnproc.eq.1) then
      write(12,'(a,3i6)') 'RESTART: iexpt,iversn,yrflag = ',
     &                              iexpt,iversn,yrflag
      write(cline,*)                nstepx,dtimex
      write(12,'(a,a)')   'RESTART: nstep,dtime = ',
     &                              cline(1:len_trim(cline))
      endif !1st tile
c
      call zaiowr3(u,      2*kdm, iu,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'u       ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(v,      2*kdm, iv,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'v       ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(dp,     2*kdm, ip,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'dp      ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(temp,   2*kdm, ip,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'temp    ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(saln,   2*kdm, ip,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'saln    ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(th3d,   2*kdm, ip,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'th3d    ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(ubavg,      3, iu,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 1,3
        do k= 0,0
          write(12,4100) 'ubavg   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(vbavg,      3, iv,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 1,3
        do k= 0,0
          write(12,4100) 'vbavg   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(pbavg,      3, ip,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 1,3
        do k= 0,0
          write(12,4100) 'pbavg   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(pbot,       1, ip,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 1,1
        do k= 0,0
          write(12,4100) 'pbot    ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(psikk,      1, ip,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 1,1
        do k= 0,0
          write(12,4100) 'psikk   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(thkk,       1, ip,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 1,1
        do k= 0,0
          write(12,4100) 'thkk    ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      endif !1st tile
      call zaiowr3(dpmixl,     2, ip,.false., xmin,xmax, 12, .true.)
      if     (mnproc.eq.1) then
      do l= 1,2
        do k= 0,0
          write(12,4100) 'dpmixl  ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      endif !1st tile
      if (icegln) then
        call zaiowr3(temice,     1, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(12,4100) 'temice  ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        endif !1st tile
        call zaiowr3(covice,     1, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(12,4100) 'covice  ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        endif !1st tile
        call zaiowr3(thkice,     1, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(12,4100) 'thkice  ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        endif !1st tile
      endif
      if (trcout) then
        call zaiowr3(tracer, kdm, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 0,0
          do k= 1,kdm
            write(12,4100) 'tracer  ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
          enddo
        enddo
        endif !1st tile
      endif
      call zaiofl(12)
      if     (mnproc.eq.1) then
      write(lp,'(a,f11.2)') ' restart created at model day',dtimex
      endif !1st tile
      call xcsync(flush_lp)
c
      return
 4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
      end subroutine restart_out

      end module mod_restart_archv

      program restart_archv
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      use mod_restart_archv
      implicit none
c
      include 'common_blocks.h'
c
c     create a HYCOM 2.0 restart file from an archive file.
c
      integer k,ios,l,layadd,ni,nstep
      real    astep(2),time0
      real*8  dtime
c
      call xcspmd
      call zaiost
c
      layadd = 0  ! disable layadd option
      call blkdat(.true.,layadd)
      trcout=.false.
c
      call geopar
c
      call restart_in
      call archv_in(nstep, dtime)
c
         u(:,:,:,2) =    u(:,:,:,1)
         v(:,:,:,2) =    v(:,:,:,1)
      temp(:,:,:,2) = temp(:,:,:,1)
      saln(:,:,:,2) = saln(:,:,:,1)
      th3d(:,:,:,2) = th3d(:,:,:,1)
        dp(:,:,:,2) =   dp(:,:,:,1)
c
      dpmixl(:,:,2) = dpmixl(:,:,1)
c
       ubavg(:,:,2) =  ubavg(:,:,1)
       ubavg(:,:,3) =  ubavg(:,:,1)
       vbavg(:,:,2) =  vbavg(:,:,1)
       vbavg(:,:,3) =  vbavg(:,:,1)
       pbavg(:,:,2) =  pbavg(:,:,1)
       pbavg(:,:,3) =  pbavg(:,:,1)
c
      call restart_out(nstep, dtime)

      end program restart_archv
