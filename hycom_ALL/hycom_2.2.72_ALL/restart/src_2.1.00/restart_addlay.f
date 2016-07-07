      module mod_restart_addlay
      contains

      subroutine restart_in(nstepx, dtimex, layadd)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
c
      integer nstepx,layadd
      real*8  dtimex
c
c     read in a restart file on unit 11.
c     leave layers 1:layadd empty (layer 1 in layer layadd+1).
c
      integer   i,i1,ios,j,k,k1,kkk
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
        write(lp,'(/ a /)') 'error - input is not a restart file'
        endif !1st tile
        call xcstop('(restart_in)')
               stop '(restart_in)'
      endif
      read( 11,'(a)') cline
      if     (mnproc.eq.1) then
      write(lp,'(a)') cline(1:len_trim(cline))
      call flush(lp)
      endif !1st tile
      i = index(cline,'=')
      read(cline(i+1:),*) nstepx,dtimex
c
      i1  = 1-nbdy
      k1  = layadd+1
      kkk = kdm-layadd
c
      call restart_in3d(u(   i1,i1,k1,1),kkk, iu, 'u       ')
      call restart_in3d(u(   i1,i1,k1,2),kkk, iu, 'u       ')
      call restart_in3d(v(   i1,i1,k1,1),kkk, iv, 'v       ')
      call restart_in3d(v(   i1,i1,k1,2),kkk, iv, 'v       ')
      call restart_in3d(dp(  i1,i1,k1,1),kkk, ip, 'dp      ')
      call restart_in3d(dp(  i1,i1,k1,2),kkk, ip, 'dp      ')
      call restart_in3d(temp(i1,i1,k1,1),kkk, ip, 'temp    ')
      call restart_in3d(temp(i1,i1,k1,2),kkk, ip, 'temp    ')
      call restart_in3d(saln(i1,i1,k1,1),kkk, ip, 'saln    ')
      call restart_in3d(saln(i1,i1,k1,2),kkk, ip, 'saln    ')
      call restart_in3d(th3d(i1,i1,k1,1),kkk, ip, 'th3d    ')
      call restart_in3d(th3d(i1,i1,k1,2),kkk, ip, 'th3d    ')
c
      call restart_in3d(ubavg,     3, iu, 'ubavg   ')
      call restart_in3d(vbavg,     3, iv, 'vbavg   ')
      call restart_in3d(pbavg,     3, ip, 'pbavg   ')
      call restart_in3d(pbot,      1, ip, 'pbot    ')
      call restart_in3d(psikk,     1, ip, 'psikk   ')
      call restart_in3d(thkk,      1, ip, 'thkk    ')
      call restart_in3d(dpmixl,    2, ip, 'dpmixl  ')
c
      if (icegln) then
        read ( 11,'(a)',iostat=ios) cline
        if     (ios.ne.0 .or. cline(1:8).ne.'temice  ') then
c
c ---     assume this is an addition of ice to the simulation.
c ---     only works with no tracers.
c
          if     (mnproc.eq.1) then
          write(lp,'(/ a /)') 'adding ice to the simulation.'
          call flush(lp)
          endif !1st tile
c
          do j= 1-nbdy,jdm+nbdy
            do i= 1-nbdy,idm+nbdy
              temice(i,j) = temp(i,j,k1,1)
              covice(i,j) = 0.0
              thkice(i,j) = 0.0
            enddo
          enddo
        else
c
c ---     reposition file for ice input
c
          rewind(11)
          do k= 1,12*kkk+16
            read (11,'(a)') cline
*           if     (mnproc.eq.1) then
*           write(lp,'(a)') cline
*           endif !1st tile
          enddo
          call restart_in3d(temice,    1, ip, 'temice  ')
          call restart_in3d(covice,    1, ip, 'covice  ')
          call restart_in3d(thkice,    1, ip, 'thkice  ')
        endif
      endif
      if (trcout) then
        call restart_in3d(tracer(i1,i1,k1),kkk, ip, 'tracer  ')
      endif
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

      end module mod_restart_addlay

      program restart_addlay
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      use mod_restart_addlay
      implicit none
c
      include 'common_blocks.h'
c
c     add zero thickness layers between layers 1 and 2 of a HYCOM restart.
c
c     used when converting a MICOM-like simulation (with layer 1 as the
c     mixed layer) to a typical HYCOM simulation with z-levels inside
c     the mixed layer.  the zero thickness layers will become z-levels
c     via hybgen remapping inside HYCOM.
c
      integer k,ios,l,layadd,ni,nstep
      real    astep(2),time0
      real*8  dtime
c
      call xcspmd
      call zaiost
c
      layadd = -1  ! enable layadd option
      call blkdat(.true.,layadd)
      trcout=.false.
c
      call geopar
c
      call restart_in( nstep, dtime, layadd)
c
      if     (layadd.gt.0) then
c
c ---   add zero thickness layers betwen layer 1 and 2
c
        l = layadd+1  ! layer 1 was input in layer l
        dp(:,:,1,:) = dp(:,:,l,:)
        do k= 1,layadd
             u(:,:,k,  :) =    u(:,:,l,:)
             v(:,:,k,  :) =    v(:,:,l,:)
          temp(:,:,k,  :) = temp(:,:,l,:)
          saln(:,:,k,  :) = saln(:,:,l,:)
          th3d(:,:,k,  :) = th3d(:,:,l,:)
            dp(:,:,k+1,:) = 0.0
        enddo
      endif
c
      call restart_out(nstep, dtime)

      end program restart_addlay
