      module mod_restart_hycom
      contains

      subroutine restart_out(nstepx, dtimex)
      use mod_xc   ! HYCOM communication interface
      use mod_za   ! HYCOM I/O interface
      use mod_rest ! HYCOM restart arrays
      implicit none
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
        call zaiowr3(tracer, 2*kdm, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 0,1
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

      end module mod_restart_hycom

      program restart
      use mod_xc   ! HYCOM communication interface
      use mod_za   ! HYCOM I/O interface
      use mod_rest ! HYCOM restart arrays
      use mod_restart_hycom
      implicit none
c
c     convert restart from HYCOM 1.0 to HYCOM 2.0 format.
c
      integer ios,layadd,ni,nstep
      real*8  dtime
c
c --- HYCOM 1.0: real is real*8.
c
      real*8    astep(2),time0
      real*8,   allocatable, dimension (:,:,:,:) :: aa4
      real*8,   allocatable, dimension (:,:,:)   :: aa3
c
      call xcspmd
      call zaiost
c
      layadd = 0  ! disable layadd option
      call blkdat(.true.,layadd)
      trcout=.false.
c
      call rest_alloc
c
      allocate( aa4(idm,jdm,kdm,2) )
      allocate( aa3(idm,jdm,9)     )
      call geopar
c
      ni=11
      open (unit=ni,file=flnmrsi,status='old',form='unformatted')
      read (ni) astep,time0,aa4
      u    = aa4
      read (ni) astep,time0,aa4
      v    = aa4
      read (ni) astep,time0,aa4
      dp   = aa4
      read (ni) astep,time0,aa4
      temp = aa4
      read (ni) astep,time0,aa4
      saln = aa4
      read (ni) astep,time0,ss4
      th3d = aa4
      read (ni) astep,time0,aa3
      ubavg(:,:,:) = aa3(:,:,1:3)
      vbavg(:,:,:) = aa3(:,:,4:6)
      pbavg(:,:,:) = aa3(:,:,7:9)
      read (ni) astep,time0,aa3(:,:,1:5)
        pbot(:,:)   = aa3(:,:,1)
       psikk(:,:)   = aa3(:,:,2)
       thkkt(:,:)   = aa3(:,:,3)
      dpmixl(:,:,:) = aa3(:,:,4:5)
      if (icegln) then
        read (ni,iostat=ios) astep,time0,aa3(:,:,1:3)
        temice(:,:) = aa3(:,:,1)
        covice(:,:) = aa3(:,:,2)
        thkice(:,:) = aa3(:,:,3)
      endif
      if (trcout) then
        read (ni) astep,time0,aa4(:,:,:,1)
        tracer(:,:,:,1) = aa4(:,:,:,1)
        tracer(:,:,:,2) = aa4(:,:,:,1)
      endif
      close (unit=ni)
      nstep=nint(astep(1))*1048576 + nint(astep(2))
      dtime=nstep/(86400.0d0/baclin)

c
      call restart_out(nstep, dtime)

      end program restart
