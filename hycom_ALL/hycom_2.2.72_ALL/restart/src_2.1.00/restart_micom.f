      module mod_restart_micom
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

      subroutine rotate(a, s)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      real a(idm,jdm),s
c
c --- rotate and scale from MICOM to HYCOM orientation.
c
      real    b(idm,jdm)
      save    b
c
*     write(lp,*) 'call rotate, s = ',s
*     call xcsync(flush_lp)
c
      b = a
      if     (s.gt.0.0) then
        call rotate_u(a,b,s)  ! u or p grid
      else
        call rotate_v(a,b,s)  ! v grid
      endif
      return
      end subroutine rotate

      subroutine rotate_u(a,b, s)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      real a(idm,jdm),b(jdm,idm),s
c
c --- rotate and scale from MICOM to HYCOM orientation.
c --- u or p grid.
c
      real       hhuge
      real*4     spval
      parameter (hhuge=0.5e33,spval=2.0**100)
c
      integer i,j
c
      do j= 1,jdm-1
        do i= 1,idm
          if     (b(jdm-j,i).lt.hhuge) then
            a(i,j) = s*b(jdm-j,i)
          else
            a(i,j) = spval
          endif
        enddo
      enddo
      do i= 1,idm
        a(i,jdm) = spval
      enddo
      return
      end subroutine rotate_u

      subroutine rotate_v(a,b, s)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      real a(idm,jdm),b(jdm,idm),s
c
c --- rotate and scale from MICOM to HYCOM orientation.
c --- v grid.
c
      real       hhuge
      real*4     spval
      parameter (hhuge=0.5e33,spval=2.0**100)
c
      integer i,j
c
      do j= 1,jdm
        do i= 1,idm
          if     (b(jdm+1-j,i).lt.hhuge) then
            a(i,j) = s*b(jdm+1-j,i)
          else
            a(i,j) = spval
          endif
        enddo
      enddo
      return
      end subroutine rotate_v

      end module mod_restart_micom

      program restart_micom
      use mod_xc   ! HYCOM communication interface
      use mod_za   ! HYCOM I/O interface
      use mod_rest ! HYCOM restart arrays
      use mod_restart_micom
      implicit none
c
c     convert restart from MICOM to HYCOM 2.0 format.
c
      integer i,ios,j,k,layadd,n,ni
      real*8  dtime
c
c --- MICOM: integer is integer*4 and real is real*8.
c
      integer*4 nstep
      real*8    time0
      real*8,   allocatable, dimension (:,:,:,:) :: aa4
      real*8,   allocatable, dimension (:,:,:)   :: aa3
c
      call xcspmd
      call zaiost
c
      layadd = 0  ! disable layadd optiuon
      call blkdat(.false.,layadd)
      trcout=.false.
c
      call rest_alloc
c
      allocate( aa4(idm,jdm,kdm,2) )
      allocate( aa3(idm,jdm,11)    )
c
      call geopar
c
      ni=11
      open (unit=ni,file=flnmrsi,status='old',form='unformatted')
      read (ni) nstep,time0,aa4
      v = aa4                   ! note switch of u and v
      write(lp,*) 'input: v'
      call xcsync(flush_lp)
      read (ni) nstep,time0,aa4
      u = aa4                   ! note switch of u and v
      write(lp,*) 'input: u'
      call xcsync(flush_lp)
      read (ni) nstep,time0,aa4
      dp = aa4
      write(lp,*) 'input: dp'
      call xcsync(flush_lp)
      read (ni) nstep,time0,aa4
      temp = aa4
      write(lp,*) 'input: temp'
      call xcsync(flush_lp)
      read (ni) nstep,time0,aa4
      saln = aa4
      write(lp,*) 'input: saln'
      call xcsync(flush_lp)
      read (ni) nstep,time0,aa3
      thmix(:,:,:) = aa3(:,:,1:2)
      vbavg(:,:,:) = aa3(:,:,3:5)   ! note switch of u and v
      ubavg(:,:,:) = aa3(:,:,6:8)   ! note switch of u and v
      pbavg(:,:,:) = aa3(:,:,9:11)
      write(lp,*) 'input: thmix,vbavg,ubavg,pbavg'
      call xcsync(flush_lp)
      read (ni) nstep,time0,aa3(:,:,1:3)
       pbot(:,:) = aa3(:,:,1)
      psikk(:,:) = aa3(:,:,2)
       thkk(:,:) = aa3(:,:,3)
      write(lp,*) 'input: pbot,psikk,thkk'
      call xcsync(flush_lp)
      if (trcout) then
        read (ni) nstep,time0,aa4(:,:,:,1)
        tracer(:,:,:,1) = aa4(:,:,:,1)
        tracer(:,:,:,2) = aa4(:,:,:,1)
      endif
      close (unit=ni)
      dtime=nstep/(86400.0d0/baclin)
c
      do k= 1,kdm
        do n= 1,2
          call rotate(   u(1,1,k,n),  0.01)  ! micom v
          call rotate(   v(1,1,k,n), -0.01)  ! micom u
          call rotate(  dp(1,1,k,n),  0.1)
          call rotate(temp(1,1,k,n),  1.0)
          call rotate(saln(1,1,k,n),  1.0)
          th3d(:,:,k,n) = sigma(k)-thbase
        enddo
        if (trcout) then
          call rotate(tracer(1,1,k),  1.0)
        endif
      enddo
c
      do n= 1,2
        dpmixl(:,:,  n) = dp(   :,:,1,n)
        th3d(  :,:,1,n) = thmix(:,:,  n)*1000.0 - thbase
        call rotate(th3d(1,1,1,n),  1.0)
      enddo
c
      do n= 1,3
        call rotate(ubavg(1,1,n),  0.01)  ! micom vbavg
        call rotate(vbavg(1,1,n), -0.01)  ! micom ubavg
        call rotate(pbavg(1,1,n),  0.1)
      enddo
c
      call rotate(pbot,  0.1)
      call rotate(psikk, 0.0001)
      thkk = thkk*1000.0 - thbase
      call rotate(thkk,  1.0)
c
      call restart_out(nstep, dtime)

      end program restart_micom
