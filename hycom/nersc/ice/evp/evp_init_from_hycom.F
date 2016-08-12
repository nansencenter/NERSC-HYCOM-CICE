      subroutine evp_init_from_hycom()
      use mod_xc
      use mod_evp , only : tmask, umask, dxt, dxt2, dxt4,
     &    dyt, dyt2, dyt4, tarea, tarear,tinyarea,puny,
     &    dxu, dyu, uarea, uarear, hte, htn, dxhy, dyhx,
     &    cxp, cyp, cxm, cym, dt, e_itst, e_jtst,ndte, 
     &    Pstar0,Cstar0,
     &    evp_ULAT=>ULAT
      use mod_raw_io
      implicit none
      include 'common_blocks.h' ! HYCOM variables
      integer :: i,j
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp
      real :: facrad
      real :: fver, cver
      logical :: ex, lfatal
      real :: c0,p0
      integer :: lmargin
      real*4, dimension(itdm,jtdm) :: globscuy, globscvx,globqlat
      real*4 :: amn,amx,spval
      integer :: i2,j2
      facrad=asin(1.)/90.



      !Create EVP t-cell mask from hycom p-mask
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         tmask(i,j) = ip(i,j)==1
      end do
      end do     
 
      !Create EVP u-cell mask from EVP t-mask
      umask=.false.
      tmp=0.
      do j=1-nbdy+1,jj+nbdy-1
      do i=1-nbdy+1,ii+nbdy-1
         umask(i,j) = 
     .      tmask(i+1,j  ) .and. 
     .      tmask(i+1,j+1) .and. 
     .      tmask(i  ,j  ) .and. 
     .      tmask(i  ,j+1) 
         if (umask(i,j)) tmp(i,j)=1.
      end do
      end do     
      call xctilr(tmp   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)

      umask=.false.
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         umask(i,j) = tmp(i,j)>.5
      end do
      end do
 
      !Create EVP t-cell dimensions from HYCOM p-cell
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         dxt (i,j)=scpx(i,j)
         dxt2(i,j)=scpx(i,j)*0.5
         dxt4(i,j)=scpx(i,j)*0.25

         dyt (i,j)=scpy(i,j)
         dyt2(i,j)=scpy(i,j)*0.5
         dyt4(i,j)=scpy(i,j)*0.25

         tarea(i,j)=dxt(i,j)*dyt(i,j)
         tarear(i,j)= 1./ tarea(i,j)
         tinyarea(i,j) = puny*tarea(i,j)

      end do
      end do     
      call xctilr(dxt     ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(dxt2    ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(dxt4    ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(dyt     ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(dyt2    ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(dyt4    ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(tarea   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(tarear  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(tinyarea( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)


      

      

 
      !Create EVP u-cell dimensions from HYCOM p-cell
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         dxu (i,j)=scqx(i,j)
         dyu (i,j)=scqy(i,j)
         uarea(i,j)=dxu(i,j)*dyu(i,j)
         uarear(i,j)= 1./ uarea(i,j)
      end do
      end do     
      call xctilr(dxu   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(dyu   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(uarea ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(uarear( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)

 
      !Create EVP eastern t-cell edge length
      call READRAW(globqlat,AMN,AMX,ITDM,JTDM,.false.,SPVAL,
     &   'regional.grid.a',4)
      call READRAW(globscuy,AMN,AMX,ITDM,JTDM,.false.,SPVAL,
     &   'regional.grid.a',15)
      call READRAW(globscvx,AMN,AMX,ITDM,JTDM,.false.,SPVAL,
     &   'regional.grid.a',17)

      lmargin=nbdy
      do j=1-lmargin,jj+lmargin
      do i=1-lmargin,ii+lmargin
         !HTE(i,j)=scuy(i+1,j)
         j2=j0+j-1
         i2=i0+i
         ! Temporary solution - works on closed domains
         i2=max(min(itdm,i2),1)
         j2=max(min(jtdm,j2),1)
         HTE(i,j)=globscuy(i2,j2)

      end do
      end do     

      !Create EVP northern t-cell edge length
      lmargin=nbdy
      do j=1-lmargin,jj+lmargin
      do i=1-lmargin,ii+lmargin
         !HTN(i,j)=scvx(i,j+1)
         j2=j0+j
         i2=i0+i-1

         ! Temporary solution - works on closed domains
         i2=max(min(itdm,i2),1)
         j2=max(min(jtdm,j2),1)
         HTN(i,j)=globscvx(i2,j2)
      end do
      end do     

      !Create EVP eastern t-cell edge length grid increment
      do j=1-nbdy+1,jj+nbdy-1
      do i=1-nbdy+1,ii+nbdy-1
         dxhy(i,j) = 0.5*(HTE(i,j) - HTE(i-1,j))
      end do
      end do     

      !Create EVP northern t-cell edge length grid increment
      do j=1-nbdy+1,jj+nbdy-1
      do i=1-nbdy+1,ii+nbdy-1
         dyhx(i,j) = 0.5*(HTN(i,j) - HTN(i,j-1))
      end do
      end do     

      !Create EVP u-cell latitudes - TODO - make this more accurate
      do j=1-nbdy+1,jj+nbdy-1
      do i=1-nbdy+1,ii+nbdy-1
C
C        evp_ULAT(i,j) = 
C    &      (plat(i  ,j+1) + plat(i+1,j  ) +
C    &       plat(i  ,j  ) + plat(i+1,j+1))*.25
C
C        ! To radians
C        evp_ULAT(i,j) = evp_ULAT(i,j) * facrad
C
         j2=j0+j
         i2=i0+i

         ! Temporary solution - works on closed domains
         i2=max(min(itdm,i2),1)
         j2=max(min(jtdm,j2),1)

         evp_ULAT(i,j) = globqlat(i2,j2) * facrad

      end do
      end do     

      do j=1-nbdy+2,jj+nbdy-2
      do i=1-nbdy+2,ii+nbdy-2
        cyp(i,j) = (1.5*HTE(i,j) - 0.5*HTE(i-1,j))
        cxp(i,j) = (1.5*HTN(i,j) - 0.5*HTN(i,j-1))
        cym(i,j) = (0.5*HTE(i,j) - 1.5*HTE(i-1,j))
        cxm(i,j) = (0.5*HTN(i,j) - 1.5*HTN(i,j-1))
      end do
      end do

      call xctilr(HTE   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_us)
      call xctilr(HTN   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_vs)
      call xctilr(dxhy  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_us)
      call xctilr(dyhx  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_us)
      call xctilr(evp_ULAT  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(cyp   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(cxp   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(cym   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(cxm   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)


      e_itst=itest
      e_jtst=jtest
      if (itest>-99 .and. jtest >-99) then
         print *,'evp_init_hyc',dxt(itest,jtest)
         print *,'evp_init_hyc',dyt(itest,jtest)
         print *,'evp_init_hyc',dxu(itest,jtest)
         print *,'evp_init_hyc',dyu(itest,jtest)
         print *,'evp_init_hyc',HTE(itest,jtest)
         print *,'evp_init_hyc',HTN(itest,jtest)
         print *,'evp_init_hyc',evp_ULAT(itest,jtest)
         print *,'evp_init_hyc',umask(itest,jtest)
         print *,'evp_init_hyc',tmask(itest,jtest)
      end if

      ! Time step for now
      ! Read init file 
      cver=1.1 ! Current version of infile
      dt=3600.
      p0=27500.
      c0=20.0
      ndte=120.
      lfatal=.false.
      inquire(exist=ex, file='infile.evp')
      if (ex) then
         open(10,file='infile.evp',action='read')
         read(10,*) fver ! Read infile version
         if (abs(fver-cver)<0.0001) then
            read(10,*) dt ! Read infile dynamics time step
            read(10,*) p0 ! Read infile plastic strength
            read(10,*) c0 ! Read infile concentration factor
            read(10,*) ndte ! Read infile number of subcycles 
            close(10)
            Pstar0=p0
            Cstar0=c0
         else
            if (mnproc==1) print *, 'Wrong version of infile.evp'
            lfatal=.true.
         end if
      else
         if (mnproc==1) print *, 'infile.evp not found'
         lfatal=.true.
      end if

      if (lfatal) then
         if (mnproc==1) then
            print *,'IO error, an example infile.evp is given below:'
            print *,' -----------------   snip   -------------------'
            write(*,'(f3.1,a)') cver,'      # File version'
            write(*,'(f6.0,a)')   dt,'   # EVP time step'
            write(*,'(f6.0,a)') p0,'   # EVP ice strength'
            write(*,'(f4.0,a)') c0,'     # EVP ice concentration factor'
            write(*,*) ndte,'     # EVP number of subcycles'
            print *,' -----------------   snip   -------------------'
         end if
         call xcstop('(evp_init_from_hycom)')
         stop '(evp_init_from_hycom)'
      end if


      end subroutine evp_init_from_hycom
