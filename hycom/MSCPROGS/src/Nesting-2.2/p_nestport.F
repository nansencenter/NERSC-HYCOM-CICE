      program nestport
      use mod_xc
      use mod_za
      implicit none
      real   , dimension(:,:), allocatable :: depths
      integer, dimension(:,:), allocatable :: ip,iu,iv
      character(len=80) :: cline
      real :: hmin, hmax, hminb, hmaxb
      integer :: i,j,ia,ja,nports,k
      integer, parameter :: maxnport=100
      integer,dimension(maxnport) :: kdport
      integer,dimension(maxnport) :: ifport, ilport
      integer,dimension(maxnport) :: jfport, jlport
      logical :: inleta, inletb, inlet
c --- Init io and arrays
      call xcspmd()
      call zaiost()
      allocate(depths(idm,jdm))
      allocate(ip(idm,jdm))
      allocate(iu(idm,jdm))
      allocate(iv(idm,jdm))
c
c --- Open and read bathymetry
      call zaiopf('regional.depth.a', 'old', 701)
      open (unit=701,file='regional.depth.b',status='old')
      call zaiord(depths,ip,.false.,hmin,hmax,701)
      read(701,'(a80)') cline ; print *,cline
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      i = index(cline,'=')
      read (cline(i+1:),*) hminb,hmaxb
      if (abs(hmin-hminb).gt.abs(hminb)*1.e-4 .or. 
     &    abs(hmax-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') 
     &      'error - .a and .b files not consistent:', 
     &      '.a,.b min = ',hmin,hminb,hmin-hminb, 
     &      '.a,.b max = ',hmax,hmaxb,hmax-hmaxb
          print *, '(forfun_nersc:plon)'
         call exit(1)
       endif
      ! Remove 1e30 points
      call zaiocl(701)
      close(701)
      where (depths>1e20) depths=0.
c
c --- Set masks
      do j=1,jdm
      do i=1,idm
         ia=mod(idm-2+i,idm)+1 ! periodic i-1
         ja=max(j-1,1)
         ip(i,j)=0.
         iu(i,j)=0.
         iv(i,j)=0.
         if (depths(i,j)>0.1) ip(i,j)=1
         if (depths(i,j)>0.1 .and. depths(ia,j)>0.1) iu(i,j)=1
         if (depths(i,j)>0.1 .and. depths(i,ja)>0.1) iv(i,j)=1
      end do
      end do
c
c --- Go along boundary and specify port input .. i=idm-1:
      nports=0
      i=idm
      do j=1,jdm
         if (ip(i,j).ne.0) then
            print *,'boundary must be land'
            call exit(1)
         end if
         inleta=iu(i-1,j-1)==1 .and. iu(i-2,j-1)==1
         inlet= iu(i-1,j  )==1 .and. iu(i-2,j  )==1
         inletb=iu(i-1,j+1)==1 .and. iu(i-2,j+1)==1
         if      (j==1.and.inlet) then
            nports=nports+1
            kdport(nports)=3
            jfport(nports)=j
            ifport(nports)=i
         else if (j==jdm .and. inlet) then
            jlport(nports)=j
            ilport(nports)=i
         else if (inlet .and. .not.inleta) then
            nports=nports+1
            kdport(nports)=3
            jfport(nports)=j
            ifport(nports)=i
         else if (inlet .and. .not.inletb) then
            jlport(nports)=j
            ilport(nports)=i
         end if
      end do
      print *
c --- Go along boundary and specify port input .. i=2:
      i=1
      do j=1,jdm
         if (ip(i,j).ne.0) then
            print *,'boundary must be land'
            call exit(1)
         end if
         inleta=iu(i+2,j-1)==1 .and. iu(i+3,j-1)==1
         inlet= iu(i+2,j  )==1 .and. iu(i+3,j  )==1
         inletb=iu(i+2,j+1)==1 .and. iu(i+3,j+1)==1
         if      (j==1.and.inlet) then
            nports=nports+1
            kdport(nports)=4
            jfport(nports)=j
            ifport(nports)=i+1
         else if (j==jdm.and.inlet) then
            jlport(nports)=j
            ilport(nports)=i+1
         else if (inlet.and. .not.inleta) then
            nports=nports+1
            kdport(nports)=4
            jfport(nports)=j
            ifport(nports)=i+1
         else if (inlet.and. .not.inletb) then
            jlport(nports)=j
            ilport(nports)=i+1
         end if
      end do
      print *
c --- Go along boundary and specify port input .. j=jdm-1:
      j=jdm
      do i=1,idm
         if (ip(i,j).ne.0) then
            print *,'boundary must be land'
            call exit(1)
         end if
         inleta=iv(i-1,j-1)==1 .and. iv(i-1,j-2)==1
         inlet= iv(i  ,j-1)==1 .and. iv(i  ,j-2)==1
         inletb=iv(i+1,j-1)==1 .and. iv(i+1,j-2)==1
         if      (i==1.and.inlet) then
            nports=nports+1
            kdport(nports)=1
            jfport(nports)=j
            ifport(nports)=i
         else if (i==idm.and.inlet) then
            jlport(nports)=j
            ilport(nports)=i
         else if (inlet .and. .not.inleta) then
            nports=nports+1
            kdport(nports)=1
            jfport(nports)=j
            ifport(nports)=i
         else if (inlet .and. .not. inletb) then
            jlport(nports)=j
            ilport(nports)=i
         end if
      end do
      print *
c --- Go along boundary and specify port input .. j=2:
      j=1
      do i=1,idm
         if (ip(i,j).ne.0) then
            print *,'boundary must be land'
            call exit(1)
         end if
         inleta=iv(i-1,j+2)==1 .and. iv(i-1,j+3)==1
         inlet= iv(i  ,j+2)==1 .and. iv(i  ,j+3)==1
         inletb=iv(i+1,j+2)==1 .and. iv(i+1,j+3)==1
         if      (i==1.and.inlet) then
            nports=nports+1
            kdport(nports)=2
            jfport(nports)=j+1
            ifport(nports)=i
         else if (i==idm.and.inlet) then
            jlport(nports)=j+1
            ilport(nports)=i
         else if (inlet.and. .not. inleta) then
            nports=nports+1
            kdport(nports)=2
            jfport(nports)=j+1
            ifport(nports)=i
         else if (inlet.and. .not. inletb) then
            jlport(nports)=j+1
            ilport(nports)=i
         end if
      end do
c 
c --- Dump the result to ports.input
      open(10,file='ports.input.nest',status='replace',form='formatted')
      write(10,'(i6,4x,a)') nports,
     &   "'nports' = number of boundary port sections"
      do k=1,nports
         write(10,'(i6,4x,a)') kdport(k),
     &      "'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)"
         write(10,'(i6,4x,a)') ifport(k),
     &      "'ifport' = first i-index"
         write(10,'(i6,4x,a)') ilport(k),
     &      "'ilport' = last  i-index"
         write(10,'(i6,4x,a)') jfport(k),
     &      "'jfport' = first j-index"
         write(10,'(i6,4x,a)') jlport(k),
     &      "'jlport' = last  j-index"
      end do
      close(10)
      print *,'nestport finished'
      print *,'number of ports: nports=',nports
      print *,'dumping to ports.input.nest'
      end program



















