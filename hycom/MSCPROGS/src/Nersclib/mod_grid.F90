module mod_grid
implicit none

! Module for reading and keeping grid details (depths,
! grid size, lon/lat etc.). One function only - get_grid

! Grid info 
real, dimension(:,:),allocatable, save :: plon,plat
real, dimension(:,:),allocatable, save :: ulon,ulat
real, dimension(:,:),allocatable, save :: vlon,vlat
real, dimension(:,:),allocatable, save :: qlon,qlat
real, dimension(:,:),allocatable, save :: scuy,scvx
real, dimension(:,:),allocatable, save :: scpy,scpx
real, dimension(:,:),allocatable, save :: depths
real, dimension(:,:),allocatable, save :: depthu, depthv
integer, dimension(:,:),allocatable, save :: ip,iu,iv

logical, save :: periodic

contains

subroutine get_grid()
   use mod_xc
   use mod_za
   use mod_parameters
   implicit none

   character(len=80)  gridid
   character(len=50)  nfile
   character(len=11)   tag7

   logical ex,ex2,exuf
   integer :: i,j, ia
   real*8, allocatable :: iodepths(:,:)
   integer :: nop=67
   logical :: exa, exb, exa2, exb2
   real, allocatable :: tmpio(:,:)
   integer, allocatable :: mask(:,:)
   real :: xmax, xmin
   real, external :: spherdist


   ! First, try getting from regional files
   inquire(file='regional.grid.a',exist=exa)
   inquire(file='regional.grid.b',exist=exb)
   inquire(file='regional.depth.a',exist=exa2)
   inquire(file='regional.depth.b',exist=exb2)

   if (exa.and.exb.and.exa2.and.exb2) then

     write(*,*)'Load grid positions from file: regional.grid.a'
     allocate(plon  (idm,jdm))
     allocate(plat  (idm,jdm))
     allocate(qlon  (idm,jdm))
     allocate(qlat  (idm,jdm))
     allocate(ulon  (idm,jdm))
     allocate(ulat  (idm,jdm))
     allocate(vlon  (idm,jdm))
     allocate(vlat  (idm,jdm))
     allocate(scpx  (idm,jdm))
     allocate(scpy  (idm,jdm))
     allocate(depths(idm,jdm))
     allocate(mask(idm,jdm)) 
     call zaiopf('regional.grid.a','old',nop)
     call zaiord(plon,mask,.false.,xmin,xmax,nop)
     call zaiord(plat,mask,.false.,xmin,xmax,nop)
     call zaiord(qlon,mask,.false.,xmin,xmax,nop)
     call zaiord(qlat,mask,.false.,xmin,xmax,nop)
     call zaiord(ulon,mask,.false.,xmin,xmax,nop)
     call zaiord(ulat,mask,.false.,xmin,xmax,nop)
     call zaiord(vlon,mask,.false.,xmin,xmax,nop)
     call zaiord(vlat,mask,.false.,xmin,xmax,nop)
     call zaiord(vlat,mask,.false.,xmin,xmax,nop) ! Not a bug, just skipped
     call zaiord(scpx,mask,.false.,xmin,xmax,nop)
     call zaiord(scpy,mask,.false.,xmin,xmax,nop)
     call zaiocl(nop)

     write(*,*)'Load depths from file: regional.depth.a'
     !print *,idm,jdm
     call zaiopf('regional.depth.a','old',nop)
     call zaiord(depths,mask,.false.,xmin,xmax,nop)
     call zaiocl(nop)

     where (depths > 0.5*huge) depths=0.

     write(*,*)'finished Load depths from file: regional.depth.a'

   ! Try getting from latlon.dat and depths file. Old, deprecated
   else 
      print *,'One of your regional.grid.[a,b] or regional.depth.[a,b] is missing, we try to read lonlat.dat'
      ! read lat-lon positions from grid
      inquire(file='./latlon.dat',exist=ex2)
      inquire(file='./Data/latlon.dat',exist=ex)
      if (ex) then
         nfile='./Data/latlon.dat'
        write(*,*)'Load grid positions from file: ./Data/latlon.dat'
      else if (ex2) then
         nfile='./latlon.dat'
        write(*,*)'Load grid positions from file: ./latlon.dat'
      else
         print *,'./Data/latlon.dat or ./latlon.dat does not exist'
         call exit(1)
      endif
      open(10,file=trim(nfile),form='formatted')
      READ(10,'(2i5)',ERR=250) idm,jdm
      250 continue

      allocate(qlat(idm,jdm))
      allocate(qlon(idm,jdm))
      allocate(ulat(idm,jdm))
      allocate(ulon(idm,jdm))
      allocate(vlat(idm,jdm))
      allocate(vlon(idm,jdm))
      allocate(plat(idm,jdm))
      allocate(plon(idm,jdm))
      allocate(scpx(idm,jdm))
      allocate(scpy(idm,jdm))
      allocate(tmpio(0:idm+1,0:jdm+1))


      read(10,'(a80)')gridid
      if (gridid(1:14) == 'CONFORMAL GRID') then
         print *,'Reading CONFORMAL GRID'
      else
         print *,'NBNBNB: This is not a conformal grid'
         print *,'Do you want to continue?'
         read (*,*)
      end if

      read(10,'(15e14.7)')((tmpio(i,j),i=0,idm+1),j=0,jdm+1)
      qlat=tmpio(1:idm,1:jdm)
      read(10,'(15e14.7)')((tmpio(i,j),i=0,idm+1),j=0,jdm+1)
      qlon=tmpio(1:idm,1:jdm)
      read(10,'(15e14.7)')((tmpio(i,j),i=0,idm+1),j=0,jdm+1)
      plat=tmpio(1:idm,1:jdm)
      read(10,'(15e14.7)')((tmpio(i,j),i=0,idm+1),j=0,jdm+1)
      plon=tmpio(1:idm,1:jdm)
      read(10,'(15e14.7)')((tmpio(i,j),i=0,idm+1),j=0,jdm+1)
      ulat=tmpio(1:idm,1:jdm)
      read(10,'(15e14.7)')((tmpio(i,j),i=0,idm+1),j=0,jdm+1)
      ulon=tmpio(1:idm,1:jdm)
      read(10,'(15e14.7)')((tmpio(i,j),i=0,idm+1),j=0,jdm+1)
      vlat=tmpio(1:idm,1:jdm)
      read(10,'(15e14.7)')((tmpio(i,j),i=0,idm+1),j=0,jdm+1)
      vlon=tmpio(1:idm,1:jdm)
      close(10)

      scpx(1:idm,1:jdm) = spherdist(ulon(1:idm  ,1:jdm  ),ulat(1:idm  ,1:jdm  ), &
                                    ulon(2:idm+1,1:jdm  ),ulat(2:idm+1,1:jdm  )   )
      scpy(1:idm,1:jdm) = spherdist(vlon(1:idm  ,1:jdm  ),vlat(1:idm  ,1:jdm  ), &
                                    vlon(1:idm  ,2:jdm+1),vlat(1:idm  ,2:jdm+1)   )

      deallocate (tmpio)


      ! Get depth mask -- from depth file 
      allocate(iodepths(idm,jdm))
      allocate(depths  (idm,jdm))
      if (idm<1000 .and. jdm<1000) then
         write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
      else
         write(tag7,'(i5.5,a,i5.5)')idm,'x',jdm
      end if
      inquire(file='depths'//trim(tag7)//'.uf',exist=ex)
      if (.not.ex) then
         print *,'depths file does not exist'
         call exit(1)
      end if 
      open (unit=10,file='depths'//trim(tag7)//'.uf',status='old',form='unformatted')
      read(10)iodepths
      close(10)
      !print *,maxval(iodepths),minval(iodepths)
      depths=iodepths
      deallocate(iodepths)

   end if



   ! Grid sizes (can also be retrieved from regional.grid))
   allocate(scuy(idm,jdm))
   allocate(scvx(idm,jdm))
   scuy(:,1:jdm-1) = spherdist(qlon(1:idm,1:jdm-1),qlat(1:idm,1:jdm-1), &
                               qlon(1:idm,2:jdm  ),qlat(1:idm,2:jdm  )   )
   scuy(:,jdm)=scuy(:,jdm-1)  ! Approx - should use periodicity if available

   scvx(1:idm-1,:) = spherdist(qlon(1:idm-1,1:jdm),qlat(1:idm-1,1:jdm  ), &
                               qlon(2:idm  ,1:jdm),qlat(2:idm  ,1:jdm))
   scvx(idm,:)=scvx(idm-1,:)! Approx

   ! Is grid periodic in i?
   periodic = any(depths(1,:)<.1) .and. any(depths(idm,:)<.1)

   ! Allocate bigrid stuff
   allocate(ip(idm,jdm))
   ip=0
   where (depths>.1 .and. depths < huge) ip=1

   ! Depths in u/v points
   ! TODO: Correct for periodic grid
   allocate(depthu(idm,jdm))
   allocate(depthv(idm,jdm))
   allocate(iu(idm,jdm))
   allocate(iv(idm,jdm))
   depthu=0.
   depthv=0.
   iu=0
   iv=0
   do j=2,jdm
   do i=1,idm
      if (periodic) then
         ia=mod(idm+i-2,idm)+1
      else
         ia=max(i-1,1)
      end if
      depthu(i,j)=min(depths(i,j),depths(ia ,j  ))
      depthv(i,j)=min(depths(i,j),depths(i  ,j-1))
      iu(i,j) = min(ip(i,j),ip(ia,j))
      iv(i,j) = min(ip(i,j),ip(i,j-1))
   end do
   end do


   write(*,*)'finished Loading grid/depth from file: regional.depth/grid.a'


end subroutine get_grid

end module mod_grid

