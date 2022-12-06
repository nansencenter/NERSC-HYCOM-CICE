! Program to create the files thkdf4.[a,b] read by HYCOM when the value is spatially variable
! Here, the thickness diffusion-4 is a function of latitude:
! thkdf4(i,j)=0.6*cos(abs(modlat(i,j)-45))+0.05
! we found some ugly numerical noise below the mixed layer in the Gulf Stream, which is well corrected when thkdf4
! is added. However, thkdf4 hampers the advection of atlantic water into the Nordic Seas and the 
! Arctic. So the compromise here is to increase the parameter value around the Gulf Stream ~45N only. 
!Created by francois.counillon@nersc.no, 2009, rephrased by laurent.bertino@nersc.no 

program thkdf4
   use mod_xc
   use mod_za
   use mod_grid, only : get_grid, depths, plon, plat, qlat, qlon
   implicit none

   character(len=24) :: filename
   integer, allocatable, dimension(:,:) :: ip
   real, allocatable, dimension(:,:) :: fild
   integer :: ii,jj
   real ::  hmin, hmax
   real pi
   data pi/3.14159265/

!   ! Initialize IO for .ab files
   CALL XCSPMD()  
   allocate(ip(idm,jdm))
   allocate(fild(idm,jdm))
   CALL ZAIOST()
   call get_grid()
   do ii=1,idm
      do jj=1,jdm
         !fild(ii,jj)=max(0.004,max(cos((plat(ii,jj)-35)*pi/90.),0)**4*0.06)+max(0.001,max(cos((plat(ii,jj)-5)*pi/90.),0)**4*0.06)
         !increased thkdf4 off brazil for Met.no system
         fild(ii,jj)=max(0.004,max(cos((plat(ii,jj)-35)*pi/90.),0.)**4*0.06)+max(0.001,max(cos((plat(ii,jj)-5)*pi/90.),0.)**4*0.15)
      end do
   end do
   print *, 'field computed',fild(140,150)
    call zaiopf('thkdf4.a', 'replace', 909)
    call zaiowr(fild(:,:),ip,.false.,hmin,hmax,909,.true.)
    open (unit=909, file='thkdf4.b',  &
         status='replace', action='write')
    write(909,'(" thkdf4:range = ",2e16.8)') hmin,hmax
    close(909)
    call zaiocl(909)
end program thkdf4
