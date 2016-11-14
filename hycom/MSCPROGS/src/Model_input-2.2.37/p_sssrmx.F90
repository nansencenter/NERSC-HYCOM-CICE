! Program to create the files sssrmx.[a,b] read by HYCOM when the value is not constant
! Here, the value follows the following equation:
! sssrmx(:,:)=0.5
!Created by Fanf

program sssrmx
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
    fild(:,:)=0.5
    call zaiopf('sssrmx.a', 'replace', 909)
    call zaiowr(fild(:,:),ip,.false.,hmin,hmax,909,.true.)
    open (unit=909, file='sssrmx.b',  &
         status='replace', action='write')
         write(909,'(a)') 'Threshold above which we do not relax'
    write(909,'(a)') ''
    write(909,'(a)') ''
    write(909,'(a)') ''
    write(909,'(a,2i5)') 'i/jdm = ',idm,jdm
    write(909,'(" sssrmx:range = ",2e16.8)') hmin,hmax
    close(909)
    call zaiocl(909)
end program sssrmx
